#!/bin/bash
# Full single-capture pipeline: download -> cellranger -> chromosome-shard velocyto -> merge -> archive
# Usage: process_capture.sh <capture_name>
# Expects capture_srr_map.tsv in the same directory as this script.
set -uo pipefail

CAPTURE=$1
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
MAPFILE=$SCRIPT_DIR/capture_srr_map.tsv

SRATOOLS=/home/huc/tools/sratoolkit.3.2.1-ubuntu64/bin
CELLRANGER=/home/huc/bin/cellranger
REF=/home/huc/projects/ref_genome/Homo_sapiens/refdata-gex-GRCh38-2024-A
GTF=$REF/genes/genes.gtf
FASTA_FAI=$REF/fasta/genome.fa.fai
PERSIST=/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/raw_repro_velocity
NCPU=${SLURM_CPUS_PER_TASK:-32}
# derive a safe cellranger memory ceiling (GB) from what SLURM actually granted,
# leaving headroom so cellranger doesn't get cgroup-OOM-killed
if [ -n "${SLURM_MEM_PER_NODE:-}" ]; then
  MEMGB=$(( SLURM_MEM_PER_NODE / 1024 - 2 ))
elif [ -n "${SLURM_MEM_PER_CPU:-}" ]; then
  MEMGB=$(( SLURM_MEM_PER_CPU * NCPU / 1024 - 2 ))
else
  MEMGB=60
fi
[ $MEMGB -lt 4 ] && MEMGB=4

WORK=$PERSIST/slurm_scratch/$CAPTURE
mkdir -p $WORK/{sra,fastq_raw,fastq,cellranger_out,shard_bam,shard_gtf,shard_loom,logs}
cd $WORK

echo "=== [$(date)] $CAPTURE : running on $(hostname), $NCPU cpus ==="

if [ -f "$PERSIST/$CAPTURE/DONE" ]; then
  echo "$CAPTURE already completed, exiting."
  exit 0
fi

RUNS=($(awk -F'\t' -v c="$CAPTURE" '$1==c{split($4,a,";"); for(i in a) print a[i]}' $MAPFILE))
echo "Runs: ${RUNS[@]}"

set -e

EXPECTED_FASTQS=$(( ${#RUNS[@]} * 2 ))
FOUND_FASTQS=$(find fastq -maxdepth 1 -iname "${CAPTURE}_S1_L*_R[12]_001.fastq.gz" 2>/dev/null | wc -l)
if [ "$FOUND_FASTQS" -eq "$EXPECTED_FASTQS" ]; then
  echo "--- [$(date)] STAGE 1-3: fastq already staged ($FOUND_FASTQS files), skipping download/dump/rename ---"
else
  echo "--- [$(date)] STAGE 1: prefetch ---"
  for srr in "${RUNS[@]}"; do
    $SRATOOLS/prefetch --max-size 100G -O sra "$srr"
  done

  echo "--- [$(date)] STAGE 2: fasterq-dump ---"
  for srr in "${RUNS[@]}"; do
    $SRATOOLS/fasterq-dump --split-files --include-technical -e $NCPU -t $WORK/tmp_fqd -O fastq_raw/$srr sra/$srr/$srr.sra
  done

  echo "--- [$(date)] STAGE 3: rename to cellranger convention ---"
  lane=0
  for srr in "${RUNS[@]}"; do
    lane=$((lane+1))
    Lstr=$(printf "L%03d" $lane)
    gzip -c fastq_raw/$srr/${srr}_3.fastq > fastq/${CAPTURE}_S1_${Lstr}_R1_001.fastq.gz &
    gzip -c fastq_raw/$srr/${srr}_4.fastq > fastq/${CAPTURE}_S1_${Lstr}_R2_001.fastq.gz &
    wait
    rm -rf fastq_raw/$srr sra/$srr
  done
fi

echo "--- [$(date)] STAGE 4: cellranger count ---"
cd cellranger_out
if [ -d "$CAPTURE" ] && [ -f "$CAPTURE/_invocation" ]; then
  echo "existing pipestance for $CAPTURE found -- resuming (not deleting), clearing stale lock if present"
  rm -f "$CAPTURE/_lock"
else
  rm -rf "$CAPTURE"
fi
$CELLRANGER count --id=$CAPTURE \
  --transcriptome=$REF \
  --fastqs=$WORK/fastq \
  --sample=$CAPTURE \
  --create-bam=true \
  --localcores=$NCPU --localmem=$MEMGB
cd $WORK
rm -rf fastq

BAM=cellranger_out/$CAPTURE/outs/possorted_genome_bam.bam
BC=cellranger_out/$CAPTURE/outs/filtered_feature_bc_matrix/barcodes.tsv.gz

echo "--- [$(date)] STAGE 5: chromosome-shard BAM + GTF ---"
SAMTOOLS=/home/huc/bin/samtools
CHROMS=$(cut -f1 $FASTA_FAI | grep -E "^chr([0-9]+|[XYM])$")
for chr in $CHROMS; do
  ( $SAMTOOLS view -@2 -b $BAM $chr > shard_bam/${chr}.bam
    $SAMTOOLS index shard_bam/${chr}.bam
    awk -v c=$chr '$1==c' $GTF > shard_gtf/${chr}.gtf ) &
  # throttle to avoid oversubscribing
  while [ "$(jobs -r | wc -l)" -ge $((NCPU/2)) ]; do wait -n; done
done
wait
# everything not in the primary chromosome set (alt contigs/scaffolds) -> one "other" shard
PRIMARY_REGION_FILE=shard_bam/primary_chroms.txt
echo "$CHROMS" | tr ' ' '\n' > $PRIMARY_REGION_FILE
$SAMTOOLS view -H $BAM | grep "^@SQ" | sed -E 's/.*SN:([^\t]+).*/\1/' | grep -vFf $PRIMARY_REGION_FILE > shard_bam/other_chroms.txt || true
if [ -s shard_bam/other_chroms.txt ]; then
  $SAMTOOLS view -@$NCPU -b $BAM $(cat shard_bam/other_chroms.txt | tr '\n' ' ') > shard_bam/other.bam
  $SAMTOOLS index shard_bam/other.bam
  grep -Ff shard_bam/other_chroms.txt $GTF > shard_gtf/other.gtf || touch shard_gtf/other.gtf
  CHROMS="$CHROMS other"
fi

echo "--- [$(date)] STAGE 6: velocyto per shard (parallel) ---"
source /home/huc/miniconda3/etc/profile.d/conda.sh
conda activate velocyto_py39
for chr in $CHROMS; do
  if [ ! -s shard_gtf/${chr}.gtf ]; then continue; fi
  ( velocyto run -b $BC -o shard_loom/$chr -@2 shard_bam/${chr}.bam shard_gtf/${chr}.gtf \
      > logs/velocyto_${chr}.log 2>&1 ) &
  while [ "$(jobs -r | wc -l)" -ge $((NCPU/2)) ]; do wait -n; done
done
wait

echo "--- [$(date)] STAGE 7: merge shard looms ---"
SHARD_LOOMS=$(find shard_loom -iname "*.loom")
python3 $SCRIPT_DIR/merge_loom_shards.py $CAPTURE $BC shard_loom/${CAPTURE}_merged.loom $SHARD_LOOMS
conda deactivate

echo "--- [$(date)] STAGE 8: archive to persistent storage ---"
mkdir -p $PERSIST/$CAPTURE/outs
cp -r cellranger_out/$CAPTURE/outs/filtered_feature_bc_matrix $PERSIST/$CAPTURE/outs/
cp cellranger_out/$CAPTURE/outs/filtered_feature_bc_matrix.h5 $PERSIST/$CAPTURE/outs/ 2>/dev/null || true
cp cellranger_out/$CAPTURE/outs/metrics_summary.csv $PERSIST/$CAPTURE/outs/
cp cellranger_out/$CAPTURE/outs/molecule_info.h5 $PERSIST/$CAPTURE/outs/ 2>/dev/null || true
cp shard_loom/${CAPTURE}_merged.loom $PERSIST/$CAPTURE/${CAPTURE}.loom
touch $PERSIST/$CAPTURE/DONE

echo "--- [$(date)] STAGE 9: cleanup scratch ---"
cd $PERSIST
rm -rf $WORK

echo "=== [$(date)] $CAPTURE COMPLETE ==="
