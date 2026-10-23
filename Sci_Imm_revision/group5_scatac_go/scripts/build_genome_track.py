import anndata as ad
import pysam
import numpy as np
import pandas as pd
import gzip

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group5_scatac_go'
FRAG_DIR = '/home/huc/projects/Ruoqiong_Cao.IAHSC/scATAC/1.mapping/cellranger_output'
GTF = '/home/huc/projects/ref_genome/Mus_musculus/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz'

GENE = 'Irf1'
CHROM = 'chr11'
GENE_START, GENE_END = 53770014, 53778374
# window chosen to cover the gene body PLUS the actual SnapATAC2 diff_test peaks
# linked to Irf1 (peak2gene.Infected_vs_Naive.opening_in_Infected.csv), which sit
# 8-43kb downstream of the gene, NOT at the promoter
WIN_START, WIN_END = 53766000, 53825000
BIN = 50  # bp resolution

DIFF_PEAKS = [
    (53786249, 53786750, 1.589456, 0.004019, 'opening'),
    (53794808, 53795309, 0.842518, 0.019041, 'opening'),
    (53805711, 53806212, 1.399979, 0.013777, 'opening'),
    (53820664, 53821165, 0.924641, 0.013777, 'opening'),
]

# ---- 1. cell barcodes per group, per sample ----
a = ad.read_h5ad(f'{OUT.replace("group5_scatac_go","")}scATAC/4.analysis/peak_calling/IAHSC_by_infection_peak_matrix_0208.h5ad'.replace('/Sci_Imm_revision/',''), backed='r') if False else None
a = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scATAC/4.analysis/peak_calling/IAHSC_by_infection_peak_matrix_0208.h5ad', backed='r')
obs = a.obs[['sample','infection_status']].copy()
obs['barcode'] = [i.split(':')[1] for i in obs.index]
print(obs.groupby(['infection_status'])['sample'].nunique())
n_cells_per_group = obs['infection_status'].value_counts()
print(n_cells_per_group)

# ---- 2. build per-group coverage in the window ----
n_bins = (WIN_END - WIN_START) // BIN + 1
coverage = {'naive': np.zeros(n_bins), 'infected': np.zeros(n_bins)}

for sample, sub in obs.groupby('sample'):
    status = sub['infection_status'].iloc[0]
    bc_set = set(sub['barcode'])
    frag_path = f'{FRAG_DIR}/{sample}/outs/fragments.tsv.gz'
    tbx = pysam.TabixFile(frag_path)
    n_frags_used = 0
    try:
        for row in tbx.fetch(CHROM, max(0,WIN_START), WIN_END):
            fields = row.split('\t')
            fchrom, fstart, fend, fbc = fields[0], int(fields[1]), int(fields[2]), fields[3]
            if fbc not in bc_set:
                continue
            n_frags_used += 1
            # count both ends (Tn5 insertion sites) into bins
            for pos in (fstart, fend):
                if WIN_START <= pos < WIN_END:
                    b = (pos - WIN_START)//BIN
                    coverage[status][b] += 1
    finally:
        tbx.close()
    print(f'{sample} ({status}): {n_frags_used} fragment-ends counted from {len(bc_set)} cells')

# normalize per group by number of cells (counts per 1000 cells, CPK-cell-like)
for status in coverage:
    ncell = n_cells_per_group[status]
    coverage[status] = coverage[status] / ncell * 1000

np.savez(f'{OUT}/genome_track_{GENE}_coverage.npz',
         naive=coverage['naive'], infected=coverage['infected'],
         win_start=WIN_START, win_end=WIN_END, bin=BIN)
print("saved coverage arrays")

# ---- 3. gene model (all exons for this gene, merged across transcripts) ----
exons = []
strand = None
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if fields[2] != 'exon':
            continue
        if fields[0] != CHROM:
            continue
        start, end = int(fields[3]), int(fields[4])
        if end < WIN_START or start > WIN_END:
            continue
        if f'gene_name "{GENE}"' not in fields[8]:
            continue
        strand = fields[6]
        exons.append((start, end))
exons = sorted(set(exons))
print(f"{len(exons)} exon intervals found for {GENE}, strand={strand}")
pd.DataFrame(exons, columns=['start','end']).to_csv(f'{OUT}/genome_track_{GENE}_exons.csv', index=False)
with open(f'{OUT}/genome_track_{GENE}_meta.txt','w') as f:
    f.write(f'{CHROM}\t{WIN_START}\t{WIN_END}\t{strand}\t{GENE_START}\t{GENE_END}\n')
print("done")
