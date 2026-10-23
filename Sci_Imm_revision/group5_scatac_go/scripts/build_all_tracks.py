import anndata as ad
import pysam
import numpy as np
import pandas as pd
import gzip
import re
import json

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group5_scatac_go'
FRAG_DIR = '/home/huc/projects/Ruoqiong_Cao.IAHSC/scATAC/1.mapping/cellranger_output'
GTF = '/home/huc/projects/ref_genome/Mus_musculus/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz'

GENES = ['Irf1','Isg15','Igtp','Nlrc5','Psmb8','Psmb9','Tap1','Bst2','Slfn9',
         'Casp12','H2-Aa','H2-Ab1','H2-Eb1','Mvb12a','Pim1','Plac8','Serpina3f','Serpina3g']
GENES_U = [g.upper() for g in GENES]
BIN = 50
PAD_MIN = 3000

# ---- gene coords from GTF ----
gene_info = {}
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if fields[2] != 'gene':
            continue
        m = re.search(r'gene_name "([^"]+)"', fields[8])
        if not m:
            continue
        gname = m.group(1).upper()
        if gname in GENES_U:
            gene_info[gname] = dict(chrom=fields[0], start=int(fields[3]), end=int(fields[4]), strand=fields[6])

# ---- exons per gene ----
exon_info = {g: [] for g in GENES_U}
with gzip.open(GTF, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.rstrip('\n').split('\t')
        if fields[2] != 'exon':
            continue
        m = re.search(r'gene_name "([^"]+)"', fields[8])
        if not m:
            continue
        gname = m.group(1).upper()
        if gname in GENES_U:
            exon_info[gname].append((int(fields[3]), int(fields[4])))
for g in exon_info:
    exon_info[g] = sorted(set(exon_info[g]))

# ---- linked opening peaks + diff stats ----
p2g = pd.read_csv(f'{OUT}/peak2gene.Infected_vs_Naive.opening_in_Infected.csv')
p2g['gene_u'] = p2g['gene'].str.upper()
diff = pd.read_csv(f'{OUT}/IAHSC.Infected_vs_Naive.diff_peaks.csv')
diff.columns = ['peak', 'log2fc', 'pval', 'padj']

peaks_per_gene = {}
for g in GENES_U:
    sub = p2g[p2g['gene_u'] == g].merge(diff, on='peak')
    plist = []
    for _, row in sub.iterrows():
        chrom, coords = row['peak'].split(':')
        s, e = coords.split('-')
        plist.append(dict(chrom=chrom, start=int(s), end=int(e), log2fc=row['log2fc'], padj=row['padj']))
    peaks_per_gene[g] = plist

# ---- windows ----
windows = {}
for g in GENES_U:
    gi = gene_info[g]
    starts = [gi['start']] + [p['start'] for p in peaks_per_gene[g]]
    ends = [gi['end']] + [p['end'] for p in peaks_per_gene[g]]
    ws = min(starts) - PAD_MIN
    we = max(ends) + PAD_MIN
    windows[g] = (gi['chrom'], ws, we)
    print(g, gi['chrom'], ws, we, we - ws, 'bp')

# ---- cell groups ----
a = ad.read_h5ad(
    '/home/huc/projects/Ruoqiong_Cao.IAHSC/scATAC/4.analysis/peak_calling/IAHSC_by_infection_peak_matrix_0208.h5ad',
    backed='r')
obs = a.obs[['sample', 'infection_status']].copy()
obs['barcode'] = [i.split(':')[1] for i in obs.index]
n_cells_per_group = obs['infection_status'].value_counts()

samples = sorted(obs['sample'].unique())
sample_bc = {}
sample_status = {}
for sample, sub in obs.groupby('sample', observed=True):
    sample_bc[sample] = set(sub['barcode'])
    sample_status[sample] = sub['infection_status'].iloc[0]

tbx_handles = {s: pysam.TabixFile(f'{FRAG_DIR}/{s}/outs/fragments.tsv.gz') for s in samples}

# ---- build coverage per gene ----
results = {}
for g in GENES_U:
    chrom, ws, we = windows[g]
    n_bins = (we - ws) // BIN + 1
    coverage = {'naive': np.zeros(n_bins), 'infected': np.zeros(n_bins)}
    for sample in samples:
        status = sample_status[sample]
        bc_set = sample_bc[sample]
        tbx = tbx_handles[sample]
        try:
            for row in tbx.fetch(chrom, max(0, ws), we):
                fields = row.split('\t')
                fstart, fend, fbc = int(fields[1]), int(fields[2]), fields[3]
                if fbc not in bc_set:
                    continue
                for pos in (fstart, fend):
                    if ws <= pos < we:
                        b = (pos - ws) // BIN
                        coverage[status][b] += 1
        except ValueError:
            pass
    for status in coverage:
        coverage[status] = coverage[status] / n_cells_per_group[status] * 1000
    results[g] = coverage
    print(g, 'done, max naive', coverage['naive'].max(), 'max infected', coverage['infected'].max())

for s in tbx_handles.values():
    s.close()

# ---- save everything needed for plotting ----
np.savez(f'{OUT}/all_gene_tracks_coverage.npz',
         **{f'{g}__naive': results[g]['naive'] for g in GENES_U},
         **{f'{g}__infected': results[g]['infected'] for g in GENES_U})

meta = {}
for g in GENES_U:
    chrom, ws, we = windows[g]
    meta[g] = dict(chrom=chrom, win_start=ws, win_end=we, bin=BIN,
                    gene_start=gene_info[g]['start'], gene_end=gene_info[g]['end'],
                    strand=gene_info[g]['strand'], exons=exon_info[g],
                    peaks=peaks_per_gene[g], display_name=GENES[GENES_U.index(g)])
with open(f'{OUT}/all_gene_tracks_meta.json', 'w') as f:
    json.dump(meta, f, indent=1)
print('saved all coverage + meta')
