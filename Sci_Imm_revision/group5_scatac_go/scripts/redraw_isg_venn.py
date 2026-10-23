import json
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
from scipy.stats import hypergeom

OUT_G5 = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group5_scatac_go'

with open(f'{OUT_G5}/atac_gene_sets.json') as f:
    gene_sets = json.load(f)
atac_genes = set(g.upper() for g in gene_sets['Infected_vs_Naive.opening_in_Infected'])

deg = pd.read_csv('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/DEG_pathway/singlecell_de.WT/IA-HSC.WT.Infected_vs_Naive.singlecell_de_DEG.tsv', sep='\t')
rna_up = set(deg.loc[(deg['log2FoldChange'] > 1) & (deg['padj'] < 0.05), 'name'].str.upper())

N = len(deg)
overlap = atac_genes & rna_up
n_atac = len(atac_genes)
n_rna = len(rna_up)
n_overlap = len(overlap)
expected = n_atac * n_rna / N
fold = n_overlap / expected
pval = hypergeom.sf(n_overlap - 1, N, n_rna, n_atac)

print(f'ATAC opening n={n_atac}, RNA up n={n_rna}, N={N}')
print(f'overlap={n_overlap}, expected={expected:.1f}, fold={fold:.2f}x, p={pval:.2e}')
print(sorted(overlap))

isg_panel = {'IRF1', 'ISG15', 'IGTP', 'NLRC5', 'PSMB8', 'PSMB9', 'TAP1', 'BST2', 'SLFN9'}
in_panel = sorted(overlap & isg_panel)
extra = sorted(overlap - isg_panel)

fig, ax = plt.subplots(figsize=(11, 6))
ax.set_xlim(0, 11)
ax.set_ylim(0, 6)
ax.axis('off')

c1 = Circle((3.3, 3), 2.6, color='#7b6fc4', alpha=0.55)
c2 = Circle((5.3, 3), 2.6, color='#e0a084', alpha=0.55)
ax.add_patch(c1)
ax.add_patch(c2)

ax.text(2.0, 3, f'{n_atac - n_overlap}', ha='center', va='center', fontsize=14)
ax.text(4.3, 3, f'{n_overlap}', ha='center', va='center', fontsize=14)
ax.text(6.6, 3, f'{n_rna - n_overlap}', ha='center', va='center', fontsize=14)

ax.text(2.3, 0.3, f'ATAC: opening\nwith infection\n(n={n_atac})', ha='center', fontsize=11)
ax.text(6.3, 0.3, f'scRNA: IA-HSC\nup with infection\n(n={n_rna})', ha='center', fontsize=11)

ax.text(0.2, 5.6, f'Fold enrichment = {fold:.2f}×  (expected {expected:.1f}, observed {n_overlap})\n'
                   f'p = {pval:.1e} (hypergeometric, N={N})', fontsize=12, ha='left')

ax.text(8.3, 5.3, f'Overlap genes (n={n_overlap}):', fontsize=13, fontweight='bold')
ax.text(8.3, 4.85, 'Interferon-stimulated gene panel (reported):', fontsize=10.5, color='#4a3d8f', fontweight='bold')
for i, chunk_start in enumerate(range(0, len(in_panel), 3)):
    chunk = in_panel[chunk_start:chunk_start + 3]
    ax.text(8.45, 4.55 - i * 0.32, ', '.join(g.capitalize() for g in chunk), fontsize=10)

y0 = 4.55 - ((len(in_panel) + 2) // 3) * 0.32 - 0.25
ax.text(8.3, y0, f'Additional overlap genes (n={len(extra)}):', fontsize=10.5, color='#b5502a', fontweight='bold')
for i, chunk_start in enumerate(range(0, len(extra), 3)):
    chunk = extra[chunk_start:chunk_start + 3]
    ax.text(8.45, y0 - 0.32 - i * 0.32, ', '.join(g.capitalize() for g in chunk), fontsize=10)

ax.set_title('ATAC-seq peaks opening with infection vs. scRNA-seq IA-HSC infection-response DEGs', fontsize=13, pad=15)

plt.tight_layout()
plt.savefig(f'{OUT_G5}/assets/atac_rna_isg_overlap_redrawn.pdf', facecolor='white', bbox_inches='tight')
plt.savefig(f'{OUT_G5}/assets/atac_rna_isg_overlap_redrawn.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')
