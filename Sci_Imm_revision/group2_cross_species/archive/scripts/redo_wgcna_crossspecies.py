import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import hypergeom

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

mouse_mod = pd.read_csv(f'{OUT}/wgcna_mouse_fullgene_gene_modules.csv', index_col=0)
human_mod = pd.read_csv(f'{OUT}/wgcna_human_fullgene_gene_modules.csv', index_col=0)
mouse_r = pd.read_csv(f'{OUT}/assets/wgcna_fullgene_module_trait_mouse_r.csv', index_col=0)
human_r = pd.read_csv(f'{OUT}/assets/wgcna_fullgene_module_trait_human_r.csv', index_col=0)

mouse_genes = set(g.upper() for g in mouse_mod.index)
human_genes = set(g.upper() for g in human_mod.index)
shared_universe = mouse_genes & human_genes
N = len(shared_universe)

results = []
for mmod, mgenes in mouse_mod.groupby('moduleColors'):
    mset = set(g.upper() for g in mgenes.index) & shared_universe
    for hmod, hgenes in human_mod.groupby('moduleColors'):
        hset = set(g.upper() for g in hgenes.index) & shared_universe
        overlap = mset & hset
        k = len(overlap)
        if k == 0:
            continue
        M, n, K = N, len(mset), len(hset)
        pval = hypergeom.sf(k - 1, M, n, K)
        expected = n * K / M
        fold = k / expected if expected > 0 else np.nan
        results.append((mmod, hmod, len(mset), len(hset), k, expected, fold, pval))

res = pd.DataFrame(results, columns=['mouse_mod', 'human_mod', 'n_mouse', 'n_human', 'overlap', 'expected', 'fold', 'p'])
sig = res[res['p'] < 0.01].sort_values('p').copy()
sig['mouse_axis_r'] = sig['mouse_mod'].map(lambda m: mouse_r.loc['ME' + m, 'axis_proj'])
sig['human_axis_r'] = sig['human_mod'].map(lambda m: human_r.loc['ME' + m, 'axis_proj'])
sig.to_csv(f'{OUT}/wgcna_fullgene_significant_crossspecies_pairs_REDONE.csv', index=False)
print(f'{len(sig)} significant pairs at p<0.01 (background N={N}, shared testable gene universe)')
print(sig.to_string())

# --- figure, same style as figure2_wgcna_correspondence.png, but honest/current numbers ---
fig, ax = plt.subplots(figsize=(9, 5.5))
x = sig['mouse_axis_r'].values
y = sig['human_axis_r'].values
size = -np.log10(sig['p'].values) * 25
concordant = (x > 0.3) & (y > 0.3)

ax.scatter(x[~concordant], y[~concordant], s=size[~concordant], color='#9aa0a6', alpha=0.7,
           edgecolor='black', linewidth=0.4, zorder=2, label='significant overlap, not axis-concordant')
if concordant.any():
    ax.scatter(x[concordant], y[concordant], s=size[concordant], color='#c0562f', alpha=0.9,
               edgecolor='black', linewidth=0.8, zorder=3, label='significant overlap AND axis-concordant')

ax.axhline(0.3, color='grey', linestyle='--', lw=0.7)
ax.axvline(0.3, color='grey', linestyle='--', lw=0.7)
ax.axhline(0, color='black', lw=0.5)
ax.axvline(0, color='black', lw=0.5)

for _, row in sig.iterrows():
    lbl = f"{row['mouse_mod']}↔{row['human_mod']}"
    ax.annotate(lbl, (row['mouse_axis_r'], row['human_axis_r']), xytext=(4, 4),
                textcoords='offset points', fontsize=8)

ax.set_xlabel("Mouse module's own correlation with the validated IA-HSC axis", fontsize=10)
ax.set_ylabel("Human module's own correlation with the validated axis", fontsize=10)
ax.set_title(f'REDONE: cross-species module overlap × axis relevance\n'
             f'(all {len(sig)} significant module pairs, corrected background N={N}, current WGCNA run)',
             fontsize=11, fontweight='bold')
ax.legend(fontsize=8, loc='lower right')
ax.set_xlim(-0.5, 0.9)
ax.set_ylim(-0.3, 0.7)

plt.tight_layout()
plt.savefig(f'{OUT}/assets/figure2_wgcna_correspondence_REDONE.pdf', facecolor='white')
plt.savefig(f'{OUT}/assets/figure2_wgcna_correspondence_REDONE.png', dpi=200, facecolor='white')
print('\nsaved figure2_wgcna_correspondence_REDONE.{pdf,png}')
print(f'\npairs with BOTH axis correlations > 0.3: {concordant.sum()}')
