import pandas as pd
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

mouse_mod = pd.read_csv(f'{OUT}/wgcna_mouse_fullgene_gene_modules.csv', index_col=0)
human_mod = pd.read_csv(f'{OUT}/wgcna_human_fullgene_gene_modules.csv', index_col=0)
mouse_full = ad.read_h5ad(f'{OUT}/mouse_fullgene_hvg6000.h5ad')
human_full = ad.read_h5ad(f'{OUT}/human_fullgene_hvg6000.h5ad')
latent = ad.read_h5ad(f'{OUT}/embed_latent50.h5ad')
mouse_subclass = latent.obs.loc[mouse_full.obs_names, 'HSC_subclass'].astype(str).values
human_cluster = latent.obs.loc[human_full.obs_names, 'hsc_clusters_lab'].astype(str).values


def score(adata, genes):
    idx = [adata.var_names.get_loc(g) for g in genes if g in adata.var_names]
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    s = X[:, idx].mean(axis=1)
    return (s - s.mean()) / s.std()


mouse_order = ['IA-HSC', 'LT-HSC', 'ST-HSC', 'Cd41+ HSC']
human_order = [f'c{i}' for i in range(10)]

mouse_targets = {'red': 'Cell cycle\n(Mitotic Anaphase)', 'firebrick': 'Chaperone + platelet\n(Hemostasis, Hsp90b1/Hspa5)'}
human_targets = {'darkgrey': 'Cytokine/IL signaling\n(r=0.906 with axis)', 'coral': 'Cell cycle\n(shares genes w/ mouse red+black)'}

fig, axes = plt.subplots(2, 2, figsize=(13, 8))

for ax, (mod, label) in zip(axes[0], mouse_targets.items()):
    genes = mouse_mod[mouse_mod['moduleColors'] == mod].index.tolist()
    s = pd.Series(score(mouse_full, genes), index=mouse_full.obs_names)
    means = s.groupby(mouse_subclass).mean().loc[mouse_order]
    colors = ['#c0562f' if g == 'IA-HSC' else '#9aa0a6' for g in mouse_order]
    ax.bar(range(4), means.values, color=colors)
    ax.set_xticks(range(4))
    ax.set_xticklabels(mouse_order, rotation=20, ha='right')
    ax.axhline(0, color='black', lw=0.6)
    ax.set_title(f'Mouse "{mod}" (n={len(genes)})\n{label}', fontsize=10, fontweight='bold')
    ax.set_ylabel('z-score')

for ax, (mod, label) in zip(axes[1], human_targets.items()):
    genes = human_mod[human_mod['moduleColors'] == mod].index.tolist()
    s = pd.Series(score(human_full, genes), index=human_full.obs_names)
    means = s.groupby(human_cluster).mean().loc[human_order]
    colors = ['#c0562f' if g in ('c3', 'c4', 'c7') else '#9aa0a6' for g in human_order]
    ax.bar(range(10), means.values, color=colors)
    ax.set_xticks(range(10))
    ax.set_xticklabels(human_order, rotation=0)
    ax.axhline(0, color='black', lw=0.6)
    ax.set_title(f'Human "{mod}" (n={len(genes)})\n{label}', fontsize=10, fontweight='bold')
    ax.set_ylabel('z-score')

plt.tight_layout()
plt.savefig(f'{OUT}/assets/wgcna_fullgene_top_modules.png', dpi=200, facecolor='white')
print('saved')
