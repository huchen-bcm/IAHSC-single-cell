import pandas as pd
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

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

ACCENT = '#eb6834'
NEUTRAL = '#b9b8b3'
INK = '#1a1a19'
INK2 = '#52514e'
PANEL_BORDER = '#d8d7d1'

plt.rcParams['font.family'] = 'DejaVu Sans'

fig = plt.figure(figsize=(11.5, 8.6))
gs = fig.add_gridspec(3, 2, height_ratios=[3.2, 0.85, 3.2], hspace=0.55, wspace=0.30,
                       top=0.80, bottom=0.11, left=0.09, right=0.97)

panel_specs = [
    (gs[0, 0], mouse_full, mouse_mod, 'firebrick', mouse_order, mouse_subclass, 'IA-HSC',
     'Mouse — chaperone/ER-stress module', "IA-HSC's #1-ranked module"),
    (gs[0, 1], human_full, human_mod, 'darkgrey', human_order, human_cluster, ('c3', 'c4', 'c7'),
     'Human — chaperone/ER-stress module', 'r = 0.906 with the validated axis'),
    (gs[2, 0], mouse_full, mouse_mod, 'red', mouse_order, mouse_subclass, None,
     'Mouse — cell-cycle module (control)', 'Not IA-HSC-specific'),
    (gs[2, 1], human_full, human_mod, 'coral', human_order, human_cluster, None,
     'Human — cell-cycle module (control)', 'Not c3/c4/c7-specific'),
]

for spec, adata, modtab, modname, order, labels, highlight, title, subtitle in panel_specs:
    ax = fig.add_subplot(spec)
    genes = modtab[modtab['moduleColors'] == modname].index.tolist()
    s = pd.Series(score(adata, genes), index=adata.obs_names)
    means = s.groupby(labels).mean().reindex(order)
    if highlight is None:
        colors = [NEUTRAL] * len(order)
    elif isinstance(highlight, tuple):
        colors = [ACCENT if g in highlight else NEUTRAL for g in order]
    else:
        colors = [ACCENT if g == highlight else NEUTRAL for g in order]

    ax.bar(range(len(order)), means.values, color=colors, width=0.62, zorder=3)
    ax.axhline(0, color=INK, lw=0.8, zorder=2)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=9.5, color=INK2)
    ax.set_ylabel('Module score (z)', fontsize=9, color=INK2)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    for spine in ('left', 'bottom'):
        ax.spines[spine].set_color(PANEL_BORDER)
    ax.tick_params(colors=INK2, labelsize=9)
    subtitle_color = ACCENT if highlight else INK2
    ax.text(0.0, 1.16, title, transform=ax.transAxes, fontsize=11, fontweight='bold', color=INK, ha='left', va='bottom')
    ax.text(0.0, 1.03, subtitle, transform=ax.transAxes, fontsize=9.3, color=subtitle_color,
            fontweight='bold' if highlight else 'normal', ha='left', va='bottom')

# --- middle annotation strip: shared genes, spanning both top columns ---
ann_ax = fig.add_subplot(gs[1, :])
ann_ax.axis('off')
ann_ax.set_xlim(0, 1)
ann_ax.set_ylim(0, 1)
ann_ax.plot([0.08, 0.92], [0.85, 0.85], color=ACCENT, lw=1.3, transform=ann_ax.transAxes)
ann_ax.plot([0.08, 0.08], [0.85, 0.55], color=ACCENT, lw=1.3, transform=ann_ax.transAxes)
ann_ax.plot([0.92, 0.92], [0.85, 0.55], color=ACCENT, lw=1.3, transform=ann_ax.transAxes)
ann_ax.text(0.5, 0.40,
            'Shared genes: HSP90B1, HSPA5, BCL2, PRKCA, RGCC, TMSB4X — enrich for ATF6/chaperone response &\n'
            'ER unfolded-protein response (independently confirmed by GO in BOTH modules’ full gene lists)',
            transform=ann_ax.transAxes, ha='center', va='center', fontsize=9.5, color=ACCENT, fontweight='bold')

# --- headline ---
fig.text(0.09, 0.965, 'A conserved chaperone/ER-stress module links mouse IA-HSC to human c3/c4/c7',
          fontsize=16, fontweight='bold', color=INK, ha='left', va='top')
fig.text(0.09, 0.925,
          'Found by an exhaustive search of all 18 mouse × 24 human WGCNA modules for pairs that are both\n'
          'significantly gene-overlapping across species and independently validated against the sysVI cross-species axis.\n'
          'A well-matched cell-cycle pair (bottom row) passes the first test but fails the second — shown for contrast.',
          fontsize=10, color=INK2, ha='left', va='top')

# --- legend ---
leg_y = 0.045
fig.patches.append(Rectangle((0.09, leg_y), 0.018, 0.018, color=ACCENT, transform=fig.transFigure))
fig.text(0.115, leg_y + 0.009, 'IA-HSC / c3, c4, c7 (the populations the reviewer asked about)',
          fontsize=9.3, va='center', color=INK)
fig.patches.append(Rectangle((0.58, leg_y), 0.018, 0.018, color=NEUTRAL, transform=fig.transFigure))
fig.text(0.605, leg_y + 0.009, 'all other subtypes / clusters', fontsize=9.3, va='center', color=INK)

plt.savefig(f'{OUT}/assets/figure_reviewer_final.png', dpi=220, facecolor='white')
print('saved')
