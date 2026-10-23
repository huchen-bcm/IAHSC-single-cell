import pandas as pd
import numpy as np
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, Rectangle
import matplotlib.font_manager as fm

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

ACCENT = '#eb6834'    # hypothesis-relevant (IA-HSC / c3,c4,c7)
NEUTRAL = '#b9b8b3'   # everything else
INK = '#1a1a19'
INK2 = '#52514e'
PANEL_BORDER = '#e1e0d9'

plt.rcParams['font.family'] = 'DejaVu Sans'

fig = plt.figure(figsize=(11.5, 7.6))
gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.32, top=0.80, bottom=0.10, left=0.08, right=0.97)

panels = [
    (gs[0, 0], mouse_full, mouse_mod, 'firebrick', mouse_order, mouse_subclass, 'IA-HSC',
     'Mouse — chaperone / platelet module', 'IA-HSC’s #1-ranked module'),
    (gs[0, 1], human_full, human_mod, 'darkgrey', human_order, human_cluster, ('c3', 'c4', 'c7'),
     'Human — cytokine / interleukin-signaling module', 'r = 0.906 with the validated axis'),
    (gs[1, 0], mouse_full, mouse_mod, 'red', mouse_order, mouse_subclass, None,
     'Mouse — cell-cycle module (control)', 'Not IA-HSC-specific'),
    (gs[1, 1], human_full, human_mod, 'coral', human_order, human_cluster, None,
     'Human — cell-cycle module (control)', 'Not c3/c4/c7-specific'),
]

axes = []
for spec, adata, modtab, modname, order, labels, highlight, title, subtitle in panels:
    ax = fig.add_subplot(spec)
    axes.append(ax)
    genes = modtab[modtab['moduleColors'] == modname].index.tolist()
    s = pd.Series(score(adata, genes), index=adata.obs_names)
    means = s.groupby(labels).mean().reindex(order)
    if highlight is None:
        colors = [NEUTRAL] * len(order)
    elif isinstance(highlight, tuple):
        colors = [ACCENT if g in highlight else NEUTRAL for g in order]
    else:
        colors = [ACCENT if g == highlight else NEUTRAL for g in order]

    bars = ax.bar(range(len(order)), means.values, color=colors, width=0.62, zorder=3)
    ax.axhline(0, color=INK, lw=0.8, zorder=2)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, fontsize=9.5, color=INK2)
    ax.set_ylabel('Module score (z)', fontsize=9, color=INK2)
    for spine in ('top', 'right'):
        ax.spines[spine].set_visible(False)
    for spine in ('left', 'bottom'):
        ax.spines[spine].set_color(PANEL_BORDER)
    ax.tick_params(colors=INK2, labelsize=9)
    ax.set_title(title, fontsize=11, fontweight='bold', color=INK, loc='left', pad=18)
    ax.text(0, 1.14, subtitle, transform=ax.transAxes, fontsize=9, color=ACCENT if highlight else INK2,
             fontweight='bold' if highlight else 'normal', ha='left')

# --- connect the two validated panels with a bracket + shared-gene annotation ---
fig.canvas.draw()
p0 = axes[0].get_position()
p1 = axes[1].get_position()
y_top = max(p0.y1, p1.y1) + 0.045
fig.lines.append(matplotlib.lines.Line2D([p0.x0 + p0.width * 0.5, p1.x0 + p1.width * 0.5], [y_top, y_top],
                                          transform=fig.transFigure, color=ACCENT, lw=1.4))
fig.lines.append(matplotlib.lines.Line2D([p0.x0 + p0.width * 0.5, p0.x0 + p0.width * 0.5], [y_top, p0.y1 + 0.045],
                                          transform=fig.transFigure, color=ACCENT, lw=1.4))
fig.lines.append(matplotlib.lines.Line2D([p1.x0 + p1.width * 0.5, p1.x0 + p1.width * 0.5], [y_top, p1.y1 + 0.045],
                                          transform=fig.transFigure, color=ACCENT, lw=1.4))
fig.text((p0.x0 + p0.width * 0.5 + p1.x0 + p1.width * 0.5) / 2, y_top + 0.012,
          'Shared genes: HSP90B1, HSPA5, BCL2, PRKCA, RGCC, TMSB4X  (3.3× enrichment, p=0.008)',
          ha='center', va='bottom', fontsize=9.5, color=ACCENT, fontweight='bold')

# --- headline title + subtitle ---
fig.text(0.08, 0.955, 'A cross-species-conserved module links mouse IA-HSC to human c3/c4/c7',
          fontsize=16, fontweight='bold', color=INK, ha='left')
fig.text(0.08, 0.925,
          'Found by an exhaustive search of all 18 mouse × 24 human WGCNA modules for pairs that are BOTH\n'
          'significantly gene-overlapping across species AND independently validated against the sysVI cross-species axis.\n'
          'A well-matched cell-cycle pair (bottom row) passes the first test but fails the second — shown for contrast.',
          fontsize=10, color=INK2, ha='left', va='top')

# legend
legend_ax = fig.add_axes([0.08, 0.025, 0.5, 0.03])
legend_ax.axis('off')
legend_ax.add_patch(Rectangle((0, 0.3), 0.03, 0.5, color=ACCENT, transform=legend_ax.transAxes))
legend_ax.text(0.04, 0.5, 'IA-HSC / c3, c4, c7 (the populations the reviewer asked about)', transform=legend_ax.transAxes,
               fontsize=9, va='center', color=INK)
legend_ax.add_patch(Rectangle((0.62, 0.3), 0.03, 0.5, color=NEUTRAL, transform=legend_ax.transAxes))
legend_ax.text(0.66, 0.5, 'all other subtypes / clusters', transform=legend_ax.transAxes,
               fontsize=9, va='center', color=INK)

plt.savefig(f'{OUT}/assets/figure_reviewer_final.png', dpi=220, facecolor='white', bbox_inches='tight')
print('saved')
