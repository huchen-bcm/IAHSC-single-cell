import anndata as ad
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

latent = ad.read_h5ad(f'{OUT}/embed_latent50.h5ad')
train = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/adata_merged_human_mouse_orthologs_hvg_integrated.0128_2026.h5ad',
                      backed='r')
assert (latent.obs_names == train.obs_names).all()
umap = train.obsm['X_umap']

ia_mask = (latent.obs['HSC_subclass'] == 'IA-HSC').values
lt_mask = (latent.obs['HSC_subclass'] == 'LT-HSC').values
d = latent.X[ia_mask].mean(axis=0) - latent.X[lt_mask].mean(axis=0)
d = d / np.linalg.norm(d)
proj = latent.X @ d

system = latent.obs['system'].values
mouse_mask = system == 'mouse'
human_mask = system == 'human'

fig = plt.figure(figsize=(13.5, 9))
gs = GridSpec(2, 2, height_ratios=[2.1, 1.2], width_ratios=[1, 1], hspace=0.32, wspace=0.18)

vmax = np.percentile(np.abs(proj), 97)

# --- top-left: mouse cells on joint UMAP, colored by projection ---
ax0 = fig.add_subplot(gs[0, 0])
ax0.scatter(umap[human_mask, 0], umap[human_mask, 1], s=2, color='#e3e1da', linewidth=0, rasterized=True)
sc0 = ax0.scatter(umap[mouse_mask, 0], umap[mouse_mask, 1], c=proj[mouse_mask], cmap='RdBu_r',
                   vmin=-vmax, vmax=vmax, s=3, linewidth=0, rasterized=True)
ax0.set_title('Mouse cells (human cells in grey)\ncolored by projection on IA-HSC↔LT-HSC axis', fontsize=11, fontweight='bold')
ax0.set_xticks([]); ax0.set_yticks([])
for spine in ax0.spines.values(): spine.set_visible(False)

# --- top-right: human cells on joint UMAP, colored by projection ---
ax1 = fig.add_subplot(gs[0, 1])
ax1.scatter(umap[mouse_mask, 0], umap[mouse_mask, 1], s=2, color='#e3e1da', linewidth=0, rasterized=True)
sc1 = ax1.scatter(umap[human_mask, 0], umap[human_mask, 1], c=proj[human_mask], cmap='RdBu_r',
                   vmin=-vmax, vmax=vmax, s=3, linewidth=0, rasterized=True)
ax1.set_title('Human cells (mouse cells in grey)\ncolored by the SAME mouse-derived axis', fontsize=11, fontweight='bold')
ax1.set_xticks([]); ax1.set_yticks([])
for spine in ax1.spines.values(): spine.set_visible(False)

cbar = fig.colorbar(sc1, ax=[ax0, ax1], orientation='horizontal', fraction=0.035, pad=0.02, aspect=40)
cbar.set_label('Projection onto mouse (IA-HSC − LT-HSC) latent axis  →  more IA-HSC-like')

# --- bottom: ridgeline/violin of projection by group, mouse + human combined ---
ax2 = fig.add_subplot(gs[1, :])
groups = []
mouse_order = ['IA-HSC', 'ST-HSC', 'LT-HSC', 'Cd41+ HSC']
for g in mouse_order:
    groups.append(('mouse:' + g, proj[mouse_mask][latent.obs.loc[mouse_mask, 'HSC_subclass'].values == g]))
human_order = ['c4', 'c7', 'c3', 'c8', 'c2', 'c5', 'c0', 'c1', 'c6', 'c9']
for g in human_order:
    groups.append(('human:' + g, proj[human_mask][latent.obs.loc[human_mask, 'hsc_clusters_lab'].astype(str).values == g]))

positions = list(range(len(groups)))
data = [g[1] for g in groups]
labels = [g[0] for g in groups]
colors = []
for name, _ in groups:
    if name in ('mouse:IA-HSC',) or name.split(':')[1] in ('c3', 'c4', 'c7'):
        colors.append('#c0562f')
    elif name.startswith('mouse'):
        colors.append('#5b7fa6')
    else:
        colors.append('#9aa0a6')

parts = ax2.violinplot(data, positions=positions, showmeans=True, showextrema=False, widths=0.8)
for i, pc in enumerate(parts['bodies']):
    pc.set_facecolor(colors[i])
    pc.set_alpha(0.75)
    pc.set_edgecolor('none')
parts['cmeans'].set_color('black')
parts['cmeans'].set_linewidth(1.2)

ax2.axvline(3.5, color='black', lw=1, linestyle='--', alpha=0.4)
ax2.axhline(0, color='black', lw=0.6)
ax2.set_xticks(positions)
ax2.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
ax2.set_ylabel('Projection onto axis')
ax2.set_title('Same axis, shown per population (mouse left of dashed line, human right)', fontsize=11, fontweight='bold')

plt.savefig(f'{OUT}/assets/axis_projection_umap.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')
