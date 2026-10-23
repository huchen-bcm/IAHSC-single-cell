import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
NMF_DIR = f'{OUT}/nmf_metaprogram'

mouse_usage = pd.read_csv(f'{NMF_DIR}/mouse_usage.csv', index_col=0)
human_usage = pd.read_csv(f'{NMF_DIR}/human_usage.csv', index_col=0)
mouse_obs = pd.read_csv(f'{NMF_DIR}/mouse_obs.csv', index_col=0)
human_obs = pd.read_csv(f'{NMF_DIR}/human_obs.csv', index_col=0)

mu = mouse_usage.div(mouse_usage.sum(axis=1), axis=0).join(mouse_obs['HSC_subclass'])
hu = human_usage.div(human_usage.sum(axis=1), axis=0).join(human_obs['hsc_clusters_lab'])

mouse_order = ['IA-HSC', 'LT-HSC', 'ST-HSC', 'Cd41+ HSC']
human_order = [f'c{i}' for i in range(10)]

mouse_mean = mu.groupby('HSC_subclass').mean().loc[mouse_order]
human_mean = hu.groupby('hsc_clusters_lab').mean().loc[human_order]

# the 4 conserved pairs, in Figure 2/3 order, with their pathway label
pairs = [
    ('mouseMP1', 'humanMP6', 'Cell Cycle /\nDNA replication'),
    ('mouseMP3', 'humanMP2', 'Ribosomal /\ntranslation'),
    ('mouseMP7', 'humanMP5', 'Neutrophil degran. /\ninnate immune'),
    ('mouseMP2', 'humanMP1', 'Chromatin\norganization'),
]

fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), gridspec_kw={'width_ratios': [4, 10]})

mouse_mat = np.array([mouse_mean[m].values for m, h, lbl in pairs])
human_mat = np.array([human_mean[h].values for m, h, lbl in pairs])

vmax = max(mouse_mat.max(), human_mat.max())

im0 = axes[0].imshow(mouse_mat, cmap='viridis', vmin=0, vmax=vmax, aspect='auto')
axes[0].set_xticks(range(4))
axes[0].set_xticklabels(mouse_order, rotation=30, ha='right', fontsize=9)
axes[0].set_yticks(range(4))
axes[0].set_yticklabels([f'{lbl}\n({m}↔{h})' for m, h, lbl in pairs], fontsize=8)
axes[0].set_title('Mouse subtype', fontsize=10.5, fontweight='bold')
for i in range(mouse_mat.shape[0]):
    for j in range(mouse_mat.shape[1]):
        v = mouse_mat[i, j]
        axes[0].text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=8,
                      color='white' if v < vmax * 0.6 else 'black')

im1 = axes[1].imshow(human_mat, cmap='viridis', vmin=0, vmax=vmax, aspect='auto')
axes[1].set_xticks(range(10))
axes[1].set_xticklabels(human_order, fontsize=9)
axes[1].set_yticks([])
axes[1].set_title('Human cluster', fontsize=10.5, fontweight='bold')
for i in range(human_mat.shape[0]):
    for j in range(human_mat.shape[1]):
        v = human_mat[i, j]
        axes[1].text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=8,
                      color='white' if v < vmax * 0.6 else 'black')
for xt, c in zip(axes[1].get_xticklabels(), human_order):
    if c in ('c3', 'c4', 'c7'):
        xt.set_color('#c0562f')
        xt.set_fontweight('bold')

plt.colorbar(im1, ax=axes[1], fraction=0.025, pad=0.02, label='normalized program usage')
fig.suptitle('Figure 4 — Where do the four conserved programs (Fig. 2–3) actually get used?',
             fontsize=12, fontweight='bold', y=1.04)
plt.tight_layout()
plt.savefig(f'{OUT}/assets/figure4_conserved_program_usage.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved figure 4')
