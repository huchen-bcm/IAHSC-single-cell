import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

sim = pd.read_csv(f'{OUT}/figure1_state_correspondence_cosine.csv', index_col=0)
mouse_order = sim.index.tolist()
human_order = sim.columns.tolist()

fig, ax = plt.subplots(figsize=(9, 3.6))
im = ax.imshow(sim.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')

for i in range(sim.shape[0]):
    for j in range(sim.shape[1]):
        v = sim.values[i, j]
        ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=9,
                 color='white' if abs(v) > 0.55 else 'black',
                 fontweight='bold' if v == sim.values[i].max() else 'normal')

ax.set_xticks(range(len(human_order)))
ax.set_xticklabels(human_order, fontsize=10)
ax.set_yticks(range(len(mouse_order)))
ax.set_yticklabels(mouse_order, fontsize=10)
for xt, c in zip(ax.get_xticklabels(), human_order):
    if c in ('c3', 'c4', 'c7'):
        xt.set_color('#c0562f')
        xt.set_fontweight('bold')
ax.set_title('Figure 1 — Mouse HSC-subtype × human-cluster correspondence\n(cosine similarity, sysVI joint latent space, 50-dim)', fontsize=11, fontweight='bold')
cbar = plt.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
cbar.set_label('cosine similarity')

plt.tight_layout()
plt.savefig(f'{OUT}/assets/figure1_state_correspondence.png', dpi=200, facecolor='white')
print('saved figure 1')
