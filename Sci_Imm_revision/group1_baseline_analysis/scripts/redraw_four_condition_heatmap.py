import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis'
df = pd.read_csv(f'{OUT}/four_condition_group_means_pseudobulk_gsva.csv', index_col=0)
cols = ['LT-HSC\nnaive', 'IA-HSC\nnaive', 'LT-HSC\ninfected', 'IA-HSC\ninfected']
df = df[cols]

vmax = np.abs(df.values).max()
fig, ax = plt.subplots(figsize=(8.5, 8.5))
im = ax.imshow(df.values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')

ax.set_xticks(range(len(cols)))
ax.set_xticklabels(cols)
ax.set_yticks(range(len(df.index)))
ax.set_yticklabels(df.index)

for i in range(df.shape[0]):
    for j in range(df.shape[1]):
        v = df.values[i, j]
        ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=10,
                color='black')

ax.set_title('Four-condition pathway/cell-cycle GSVA scores\n(scRNA pseudobulk, same algorithm as bulk)', fontsize=12)
for spine in ax.spines.values():
    spine.set_visible(True)
plt.tight_layout()
plt.savefig(f'{OUT}/assets/four_condition_heatmap_pseudobulk_gsva_redrawn.pdf', facecolor='white')
plt.savefig(f'{OUT}/assets/four_condition_heatmap_pseudobulk_gsva_redrawn.png', dpi=200, facecolor='white')
print('saved')
