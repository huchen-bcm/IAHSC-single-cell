import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results'
hsc_lab = pd.read_csv('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/GEO_data/sc_pdata_phate_ndim2_knn6_t10_gamma1.csv', index_col=0)

samples = [('S20_Td0_BCG', 'S20 (BCG) Td0'), ('S20_Tm3_BCG', 'S20 (BCG) Tm3')]
c_order = [f'c{i}' for i in range(10)]

fig, axes = plt.subplots(1, len(samples), figsize=(6.5 * len(samples), 5))
if len(samples) == 1:
    axes = [axes]

for ax, (label, title) in zip(axes, samples):
    tp = pd.read_csv(f'{OUT}/{label}_to_terminal_states.csv')
    tp = tp.set_index('Cell ID')
    merged = tp.join(hsc_lab['hsc_clusters_lab'], how='inner')
    means = merged.groupby('hsc_clusters_lab').mean(numeric_only=True).reindex(c_order)

    im = ax.imshow(means.values, cmap='viridis', vmin=0, vmax=means.values.max(), aspect='auto')
    ax.set_xticks(range(means.shape[1]))
    ax.set_xticklabels(means.columns, rotation=45, ha='right', fontsize=8)
    ax.set_yticks(range(len(c_order)))
    ax.set_yticklabels(c_order, fontsize=9)
    for yt, c in zip(ax.get_yticklabels(), c_order):
        if c in ('c2', 'c3', 'c4', 'c7'):
            yt.set_color('#c0562f')
            yt.set_fontweight('bold')
    for i in range(means.shape[0]):
        for j in range(means.shape[1]):
            v = means.values[i, j]
            if v > 0.05:
                ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=7,
                         color='white' if v > means.values.max() * 0.55 else 'black')
    ax.set_title(title, fontsize=11, fontweight='bold')

fig.suptitle('Fate probability by HSC subcluster (c0–c9) — CellRank reproduction, one BCG donor',
             fontsize=12.5, fontweight='bold', y=1.03)
plt.tight_layout()
plt.savefig(f'{OUT}/fate_bias_by_subcluster_S20.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')
