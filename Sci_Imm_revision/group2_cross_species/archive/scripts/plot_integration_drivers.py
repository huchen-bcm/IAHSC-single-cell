import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from adjustText import adjust_text

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

res = pd.read_csv(f'{OUT}/integration_driver_genes_full.csv')
go_pos = pd.read_csv(f'{OUT}/integration_driver_GO_positive_IA_c347.csv')

fig = plt.figure(figsize=(13, 5.2))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.1], wspace=0.28)

# --- Panel A: axis validation bar chart ---
ax0 = fig.add_subplot(gs[0])
projs = {'c4': 1.705, 'c7': 1.642, 'c3': 1.456, 'c8': 0.743, 'c2': 0.732,
         'c5': -0.379, 'c0': -1.055, 'c1': -1.397, 'c6': -3.482, 'c9': -4.657}
order = sorted(projs, key=lambda k: -projs[k])
colors = ['#c0562f' if k in ('c3', 'c4', 'c7') else '#9aa0a6' for k in order]
ax0.barh(range(len(order)), [projs[k] for k in order], color=colors)
ax0.set_yticks(range(len(order)))
ax0.set_yticklabels(order)
ax0.invert_yaxis()
ax0.axvline(0, color='black', lw=0.8)
ax0.set_xlabel('Projection onto mouse (IA-HSC − LT-HSC) latent axis')
ax0.set_title('Axis defined ONLY from mouse data\ncorrectly separates c3/c4/c7 in human data', fontsize=11, fontweight='bold')

# --- Panel B: gene scatter ---
ax1 = fig.add_subplot(gs[1])
concordant = np.sign(res['mouse_corr']) == np.sign(res['human_corr'])
ax1.scatter(res.loc[~concordant, 'mouse_corr'], res.loc[~concordant, 'human_corr'],
            s=6, color='#cccccc', alpha=0.5, linewidth=0, label='discordant')
ax1.scatter(res.loc[concordant, 'mouse_corr'], res.loc[concordant, 'human_corr'],
            s=8, color='#c0562f', alpha=0.6, linewidth=0, label='concordant')
ax1.axhline(0, color='black', lw=0.5)
ax1.axvline(0, color='black', lw=0.5)
ax1.set_xlabel("Spearman r with axis, mouse cells")
ax1.set_ylabel("Spearman r with axis, human cells")
ax1.set_title('Per-gene correlation with the IA-HSC/c3-4-7 axis,\ncomputed independently per species', fontsize=11, fontweight='bold')

top_genes = res[concordant].reindex(res[concordant]['min_abs_corr'].abs().sort_values(ascending=False).index).head(16)
texts = []
for _, row in top_genes.iterrows():
    texts.append(ax1.text(row['mouse_corr'], row['human_corr'], row['gene'], fontsize=7.5))
adjust_text(texts, ax=ax1, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
ax1.legend(fontsize=8, loc='lower right')

plt.tight_layout()
plt.savefig(f'{OUT}/assets/integration_driver_genes.png', dpi=200, facecolor='white')
print('saved')
