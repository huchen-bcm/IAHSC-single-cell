import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

mouse_axis = {'black': 0.803, 'snow': 0.755, 'dimgrey': 0.116, 'lightcoral': -0.364,
              'rosybrown': -0.110, 'gainsboro': -0.027, 'silver': -0.028}
human_axis = {'lightcoral': 0.345, 'silver': 0.120, 'gainsboro': 0.891, 'snow': 0.291,
              'darkgrey': 0.597, 'rosybrown': 0.764, 'black': 0.213, 'brown': 0.385, 'dimgrey': 0.157}

DISPLAY = {'black': 'black', 'snow': 'snow', 'dimgrey': 'dimgrey', 'lightcoral': 'lightcoral',
           'rosybrown': 'rosybrown', 'gainsboro': 'gainsboro', 'silver': 'silver', 'darkgrey': 'darkgrey',
           'brown': 'brown'}

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, axis_dict, title in [(axes[0], mouse_axis, 'Mouse WGCNA modules'), (axes[1], human_axis, 'Human WGCNA modules')]:
    items = sorted(axis_dict.items(), key=lambda kv: -kv[1])
    names = [k for k, v in items]
    vals = [v for k, v in items]
    bar_colors = ['#c0562f' if abs(v) > 0.4 else '#9aa0a6' for v in vals]
    ax.barh(range(len(names)), vals, color=bar_colors)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names)
    ax.invert_yaxis()
    ax.axvline(0, color='black', lw=0.7)
    ax.set_xlabel('Spearman r: module eigengene vs. axis projection')
    ax.set_title(title, fontsize=11, fontweight='bold')
    ax.set_xlim(-1, 1)

axes[0].annotate('↔ human darkgrey\n(103 shared genes,\n1.23× enriched, p=0.003)', xy=(0.755, 1), xytext=(0.2, 2.3),
                  fontsize=8, color='#8a3a0a', arrowprops=dict(arrowstyle='->', color='#8a3a0a', lw=1))
axes[1].annotate('↔ mouse snow', xy=(0.597, 4), xytext=(0.15, 5.3),
                  fontsize=8, color='#8a3a0a', arrowprops=dict(arrowstyle='->', color='#8a3a0a', lw=1))

plt.tight_layout()
plt.savefig(f'{OUT}/assets/wgcna_module_axis_correlation.png', dpi=200, facecolor='white')
print('saved')
