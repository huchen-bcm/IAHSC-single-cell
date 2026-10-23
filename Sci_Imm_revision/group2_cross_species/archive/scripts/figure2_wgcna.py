import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

mouse_axis = {'darkgrey': -0.295, 'red': 0.785, 'snow': -0.178, 'black': 0.088, 'lightgrey': -0.221,
              'mistyrose': -0.363, 'rosybrown': -0.373, 'firebrick': 0.735, 'gainsboro': -0.290}
human_axis = {'lightgrey': 0.489, 'coral': 0.228, 'lightsalmon': 0.092, 'sienna': 0.260, 'darkgrey': 0.906, 'snow': 0.007}

df = pd.read_csv(f'{OUT}/wgcna_fullgene_significant_crossspecies_pairs.csv')
df['mouse_axis_r'] = df['mouse_mod'].map(mouse_axis)
df['human_axis_r'] = df['human_mod'].map(human_axis)
df = df.sort_values('p')

fig, ax = plt.subplots(figsize=(9, 5.5))

x = df['mouse_axis_r'].values
y = df['human_axis_r'].values
size = -np.log10(df['p'].values) * 25
concordant = (x > 0.3) & (y > 0.3)

ax.scatter(x[~concordant], y[~concordant], s=size[~concordant], color='#9aa0a6', alpha=0.7,
           edgecolor='black', linewidth=0.4, zorder=2, label='significant overlap, not axis-concordant')
ax.scatter(x[concordant], y[concordant], s=size[concordant], color='#c0562f', alpha=0.9,
           edgecolor='black', linewidth=0.8, zorder=3, label='significant overlap AND axis-concordant')

ax.axhline(0.3, color='grey', linestyle='--', lw=0.7)
ax.axvline(0.3, color='grey', linestyle='--', lw=0.7)
ax.axhline(0, color='black', lw=0.5)
ax.axvline(0, color='black', lw=0.5)

for _, row in df.iterrows():
    lbl = f"{row['mouse_mod']}↔{row['human_mod']}"
    dx, dy = 0.02, 0.02
    if row['mouse_mod'] == 'firebrick' and row['human_mod'] == 'darkgrey':
        ax.annotate(lbl, (row['mouse_axis_r'], row['human_axis_r']), xytext=(row['mouse_axis_r'] - 0.28, row['human_axis_r'] - 0.05),
                    fontsize=9, fontweight='bold', color='#8a3a0a')
    elif row['mouse_mod'] == 'red' and row['human_mod'] == 'coral':
        ax.annotate(lbl, (row['mouse_axis_r'], row['human_axis_r']), xytext=(row['mouse_axis_r'] - 0.15, row['human_axis_r'] + 0.05),
                    fontsize=8, color='#555555')

ax.set_xlabel("Mouse module's own correlation with the validated IA-HSC axis", fontsize=10)
ax.set_ylabel("Human module's own correlation with the validated axis", fontsize=10)
ax.set_title('Figure 2 (WGCNA version) — cross-species module overlap × axis relevance\n(all 13 significant module pairs from the exhaustive 18×24 search)',
             fontsize=11.5, fontweight='bold')
ax.legend(fontsize=8, loc='lower right')
ax.set_xlim(-0.5, 0.9)
ax.set_ylim(-0.1, 1.0)

plt.tight_layout()
plt.savefig(f'{OUT}/assets/figure2_wgcna_correspondence.png', dpi=200, facecolor='white')
print('saved')
