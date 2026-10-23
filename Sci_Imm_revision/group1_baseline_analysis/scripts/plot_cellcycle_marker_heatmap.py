import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis'

df = pd.read_csv(f'{OUT}/individual_proliferation_marker_stats.csv')

# canonical cell-cycle / proliferation marker panel already established in this project
# (same 10 genes used for the bulk cohort 1 gene-level heatmap, README section 9)
marker_order = ['Mki67', 'Pcna', 'Mcm2', 'Ccna2', 'Cdk1', 'Top2a', 'Ccnb1', 'Cdc20', 'Aurka', 'Birc5']

long = df.melt(id_vars=['HSC_subclass', 'marker'], value_vars=['naive_mean', 'infected_mean'],
               var_name='Treatment', value_name='mean_expr')
long['Treatment'] = long['Treatment'].map({'naive_mean': 'Naive', 'infected_mean': 'Infected'})
long['col'] = long['HSC_subclass'] + ' ' + long['Treatment']

mat = long.pivot(index='marker', columns='col', values='mean_expr')
mat = mat.loc[marker_order, ['LT-HSC Naive', 'LT-HSC Infected', 'IA-HSC Naive', 'IA-HSC Infected']]
mat.columns = ['LT-HSC\nnaive', 'LT-HSC\ninfected', 'IA-HSC\nnaive', 'IA-HSC\ninfected']

vmax = mat.values.max()
fig, ax = plt.subplots(figsize=(5.5, max(1.6, 0.4 * len(mat) + 1)))
ax.imshow(mat.values, cmap='Reds', vmin=0, vmax=vmax, aspect='auto')
ax.set_xticks(range(len(mat.columns)))
ax.set_xticklabels(mat.columns, fontsize=8)
ax.set_yticks(range(len(mat.index)))
ax.set_yticklabels(mat.index, fontsize=8, style='italic')
for i in range(mat.shape[0]):
    for j in range(mat.shape[1]):
        v = mat.values[i, j]
        text_color = 'white' if v / vmax > 0.55 else 'black'
        ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=7, color=text_color)
ax.set_title('Canonical cell-cycle marker expression\n(scRNA, WT only, mean log-normalized expression)', fontsize=9)
for spine in ax.spines.values():
    spine.set_visible(True)
plt.tight_layout()
plt.savefig(f'{OUT}/assets/cellcycle_marker_gene_heatmap_scRNA.pdf', facecolor='white')
plt.savefig(f'{OUT}/assets/cellcycle_marker_gene_heatmap_scRNA.png', dpi=200, facecolor='white')
print('saved')
print(mat)
