import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results'

dfs = []
for label, tp in [('S20_Td0_BCG', 'Td0'), ('S20_Tm3_BCG', 'Tm3')]:
    d = pd.read_csv(f'{OUT}/{label}_fate_by_subcluster_percell.csv', index_col=0)
    d['timepoint'] = tp
    dfs.append(d)
df = pd.concat(dfs)

term_cols = ['CMP_a', 'CMP_b', 'CMP_c', 'GMP_a', 'GMP_b', 'MEP_a', 'MEP_b', 'MEP_c', 'MLP_a', 'MLP_b', 'PreBNK']
long = df.reset_index().melt(id_vars=['hsc_clusters_lab', 'timepoint'], value_vars=term_cols,
                              var_name='terminal_state', value_name='fate_probability')

c_order = [f'c{i}' for i in range(10)]
long['hsc_clusters_lab'] = pd.Categorical(long['hsc_clusters_lab'], categories=c_order, ordered=True)

sns.set_style('whitegrid')
g = sns.catplot(
    data=long, x='terminal_state', y='fate_probability', hue='timepoint',
    col='hsc_clusters_lab', col_wrap=5, kind='box', col_order=c_order,
    height=2.6, aspect=1.15, fliersize=1.5, linewidth=0.8,
    palette={'Td0': '#5b7fa6', 'Tm3': '#c0562f'},
    showfliers=False,
)
g.set_xticklabels(term_cols, rotation=90, fontsize=7)
g.set_titles('{col_name}')
g.set_axis_labels('', 'fate probability')

for ax, c in zip(g.axes.flat, c_order):
    if c in ('c2', 'c3', 'c4', 'c7'):
        ax.set_title(c, color='#c0562f', fontweight='bold')
        for spine in ax.spines.values():
            spine.set_edgecolor('#c0562f')
            spine.set_linewidth(1.5)

g.fig.suptitle('Fate probability distribution by terminal state, faceted by HSC subcluster (c0–c9)\nS20 (BCG donor), Td0 vs Tm3 — orange outline = c2/c3/c4/c7',
               fontsize=12.5, fontweight='bold', y=1.04)
g.fig.tight_layout()
g.savefig(f'{OUT}/fate_boxplot_wrap_S20.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')
