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

# Sun et al. Fig 3D grouping: green = CMPs/GMPs, blue = MLPs/PreBNK, orange = MEPs
GROUPS = {
    'CMPs/GMPs': ['CMP_a', 'CMP_b', 'CMP_c', 'GMP_a', 'GMP_b'],
    'MLPs/PreBNK': ['MLP_a', 'MLP_b', 'PreBNK'],
    'MEPs': ['MEP_a', 'MEP_b', 'MEP_c'],
}
GROUP_COLOR = {'CMPs/GMPs': '#3a8f4a', 'MLPs/PreBNK': '#3a6fa8', 'MEPs': '#d98a2b'}

for g, cols in GROUPS.items():
    df[g] = df[cols].sum(axis=1)

c_order = [f'c{i}' for i in range(10)]
df['hsc_clusters_lab'] = pd.Categorical(df['hsc_clusters_lab'], categories=c_order, ordered=True)

long = df.reset_index().melt(id_vars=['hsc_clusters_lab', 'timepoint'], value_vars=list(GROUPS.keys()),
                              var_name='fate_group', value_name='probability')
long['fate_group'] = pd.Categorical(long['fate_group'], categories=list(GROUPS.keys()), ordered=True)

sns.set_style('whitegrid')
g = sns.catplot(
    data=long, x='fate_group', y='probability', hue='timepoint',
    row='hsc_clusters_lab', row_order=c_order, kind='box',
    height=1.5, aspect=4.2, fliersize=1.5, linewidth=0.8,
    palette={'Td0': '#5b7fa6', 'Tm3': '#c0562f'},
    showfliers=False, legend_out=True,
)
g.set_axis_labels('', '')
g.set_titles('')
g.fig.supylabel('fate probability', fontsize=11, x=0.01)

for ax, c in zip(g.axes.flat, c_order):
    label_color = '#c0562f' if c in ('c2', 'c3', 'c4', 'c7') else '#333333'
    fw = 'bold' if c in ('c2', 'c3', 'c4', 'c7') else 'normal'
    ax.text(-0.16, 0.5, c, transform=ax.transAxes, ha='right', va='center',
            fontsize=13, color=label_color, fontweight=fw)
    ax.set_ylim(0, 1)
    if c in ('c2', 'c3', 'c4', 'c7'):
        for spine in ax.spines.values():
            spine.set_edgecolor('#c0562f')
            spine.set_linewidth(1.3)

# color the bottom-row x tick labels to match Sun's Fig 3D group colors
bottom_ax = g.axes.flat[-1]
bottom_ax.set_xticklabels(list(GROUPS.keys()))
for tick, group in zip(bottom_ax.get_xticklabels(), GROUPS.keys()):
    tick.set_color(GROUP_COLOR[group])
    tick.set_fontweight('bold')

g.fig.suptitle('Fate probability by state-group (Sun et al. Fig. 3D grouping), row-split by HSC subtype (c0–c9)\nS20 (BCG donor) — rows for c2/c3/c4/c7 outlined in orange',
               fontsize=12.5, fontweight='bold', y=1.01)
g.fig.tight_layout()
g.savefig(f'{OUT}/fate_grouped_rowsplit_S20.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')
