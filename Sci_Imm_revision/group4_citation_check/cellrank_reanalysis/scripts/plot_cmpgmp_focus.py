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

GROUPS = {
    'CMPs/GMPs': ['CMP_a', 'CMP_b', 'CMP_c', 'GMP_a', 'GMP_b'],
    'MLPs/PreBNK': ['MLP_a', 'MLP_b', 'PreBNK'],
    'MEPs': ['MEP_a', 'MEP_b', 'MEP_c'],
}
for g, cols in GROUPS.items():
    df[g] = df[cols].sum(axis=1)

medians_for_order = df.groupby('hsc_clusters_lab', observed=True)['CMPs/GMPs'].median().sort_values(ascending=False)
c_order = medians_for_order.index.tolist()
df['hsc_clusters_lab'] = pd.Categorical(df['hsc_clusters_lab'], categories=c_order, ordered=True)

fig, ax = plt.subplots(figsize=(9, 5))
colors = ['#c0562f' if c in ('c2', 'c3', 'c4', 'c7') else '#9aa0a6' for c in c_order]

sns.boxplot(data=df, x='hsc_clusters_lab', y='CMPs/GMPs', order=c_order, ax=ax,
            palette=colors, fliersize=1.5, linewidth=1, showfliers=False)

# overlay group medians as text + a reference line for the overall median across all clusters
overall_median = df['CMPs/GMPs'].median()
ax.axhline(overall_median, color='black', linestyle='--', lw=0.8, alpha=0.6, zorder=0)
ax.text(9.6, overall_median + 0.02, f'all-cluster median = {overall_median:.2f}', fontsize=8, color='#333333', ha='right')

medians = df.groupby('hsc_clusters_lab', observed=True)['CMPs/GMPs'].median().reindex(c_order)
for i, (c, m) in enumerate(medians.items()):
    ax.text(i, 1.03, f'{m:.2f}', ha='center', fontsize=9, fontweight='bold' if c in ('c2', 'c3', 'c4', 'c7') else 'normal',
            color='#c0562f' if c in ('c2', 'c3', 'c4', 'c7') else '#555555')

for xt, c in zip(ax.get_xticklabels(), c_order):
    if c in ('c2', 'c3', 'c4', 'c7'):
        xt.set_color('#c0562f')
        xt.set_fontweight('bold')

ax.set_ylim(0, 1.1)
ax.set_ylabel('CMPs/GMPs fate probability (per cell)')
ax.set_xlabel('HSC subcluster')
ax.set_title('CMPs/GMPs fate bias by HSC subcluster, both timepoints pooled\nS20 (BCG donor) — c2/c3/c4/c7 in orange, median labeled above each box',
              fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig(f'{OUT}/cmpgmp_focus_S20.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')

print(medians.round(3))
