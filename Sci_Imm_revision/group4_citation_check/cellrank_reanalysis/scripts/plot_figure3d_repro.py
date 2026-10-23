import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results'
META = '/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/GEO_data/sc_pdata_phate_ndim2_knn6_t10_gamma1.csv'
meta = pd.read_csv(META, index_col=0)

GROUPS = {
    'CMPs/GMPs': ['CMP_a', 'CMP_b', 'CMP_c', 'GMP_a', 'GMP_b'],
    'MLPs/PreBNK': ['MLP_a', 'MLP_b', 'PreBNK'],
    'MEPs': ['MEP_a', 'MEP_b', 'MEP_c'],
}
GROUP_COLOR = {'CMPs/GMPs': '#3a8f4a', 'MLPs/PreBNK': '#3a6fa8', 'MEPs': '#d98a2b'}

samples = [
    ('CTL', 'Td0', 'S18_Td0_CTL'),
    ('CTL', 'Tm3', 'S18_Tm3_CTL'),
    ('BCG', 'Td0', 'S20_Td0_BCG'),
    ('BCG', 'Tm3', 'S20_Tm3_BCG'),
]

fig, axes = plt.subplots(2, 2, figsize=(10, 10))
row_map = {'CTL': 0, 'BCG': 1}
col_map = {'Td0': 0, 'Tm3': 1}

for vacc, tp, label in samples:
    ax = axes[row_map[vacc], col_map[tp]]
    tpdf = pd.read_csv(f'{OUT}/{label}_to_terminal_states.csv')
    tpdf = tpdf.set_index('Cell ID')
    for g, cols in GROUPS.items():
        tpdf[g] = tpdf[cols].sum(axis=1)
    tpdf['predicted_fate'] = tpdf[list(GROUPS.keys())].idxmax(axis=1)

    m = meta.loc[meta.index.intersection(tpdf.index), ['hsc_PHATE1', 'hsc_PHATE2']]
    plotdf = m.join(tpdf['predicted_fate'])

    for g, color in GROUP_COLOR.items():
        sub = plotdf[plotdf['predicted_fate'] == g]
        ax.scatter(sub['hsc_PHATE1'], sub['hsc_PHATE2'], s=6, color=color, alpha=0.75, linewidth=0, label=g)

    ax.set_title(f'{vacc} — {tp}  (n={len(plotdf)}, donor {label.split("_")[0]})', fontsize=11, fontweight='bold')
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)

handles = [plt.Line2D([0], [0], marker='o', color='w', markerfacecolor=c, markersize=8, label=g)
           for g, c in GROUP_COLOR.items()]
fig.legend(handles=handles, loc='lower center', ncol=3, bbox_to_anchor=(0.5, -0.02), fontsize=10, frameon=False)

fig.suptitle('Reproducing Fig. 3D: predicted terminal fate on the HSC PHATE map\n(rows: placebo/CTL vs. BCG; columns: Td0 vs. Tm3 — one donor per group)',
             fontsize=13, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(f'{OUT}/figure3d_repro.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')

for vacc, tp, label in samples:
    tpdf = pd.read_csv(f'{OUT}/{label}_to_terminal_states.csv').set_index('Cell ID')
    for g, cols in GROUPS.items():
        tpdf[g] = tpdf[cols].sum(axis=1)
    tpdf['predicted_fate'] = tpdf[list(GROUPS.keys())].idxmax(axis=1)
    print(label, tpdf['predicted_fate'].value_counts(normalize=True).round(3).to_dict())
