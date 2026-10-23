import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results'
HSC_META = '/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/GEO_data/sc_pdata_phate_ndim2_knn6_t10_gamma1.csv'
hsc = pd.read_csv(HSC_META, index_col=0)

GROUPS = {
    'CMPs/GMPs': ['CMP_a','CMP_b','CMP_c','GMP_a','GMP_b'],
    'MLPs/PreBNK': ['MLP_a','MLP_b','PreBNK'],
    'MEPs': ['MEP_a','MEP_b','MEP_c'],
}
COLORS = {'CMPs/GMPs':'#2ca02c', 'MLPs/PreBNK':'#1f3a93', 'MEPs':'#e08214'}

# rows = CTL / BCG (both at D90/Tm3); columns = two different example donors per row
panels = [
    ('CTL', 'S11_Tm3_CTL', 'S11'),
    ('CTL', 'S18_Tm3_CTL', 'S18'),
    ('BCG', 'S1_Tm3_BCG',  'S1'),
    ('BCG', 'S20_Tm3_BCG', 'S20'),
]

def load(label):
    df = pd.read_csv(f'{OUT}/{label}_to_terminal_states.csv').set_index('Cell ID')
    for g, cols in GROUPS.items():
        cols_present = [c for c in cols if c in df.columns]
        df[g] = df[cols_present].sum(axis=1) if cols_present else 0.0
    df['pred_group'] = df[list(GROUPS.keys())].idxmax(axis=1)
    j = df.join(hsc[['hsc_PHATE1','hsc_PHATE2','clust_names']], how='inner')
    return j[j['clust_names'].isin(['HSC_a','HSC_b'])]

fig, axes = plt.subplots(2, 2, figsize=(9, 8.5))
row_order = ['CTL', 'BCG']
for row_i, grp in enumerate(row_order):
    cols = [p for p in panels if p[0] == grp]
    for col_i, (_, label, donor) in enumerate(cols):
        ax = axes[row_i, col_i]
        j = load(label)
        for g, c in COLORS.items():
            sub = j[j['pred_group'] == g]
            ax.scatter(sub['hsc_PHATE1'], sub['hsc_PHATE2'], s=6, c=c, alpha=0.75, linewidths=0)
        ax.set_title(f'{donor} ({grp}) D90, n={len(j)}', fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
    axes[row_i, 0].set_ylabel(f'{grp} D90', fontsize=13, fontweight='bold')

handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=8, label=g) for g,c in COLORS.items()]
fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.02))
plt.tight_layout()
plt.savefig(f'{OUT}/figure3d_repro_v2.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved figure3d_repro_v2.png')
