import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT2 = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results'
META = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/FINAL_phates_meta.csv'
HSC_META = '/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/GEO_data/sc_pdata_phate_ndim2_knn6_t10_gamma1.csv'
meta = pd.read_csv(META, index_col=0)
hsc = pd.read_csv(HSC_META, index_col=0)

donor_group = {
    'S1': 'BCG', 'S5': 'BCG', 'S6': 'BCG', 'S7': 'CTL', 'S8': 'BCG', 'S10': 'BCG',
    'S11': 'CTL', 'S13': 'CTL', 'S14': 'BCG', 'S15': 'BCG', 'S16': 'BCG', 'S17': 'BCG',
    'S18': 'CTL', 'S20': 'BCG', 'S21': 'BCG',
}

GROUPS = {
    'CMPs/GMPs': ['CMP_a','CMP_b','CMP_c','GMP_a','GMP_b'],
    'MLPs/PreBNK': ['MLP_a','MLP_b','PreBNK'],
    'MEPs': ['MEP_a','MEP_b','MEP_c'],
}
FATE_COLORS = {'CMPs/GMPs':'#2ca02c', 'MLPs/PreBNK':'#1f3a93', 'MEPs':'#e08214'}

# ============ FIGURE 1: clean Figure 3C style, BCG left, CTL right, no p-value ============
def get_d90(OUT):
    rows = []
    for donor, group in donor_group.items():
        label = f'{donor}_Tm3_{group}'
        tpdf = pd.read_csv(f'{OUT}/{label}_to_terminal_states.csv').set_index('Cell ID')
        most_likely = tpdf.idxmax(axis=1)
        hsc_mask = meta.loc[meta.index.intersection(tpdf.index), 'clust_names'].isin(['HSC_a', 'HSC_b'])
        hsc_cells = hsc_mask[hsc_mask].index
        ml_hsc = most_likely.loc[most_likely.index.intersection(hsc_cells)]
        pct = (ml_hsc == 'CMP_b').mean() * 100
        rows.append((donor, group, pct))
    return pd.DataFrame(rows, columns=['donor','group','pct'])

df = get_d90(OUT2)
bcg = df[df['group']=='BCG']['pct'].values
ctl = df[df['group']=='CTL']['pct'].values

fig, ax = plt.subplots(figsize=(3.4, 4.4))
colors = {'BCG': '#FFD320', 'CTL': '#004586'}
box = ax.boxplot([bcg, ctl], positions=[0,1], widths=0.5, patch_artist=True, showfliers=False)
for patch, grp in zip(box['boxes'], ['BCG','CTL']):
    patch.set_facecolor(colors[grp]); patch.set_alpha(0.55); patch.set_edgecolor('black')
for med in box['medians']:
    med.set_color('black')
np.random.seed(0)
for i, vals in enumerate([bcg, ctl]):
    jitter = np.random.uniform(-0.12, 0.12, size=len(vals))
    grp = ['BCG','CTL'][i]
    ax.scatter(np.full(len(vals), i) + jitter, vals, color=colors[grp], s=45,
               edgecolor='black', linewidth=0.6, zorder=3)
ax.set_xticks([0,1]); ax.set_xticklabels(['BCG', 'CTL'])
ax.set_ylabel('% HSCs with CMP c2 bias (D90)')
ax.set_ylim(bottom=-1)
plt.tight_layout()
plt.savefig(f'{OUT2}/figure3c_clean.pdf', dpi=220, facecolor='white', bbox_inches='tight')
print('saved figure3c_clean.pdf')

# ============ FIGURE 2: Figure 3D-style grid, row1=CTL D90 (all 4), row2=BCG D90 (5 picked) ============
ctl_donors = ['S7','S11','S13','S18']
bcg_donors = ['S1','S5','S6','S8','S10']

def load_phate(label):
    d = pd.read_csv(f'{OUT2}/{label}_to_terminal_states.csv').set_index('Cell ID')
    for g, cols in GROUPS.items():
        cols_present = [c for c in cols if c in d.columns]
        d[g] = d[cols_present].sum(axis=1) if cols_present else 0.0
    d['pred_group'] = d[list(GROUPS.keys())].idxmax(axis=1)
    j = d.join(hsc[['hsc_PHATE1','hsc_PHATE2','clust_names']], how='inner')
    return j[j['clust_names'].isin(['HSC_a','HSC_b'])]

ncols = 5
fig, axes = plt.subplots(2, ncols, figsize=(3.1*ncols, 6.6))
for row_i, (grp_name, donors) in enumerate([('CTL D90', ctl_donors), ('BCG D90', bcg_donors)]):
    for col_i in range(ncols):
        ax = axes[row_i, col_i]
        if col_i >= len(donors):
            ax.axis('off'); continue
        donor = donors[col_i]
        label = f'{donor}_Tm3_{donor_group[donor]}'
        j = load_phate(label)
        for g, c in FATE_COLORS.items():
            sub = j[j['pred_group']==g]
            ax.scatter(sub['hsc_PHATE1'], sub['hsc_PHATE2'], s=6, c=c, alpha=0.75, linewidths=0)
        ax.set_title(donor, fontsize=11)
        ax.set_xticks([]); ax.set_yticks([])
    axes[row_i, 0].set_ylabel(grp_name, fontsize=13, fontweight='bold')

handles = [plt.Line2D([0],[0], marker='o', color='w', markerfacecolor=c, markersize=8, label=g) for g,c in FATE_COLORS.items()]
fig.legend(handles=handles, loc='lower center', ncol=3, frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.02))
plt.tight_layout()
plt.savefig(f'{OUT2}/figure3d_5per_group.pdf', dpi=180, facecolor='white', bbox_inches='tight')
print('saved figure3d_5per_group.pdf')

# ============ FIGURE 3: P(green=CMPs/GMPs) per cell, boxplot by HSC subtype c0-c9 ============
GREEN_COLS = ['CMP_a','CMP_b','CMP_c','GMP_a','GMP_b']
all_rows = []
for donor, group in donor_group.items():
    for tp in ['Td0','Tm3']:
        label = f'{donor}_{tp}_{group}'
        try:
            d = pd.read_csv(f'{OUT2}/{label}_to_terminal_states.csv').set_index('Cell ID')
        except FileNotFoundError:
            continue
        cols_present = [c for c in GREEN_COLS if c in d.columns]
        green_prob = d[cols_present].sum(axis=1) if cols_present else pd.Series(0.0, index=d.index)
        j = pd.DataFrame({'green_prob': green_prob}).join(hsc[['hsc_clusters_lab','clust_names']], how='inner')
        j = j[j['clust_names'].isin(['HSC_a','HSC_b'])]
        all_rows.append(j)

pooled = pd.concat(all_rows)
all_subtypes = pooled['hsc_clusters_lab'].dropna().unique()
subtypes = sorted(all_subtypes, key=lambda s: -np.median(pooled.loc[pooled['hsc_clusters_lab']==s, 'green_prob']))
data = [pooled.loc[pooled['hsc_clusters_lab']==s, 'green_prob'].values for s in subtypes]
medians = [np.median(x) for x in data]

fig, ax = plt.subplots(figsize=(8, 4.5))
box = ax.boxplot(data, positions=range(len(subtypes)), widths=0.6, patch_artist=True, showfliers=False)
for patch in box['boxes']:
    patch.set_facecolor('#8fbf7f'); patch.set_alpha(0.7); patch.set_edgecolor('black')
for med in box['medians']:
    med.set_color('black')
ax.axhline(np.median(pooled['green_prob']), color='grey', linestyle='--', linewidth=1)
for i, m in enumerate(medians):
    ax.text(i, m+0.03, f'{m:.2f}', ha='center', fontsize=8.5)
ax.set_xticks(range(len(subtypes))); ax.set_xticklabels(subtypes)
ax.set_xlabel('HSC subtype')
ax.set_ylabel('P(CMP c1/c2/c3 + GMP c1/c2)\nper cell')
ax.set_ylim(0, 1.05)
plt.tight_layout()
plt.savefig(f'{OUT2}/hsc_subtype_green_bias_boxplot.pdf', dpi=200, facecolor='white', bbox_inches='tight')
print('saved hsc_subtype_green_bias_boxplot.pdf')
print(pooled.groupby('hsc_clusters_lab')['green_prob'].agg(['median','mean','count']))
