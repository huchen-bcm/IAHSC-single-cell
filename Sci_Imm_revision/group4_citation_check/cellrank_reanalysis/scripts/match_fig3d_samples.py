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

donor_group = {
    'S1':'BCG','S5':'BCG','S6':'BCG','S7':'CTL','S8':'BCG','S10':'BCG',
    'S11':'CTL','S13':'CTL','S14':'BCG','S15':'BCG','S16':'BCG','S17':'BCG',
    'S18':'CTL','S20':'BCG','S21':'BCG',
}

samples = []
for donor, grp in donor_group.items():
    label = f'{donor}_Tm3_{grp}'
    try:
        df = pd.read_csv(f'{OUT}/{label}_to_terminal_states.csv').set_index('Cell ID')
    except FileNotFoundError:
        continue
    for g, cols in GROUPS.items():
        cols_present = [c for c in cols if c in df.columns]
        df[g] = df[cols_present].sum(axis=1) if cols_present else 0.0
    df['pred_group'] = df[list(GROUPS.keys())].idxmax(axis=1)
    joined = df.join(hsc[['hsc_PHATE1','hsc_PHATE2','clust_names']], how='inner')
    joined = joined[joined['clust_names'].isin(['HSC_a','HSC_b'])]
    samples.append((donor, grp, joined))

samples.sort(key=lambda x: (x[1], x[0]))

n = len(samples)
ncols = 4
nrows = int(np.ceil(n/ncols))
fig, axes = plt.subplots(nrows, ncols, figsize=(4*ncols, 3.6*nrows))
axes = axes.flatten()
for i, (donor, grp, joined) in enumerate(samples):
    ax = axes[i]
    for g, c in COLORS.items():
        sub = joined[joined['pred_group']==g]
        ax.scatter(sub['hsc_PHATE1'], sub['hsc_PHATE2'], s=6, c=c, alpha=0.75, linewidths=0)
    frac = (joined['pred_group'].value_counts(normalize=True)*100).round(1)
    fr_str = ', '.join(f"{g}:{frac.get(g,0):.0f}%" for g in COLORS)
    ax.set_title(f'{donor} ({grp}) D90, n={len(joined)}\n{fr_str}', fontsize=8.5)
    ax.set_xticks([]); ax.set_yticks([])
for j in range(i+1, len(axes)):
    axes[j].axis('off')
plt.tight_layout()
plt.savefig(f'{OUT}/fig3d_all_d90_candidates.png', dpi=170, facecolor='white', bbox_inches='tight')

summary = pd.DataFrame([
    {'donor': d, 'group': g, 'n_hsc': len(j),
     **{grp: round((j['pred_group']==grp).mean()*100,1) for grp in COLORS}}
    for d, g, j in samples
])
print(summary.to_string(index=False))
summary.to_csv(f'{OUT}/fig3d_all_d90_candidates_summary.csv', index=False)
