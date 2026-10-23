import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results'
META = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/FINAL_phates_meta.csv'
meta = pd.read_csv(META, index_col=0)

# matches Sun et al.'s own donor list & group assignment (S2/S3/S12 already excluded)
donor_group = {
    'S1': 'BCG', 'S5': 'BCG', 'S6': 'BCG', 'S7': 'CTL', 'S8': 'BCG', 'S10': 'BCG',
    'S11': 'CTL', 'S13': 'CTL', 'S14': 'BCG', 'S15': 'BCG', 'S16': 'BCG', 'S17': 'BCG',
    'S18': 'CTL', 'S20': 'BCG', 'S21': 'BCG',
}

rows = []
for donor, group in donor_group.items():
    for tp in ['Td0', 'Tm3']:
        label = f'{donor}_{tp}_{group}'
        try:
            tpdf = pd.read_csv(f'{OUT}/{label}_to_terminal_states.csv').set_index('Cell ID')
        except FileNotFoundError:
            print(f'MISSING: {label}')
            continue
        most_likely = tpdf.idxmax(axis=1)
        hsc_mask = meta.loc[meta.index.intersection(tpdf.index), 'clust_names'].isin(['HSC_a', 'HSC_b'])
        hsc_cells = hsc_mask[hsc_mask].index
        ml_hsc = most_likely.loc[most_likely.index.intersection(hsc_cells)]
        pct = (ml_hsc == 'CMP_b').mean() * 100
        rows.append((donor, group, tp, pct, len(ml_hsc)))

df = pd.DataFrame(rows, columns=['donor', 'group', 'timepoint', 'pct_cmpb', 'n_hsc'])
df.to_csv(f'{OUT}/figure3c_full_percell_summary.csv', index=False)
print(df.to_string())

d90 = df[df['timepoint'] == 'Tm3'].copy()
bcg = d90[d90['group'] == 'BCG']['pct_cmpb'].values
ctl = d90[d90['group'] == 'CTL']['pct_cmpb'].values
print(f'\nD90: BCG n={len(bcg)}, median={np.median(bcg):.2f}%, CTL n={len(ctl)}, median={np.median(ctl):.2f}%')

stat, pval = mannwhitneyu(bcg, ctl, alternative='two-sided')
print(f'Mann-Whitney (Wilcoxon rank-sum) p-value = {pval:.4f}')

fig, ax = plt.subplots(figsize=(3, 4))
box = ax.boxplot([bcg, ctl], positions=[0, 1], widths=0.5, patch_artist=True, showfliers=False)
colors = ['#FFD320', '#004586']
for patch, c in zip(box['boxes'], colors):
    patch.set_facecolor(c)
    patch.set_alpha(0.5)
np.random.seed(0)
for i, vals in enumerate([bcg, ctl]):
    jitter = np.random.uniform(-0.15, 0.15, size=len(vals))
    ax.scatter(np.full(len(vals), i) + jitter, vals, color='black', s=40, zorder=3,
               facecolor=colors[i], edgecolor='black')
ax.set_xticks([0, 1])
ax.set_xticklabels(['BCG', 'CTL'])
ax.set_ylabel('% CMP_b (=CMP c2) biased HSCs')
ax.set_title(f'CMP_b bias Day 90 (full donor set)\np = {pval:.3f}', fontsize=10.5, fontweight='bold')
ax.axhline(0, color='brown', lw=0.8)
plt.tight_layout()
plt.savefig(f'{OUT}/figure3c_full_repro.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved figure')
