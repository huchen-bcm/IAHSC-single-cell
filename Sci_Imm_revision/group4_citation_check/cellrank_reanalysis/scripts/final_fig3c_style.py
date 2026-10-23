import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT2 = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results'
OUT1 = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results_v1'
META = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/FINAL_phates_meta.csv'
meta = pd.read_csv(META, index_col=0)

donor_group = {
    'S1': 'BCG', 'S5': 'BCG', 'S6': 'BCG', 'S7': 'CTL', 'S8': 'BCG', 'S10': 'BCG',
    'S11': 'CTL', 'S13': 'CTL', 'S14': 'BCG', 'S15': 'BCG', 'S16': 'BCG', 'S17': 'BCG',
    'S18': 'CTL', 'S20': 'BCG', 'S21': 'BCG',
}

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

fig, axes = plt.subplots(1, 2, figsize=(6.4, 4.2), sharey=False)
colors = {'BCG': '#FFD320', 'CTL': '#004586'}

for ax, (title, OUT) in zip(axes, [('CellRank 2.0.6\n(modern)', OUT2), ('CellRank 1.5.1\n(period-matched)', OUT1)]):
    df = get_d90(OUT)
    bcg = df[df['group']=='BCG']['pct'].values
    ctl = df[df['group']=='CTL']['pct'].values
    _, p = mannwhitneyu(bcg, ctl, alternative='two-sided')

    box = ax.boxplot([ctl, bcg], positions=[0,1], widths=0.5, patch_artist=True, showfliers=False)
    for patch, grp in zip(box['boxes'], ['CTL','BCG']):
        patch.set_facecolor(colors[grp]); patch.set_alpha(0.55); patch.set_edgecolor('black')
    for med in box['medians']:
        med.set_color('black')
    np.random.seed(0)
    for i, vals in enumerate([ctl, bcg]):
        jitter = np.random.uniform(-0.12, 0.12, size=len(vals))
        grp = ['CTL','BCG'][i]
        ax.scatter(np.full(len(vals), i) + jitter, vals, color=colors[grp], s=45,
                   edgecolor='black', linewidth=0.6, zorder=3)
    ax.set_xticks([0,1]); ax.set_xticklabels(['Placebo', 'BCG'])
    ax.set_ylabel('% HSCs with CMP c2 bias (D90)')
    sig = f'p = {p:.3f}' + (' (paper: p=0.018*)' if False else '')
    ax.set_title(f'{title}\n{sig}', fontsize=10.5)
    ax.set_ylim(bottom=-1)

plt.suptitle('Reproduction of Figure 3C (CMP c2 fate bias, BCG vs. placebo, D90)', fontsize=11.5, y=1.02)
plt.tight_layout()
plt.savefig(f'{OUT2}/figure3c_final_both_versions.png', dpi=220, facecolor='white', bbox_inches='tight')
print('saved')
