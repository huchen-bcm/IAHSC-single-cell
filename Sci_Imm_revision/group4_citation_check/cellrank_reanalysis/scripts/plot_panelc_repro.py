import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results'
META = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/FINAL_phates_meta.csv'
meta = pd.read_csv(META, index_col=0)

rows = []
for label, donor, group, tp in [
    ('S20_Td0_BCG', 'S20', 'BCG', 'Td0'),
    ('S20_Tm3_BCG', 'S20', 'BCG', 'Tm3'),
    ('S18_Td0_CTL', 'S18', 'CTL', 'Td0'),
    ('S18_Tm3_CTL', 'S18', 'CTL', 'Tm3'),
]:
    tpdf = pd.read_csv(f'{OUT}/{label}_to_terminal_states.csv').set_index('Cell ID')
    most_likely = tpdf.idxmax(axis=1)
    hsc_mask = meta.loc[meta.index.intersection(tpdf.index), 'clust_names'].isin(['HSC_a', 'HSC_b'])
    hsc_cells = hsc_mask[hsc_mask].index
    ml_hsc = most_likely.loc[most_likely.index.intersection(hsc_cells)]
    pct = (ml_hsc == 'CMP_b').mean() * 100
    rows.append((donor, group, tp, pct))

df = pd.DataFrame(rows, columns=['donor', 'group', 'timepoint', 'pct_cmpb'])
d90 = df[df['timepoint'] == 'Tm3']

fig, ax = plt.subplots(figsize=(3.2, 4))
colors = {'BCG': '#FFD320', 'CTL': '#004586'}
for _, row in d90.iterrows():
    ax.scatter(row['group'], row['pct_cmpb'], s=140, color=colors[row['group']], edgecolor='black', zorder=3)
    ax.annotate(row['donor'], (row['group'], row['pct_cmpb']), textcoords='offset points', xytext=(10, 0), fontsize=9)

ax.set_ylabel('% CMP_b (= CMP c2) biased HSCs')
ax.set_title('CMP_b bias at Day 90\n(1 donor per group so far)', fontsize=11, fontweight='bold')
ax.set_xlim(-0.5, 1.5)
ax.set_ylim(0, max(d90['pct_cmpb']) * 1.3)
plt.tight_layout()
plt.savefig(f'{OUT}/panelc_repro_d90.png', dpi=200, facecolor='white', bbox_inches='tight')
print(df)
print('saved')
