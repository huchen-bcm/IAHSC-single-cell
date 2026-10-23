import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu

META = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/FINAL_phates_meta.csv'
meta = pd.read_csv(META, index_col=0)

donor_group = {
    'S1': 'BCG', 'S5': 'BCG', 'S6': 'BCG', 'S7': 'CTL', 'S8': 'BCG', 'S10': 'BCG',
    'S11': 'CTL', 'S13': 'CTL', 'S14': 'BCG', 'S15': 'BCG', 'S16': 'BCG', 'S17': 'BCG',
    'S18': 'CTL', 'S20': 'BCG', 'S21': 'BCG',
}

def pct_cmpb(path, meta):
    tpdf = pd.read_csv(path).set_index('Cell ID')
    most_likely = tpdf.idxmax(axis=1)
    hsc_mask = meta.loc[meta.index.intersection(tpdf.index), 'clust_names'].isin(['HSC_a', 'HSC_b'])
    hsc_cells = hsc_mask[hsc_mask].index
    ml_hsc = most_likely.loc[most_likely.index.intersection(hsc_cells)]
    return (ml_hsc == 'CMP_b').mean() * 100, len(ml_hsc)

for version, OUT in [('CellRank 2.x', '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results'),
                      ('CellRank 1.x', '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results_v1')]:
    rows = []
    for donor, group in donor_group.items():
        try:
            pct_d0, n0 = pct_cmpb(f'{OUT}/{donor}_Td0_{group}_to_terminal_states.csv', meta)
            pct_d90, n90 = pct_cmpb(f'{OUT}/{donor}_Tm3_{group}_to_terminal_states.csv', meta)
        except FileNotFoundError as e:
            print('MISSING', e)
            continue
        shift = pct_d90 - pct_d0
        rows.append((donor, group, pct_d0, pct_d90, shift, n0, n90))
    df = pd.DataFrame(rows, columns=['donor','group','pct_d0','pct_d90','shift','n0','n90'])
    print(f'\n=== {version} ===')
    print(df.to_string(index=False))
    bcg = df[df['group']=='BCG']['shift'].values
    ctl = df[df['group']=='CTL']['shift'].values
    stat, p = mannwhitneyu(bcg, ctl, alternative='two-sided')
    print(f'shift score: BCG median={np.median(bcg):.2f} (n={len(bcg)}), CTL median={np.median(ctl):.2f} (n={len(ctl)}), Mann-Whitney p={p:.4f}')
