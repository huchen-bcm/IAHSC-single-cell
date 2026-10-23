import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
w = pd.read_csv(f'{OUT}/gsea_hallmark_CCcorrected_withTirosh.csv').set_index('Term')
wo = pd.read_csv(f'{OUT}/gsea_hallmark_CCcorrected_noTirosh.csv').set_index('Term')

terms = w.sort_values('FDR q-val').index.tolist()
w = w.loc[terms]
wo = wo.reindex(terms)

y = np.arange(len(terms))
fig, ax = plt.subplots(figsize=(8.5, 11))

ax.barh(y + 0.2, w['NES'], height=0.38, color='#4a6fa5', edgecolor='black', linewidth=0.4, label='with Tirosh S/G2M genes (n=2526)')
ax.barh(y - 0.2, wo['NES'], height=0.38, color='#c0562f', edgecolor='black', linewidth=0.4, label='without Tirosh S/G2M genes (n=2502)')

for yi, term in zip(y, terms):
    if w.loc[term, 'FDR q-val'] < 0.05:
        ax.text(w.loc[term, 'NES'] + 0.03, yi + 0.2, '*', fontsize=13, va='center', fontweight='bold')
    if not np.isnan(wo.loc[term, 'FDR q-val']) and wo.loc[term, 'FDR q-val'] < 0.05:
        ax.text(wo.loc[term, 'NES'] + 0.03, yi - 0.2, '*', fontsize=13, va='center', fontweight='bold')

ax.set_yticks(y)
ax.set_yticklabels(terms, fontsize=9)
ax.invert_yaxis()
ax.axvline(0, color='black', lw=0.8)
ax.set_xlabel('Normalized Enrichment Score (NES)')
ax.set_title('Effect of removing the 95-gene Tirosh S/G2M panel\nbefore pre-ranked GSEA (Hallmark, cell-cycle-corrected axis)\n* = FDR<0.05',
             fontsize=11.5, fontweight='bold')
ax.legend(loc='lower right', fontsize=9)

plt.tight_layout()
plt.savefig(f'{OUT}/assets/gsea_tirosh_removal_comparison.pdf', facecolor='white', bbox_inches='tight')
plt.savefig(f'{OUT}/assets/gsea_tirosh_removal_comparison.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')
