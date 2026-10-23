import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
res = pd.read_csv(f'{OUT}/gsea_prerank_hallmark_crossspecies.csv')
res = res.sort_values('NES', ascending=False).reset_index(drop=True)

fig, ax = plt.subplots(figsize=(8, 11))
colors = ['#c0562f' if q < 0.25 else '#9aa0a6' for q in res['FDR q-val']]
ax.barh(range(len(res)), res['NES'], color=colors, edgecolor='black', linewidth=0.4)
ax.set_yticks(range(len(res)))
ax.set_yticklabels(res['Term'], fontsize=9)
ax.invert_yaxis()
ax.axvline(0, color='black', lw=0.8)
ax.set_xlabel('Normalized Enrichment Score (NES)')
ax.set_title('Pre-ranked GSEA: MSigDB Hallmark gene sets\nranked by mouse_corr x human_corr (cross-species axis concordance)\n'
             'colored: FDR<0.25 (standard GSEA threshold)', fontsize=11, fontweight='bold')

handles = [plt.Rectangle((0, 0), 1, 1, color='#c0562f', label='FDR < 0.25'),
           plt.Rectangle((0, 0), 1, 1, color='#9aa0a6', label='not significant')]
ax.legend(handles=handles, loc='lower right', fontsize=9)

plt.tight_layout()
plt.savefig(f'{OUT}/assets/gsea_prerank_hallmark_crossspecies.pdf', facecolor='white', bbox_inches='tight')
plt.savefig(f'{OUT}/assets/gsea_prerank_hallmark_crossspecies.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')
