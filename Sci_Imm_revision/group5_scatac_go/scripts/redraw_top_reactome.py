import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group5_scatac_go'

def top_reactome(path, n=8):
    df = pd.read_csv(path)
    df = df[df['Gene_set'] == 'Reactome_Pathways_2024'].copy()
    df = df.sort_values('Adjusted P-value').head(n)
    df['nlp'] = -np.log10(df['Adjusted P-value'])
    return df.iloc[::-1]  # reverse for barh (largest on top)

d1 = top_reactome(f'{OUT}/enrichr.IAHSC_vs_OtherHSC.opening_in_IAHSC.csv')
d2 = top_reactome(f'{OUT}/enrichr.Infected_vs_Naive.opening_in_Infected.csv')

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))
thresh = -np.log10(0.05)

axes[0].barh(range(len(d1)), d1['nlp'], color='#4b3f9e')
axes[0].set_yticks(range(len(d1)))
axes[0].set_yticklabels(d1['Term'])
axes[0].axvline(thresh, color='grey', linestyle='--', lw=1)
axes[0].set_xlabel('-log10(FDR)')
axes[0].set_title('IA-HSC vs OtherHSC\ngenes near opening peaks', fontsize=12)

axes[1].barh(range(len(d2)), d2['nlp'], color='#d9564a')
axes[1].set_yticks(range(len(d2)))
axes[1].set_yticklabels(d2['Term'])
axes[1].axvline(thresh, color='grey', linestyle='--', lw=1)
axes[1].set_xlabel('-log10(FDR)')
axes[1].set_title('Infected vs Naive (IA-HSC)\ngenes near opening peaks', fontsize=12)

fig.suptitle('Top Reactome pathway terms (genes near opening peaks only)', fontsize=13, color='#555555', y=1.03)
plt.tight_layout()
plt.savefig(f'{OUT}/assets/top_reactome_terms_opening_redrawn.pdf', facecolor='white', bbox_inches='tight')
plt.savefig(f'{OUT}/assets/top_reactome_terms_opening_redrawn.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')
