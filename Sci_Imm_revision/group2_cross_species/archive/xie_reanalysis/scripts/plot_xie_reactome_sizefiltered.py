import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

XIE = '/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/Xie'
TERM_SIZE_MAX = 100
N_GENES_TO_AREA = 15

DATASETS = [
    ('r0', 'r>0 (n=531)', 'xie_reactome_r0.csv'),
    ('r005', 'r>0.05 (n=80)', 'xie_reactome_r005.csv'),
]


def load(fname, n=20):
    df = pd.read_csv(f'{XIE}/{fname}')
    df['n_genes'] = df['Overlap'].apply(lambda s: int(s.split('/')[0]))
    df['term_size'] = df['Overlap'].apply(lambda s: int(s.split('/')[1]))
    df = df[(df['term_size'] < TERM_SIZE_MAX) & (df['Adjusted P-value'] < 0.05)]
    df = df.sort_values('Adjusted P-value').head(n).copy()
    df['nlp'] = -np.log10(df['Adjusted P-value'])
    return df.iloc[::-1]


def pick_legend_values(n_genes):
    lo, hi = int(n_genes.min()), int(n_genes.max())
    if hi <= 5:
        cand = [1, 2, 3, 4, 5]
    elif hi <= 20:
        cand = [1, 5, 10, 20]
    elif hi <= 100:
        cand = [5, 25, 50, 100]
    else:
        cand = [5, 25, 100, 300]
    return [v for v in cand if lo <= v <= hi] or [lo, hi]


fig, axes = plt.subplots(1, 2, figsize=(14, 6.5))
color = '#5a9e6f'

for ax, (tag, label, fname) in zip(axes, DATASETS):
    df = load(fname)
    if len(df) == 0:
        ax.text(0.5, 0.5, 'no terms survive filter', ha='center', va='center', transform=ax.transAxes)
        ax.set_title(label, fontsize=12, fontweight='bold')
        continue
    ax.scatter(df['nlp'], range(len(df)), s=df['n_genes'] * N_GENES_TO_AREA, c=color,
               alpha=0.8, edgecolor='black', linewidth=0.5)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df['Term'], fontsize=8.5)
    ax.axvline(-np.log10(0.05), color='grey', linestyle='--', lw=0.8)
    ax.set_xlabel('-log10(FDR)')
    ax.set_title(label, fontsize=12, fontweight='bold')
    for v in pick_legend_values(df['n_genes']):
        ax.scatter([], [], s=v * N_GENES_TO_AREA, c=color, alpha=0.8, edgecolor='black',
                   linewidth=0.5, label=f'{v} genes')
    ax.legend(title='gene overlap', loc='lower right', fontsize=7.5, title_fontsize=8)

fig.suptitle(f'Xie et al. integration: Reactome enrichment, genes concordantly correlated with the\n'
             f'mouse-derived IA-HSC axis (terms with Reactome pathway size < {TERM_SIZE_MAX} genes only)',
             fontsize=13, fontweight='bold', y=1.05)
plt.tight_layout()
plt.savefig(f'{XIE}/xie_reactome_sizefiltered.pdf', facecolor='white', bbox_inches='tight')
plt.savefig(f'{XIE}/xie_reactome_sizefiltered.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')
