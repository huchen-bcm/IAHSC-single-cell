import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

# single global scale used identically across every panel, every threshold:
# dot area (pt^2) = N_GENES_TO_AREA * (number of overlapping genes)
# (this constant is fixed across ALL panels/thresholds below, so a dot of a given
# size always means the same gene count no matter which panel you're looking at;
# only the LEGEND reference values shown are picked per-panel to match that
# panel's own data range, so the legend is actually readable)
N_GENES_TO_AREA = 15

DATASETS = [
    ('r0', 'mouse r>0 AND human r>0, n=1746 genes',
     'integration_driver_GO_positive_IA_c347.csv', 'integration_driver_Reactome_positive_IA_c347.csv'),
    ('r01', 'mouse r>0.1 AND human r>0.1, n=305 genes',
     'integration_driver_GO_positive_r01.csv', 'integration_driver_Reactome_positive_r01.csv'),
    ('r02', 'mouse r>0.2 AND human r>0.2, n=12 genes',
     'integration_driver_GO_positive_r02.csv', 'integration_driver_Reactome_positive_r02.csv'),
]


def load(fname, n=15):
    df = pd.read_csv(f'{OUT}/{fname}').sort_values('Adjusted P-value').head(n).copy()
    df['nlp'] = -np.log10(df['Adjusted P-value'])
    df['n_genes'] = df['Overlap'].apply(lambda s: int(s.split('/')[0]))
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


def add_size_legend(ax, color, n_genes):
    for n in pick_legend_values(n_genes):
        ax.scatter([], [], s=n * N_GENES_TO_AREA, c=color, alpha=0.8, edgecolor='black',
                   linewidth=0.5, label=f'{n} genes')
    ax.legend(title='gene overlap', loc='lower right', fontsize=8, title_fontsize=8.5)


for tag, subtitle, go_file, re_file in DATASETS:
    go = load(go_file)
    re = load(re_file)

    fig, axes = plt.subplots(1, 2, figsize=(17, 6.5))
    for ax, df, color, title in [
        (axes[0], go, '#4a6fa5', 'GO Biological Process'),
        (axes[1], re, '#5a9e6f', 'Reactome'),
    ]:
        ax.scatter(df['nlp'], range(len(df)), s=df['n_genes'] * N_GENES_TO_AREA, c=color,
                   alpha=0.8, edgecolor='black', linewidth=0.5)
        ax.set_yticks(range(len(df)))
        ax.set_yticklabels(df['Term'], fontsize=8.5)
        ax.axvline(-np.log10(0.05), color='grey', linestyle='--', lw=0.8)
        ax.set_xlabel('-log10(FDR)')
        ax.set_title(title, fontsize=12, fontweight='bold')
        add_size_legend(ax, color, df['n_genes'])

    fig.suptitle(f'Enrichment of genes concordantly correlated with the IA-HSC/c3-4-7 axis\n({subtitle})',
                 fontsize=13, fontweight='bold', y=1.04)
    plt.tight_layout()
    plt.savefig(f'{OUT}/assets/integration_driver_{tag}_go_reactome.pdf', facecolor='white', bbox_inches='tight')
    plt.savefig(f'{OUT}/assets/integration_driver_{tag}_go_reactome.png', dpi=200, facecolor='white', bbox_inches='tight')
    print(f'saved {tag}')
