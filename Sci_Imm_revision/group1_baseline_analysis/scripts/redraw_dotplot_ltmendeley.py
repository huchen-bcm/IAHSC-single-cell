import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis'
df = pd.read_csv(f'{OUT}/LTmendeley_vs_Fig4A_logFC_merged.csv')
df = df[(df['LT_Mendeley_logFC'].abs() < 10) & (df['Fig4A_IAvsLT_logFC'].abs() < 10)].copy()

x = df['LT_Mendeley_logFC'].values
y = df['Fig4A_IAvsLT_logFC'].values
from scipy.stats import spearmanr
rho, p = spearmanr(x, y)
slope, intercept = np.polyfit(x, y, 1)
r = np.corrcoef(x, y)[0, 1]

fig, ax = plt.subplots(figsize=(7, 7))
sig = df['sig_bulk'].astype(bool)
ax.scatter(x[~sig], y[~sig], s=6, color='#d9d5c9', alpha=0.6, linewidth=0, label='other', zorder=1)
ax.scatter(x[sig], y[sig], s=28, color='#6a51a3', alpha=0.75, edgecolor='none', label='Fig.4A identity gene', zorder=2)

xs = np.linspace(-10, 10, 100)
ax.plot(xs, slope * xs + intercept, '--', color='#d95f02', lw=1.8, label=f'OLS fit (slope={slope:.2f})', zorder=3)

ax.axhline(0, color='grey', lw=0.7, zorder=0)
ax.axvline(0, color='grey', lw=0.7, zorder=0)
ax.set_xlim(-10, 10)
ax.set_ylim(-10, 10)
ax.set_xlabel('log2FC: LT-HSC infected vs naive (Mendeley bulk)')
ax.set_ylabel('log2FC: IA-HSC vs LT-HSC, both infected (Fig.4A, zhijiany bulk)')
ax.set_title(f'n={len(df)} genes (|log2FC|<10 both axes)\nSpearman ρ={rho:.2f} (p={p:.1e}), OLS r={r:.2f}\nFig.4B leading-edge genes labeled', fontsize=11)

categories = {
    'Response to bacterium': ('#2ca08a', ['Plac8', 'Mpeg1', 'Tent5a', 'Mt2', 'Tifab']),
    'Myeloid cell activation': ('#e8836a', ['Myo1f', 'Ifngr1', 'Sphk1', 'Cd48', 'Cebpa', 'Cst7', 'Tyrobp', 'Rab44', 'Ms4a2', 'Fcgr3']),
    'Monocyte chemotaxis': ('#d4b106', ['Flt1', 'Ccl9', 'Ccr1', 'Ctsg', 'Cx3cr1']),
    'Regulation of cytokine response': ('#1a1a1a', ['Alox5', 'Mefv', 'Elane', 'Il17rb', 'Cd6']),
    'Type II IFN production': ('#a020a0', ['Irf8', 'Lgals9', 'Pglyrp1', 'Klrk1', 'Ccr2']),
}

gene_lookup = {g.upper(): g for g in df['Gene.name']}
legend_handles = []
for cat, (color, genes) in categories.items():
    legend_handles.append(ax.scatter([], [], color=color, s=60, edgecolor='black', linewidth=0.6, label=cat))
    for g in genes:
        if g.upper() in gene_lookup:
            row = df[df['Gene.name'] == gene_lookup[g.upper()]].iloc[0]
            ax.scatter(row['LT_Mendeley_logFC'], row['Fig4A_IAvsLT_logFC'], color=color, s=70,
                       edgecolor='black', linewidth=0.7, zorder=4)
            ax.annotate(g, (row['LT_Mendeley_logFC'], row['Fig4A_IAvsLT_logFC']),
                        xytext=(6, 4), textcoords='offset points', fontsize=9.5, fontweight='bold')

leg1 = ax.legend(handles=legend_handles, title='Fig.4B category', loc='upper left', fontsize=8.5, title_fontsize=9)
ax.add_artist(leg1)

other_handle = ax.scatter([], [], s=15, color='#d9d5c9', label='other')
sig_handle = ax.scatter([], [], s=40, color='#6a51a3', label='Fig.4A identity gene')
fit_line, = ax.plot([], [], '--', color='#d95f02', lw=1.8, label=f'OLS fit (slope={slope:.2f})')
ax.legend(handles=[other_handle, sig_handle, fit_line], loc='lower right', fontsize=8.5)

plt.tight_layout()
plt.savefig(f'{OUT}/assets/dotplot_LTmendeley_vs_Fig4A_logFC_labeled_redrawn.pdf', facecolor='white')
plt.savefig(f'{OUT}/assets/dotplot_LTmendeley_vs_Fig4A_logFC_labeled_redrawn.png', dpi=200, facecolor='white')
print('saved', len(df), rho, p, r, slope)
