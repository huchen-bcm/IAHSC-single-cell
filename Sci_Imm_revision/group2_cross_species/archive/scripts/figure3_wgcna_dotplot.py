import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
r = pickle.load(open(f'{OUT}/figure3_wgcna_enrichment.pkl', 'rb'))

groups = [
    ('mouse_firebrick', 'mouse firebrick\n(IA-HSC top module)', '#c0562f', 3),
    ('human_darkgrey', 'human darkgrey\n(r=0.906, c7>c4>c3)', '#8a3a0a', 3),
    ('mouse_red', 'mouse red\n(not axis-concordant)', '#4a6fa5', 3),
    ('human_coral', 'human coral\n(not axis-concordant)', '#9aa0a6', 3),
]

rows = []
for key, label, color, n in groups:
    df = r[key].head(n)
    for _, row in df.iterrows():
        term = row['Term']
        if len(term) > 42:
            term = term[:41] + '…'
        overlap_n = int(row['Overlap'].split('/')[0])
        rows.append((term, -np.log10(row['Adjusted P-value']), overlap_n, color, label))

fig, ax = plt.subplots(figsize=(9, 0.5 * len(rows) + 2.2))
y = np.arange(len(rows))[::-1]

for yi, (term, nlp, overlap_n, color, label) in zip(y, rows):
    ax.scatter(nlp, yi, s=overlap_n * 12, color=color, alpha=0.85, edgecolor='black', linewidth=0.5, zorder=3)

ax.set_yticks(y)
ax.set_yticklabels([row[0] for row in rows], fontsize=9)
ax.set_xlabel('-log10(adjusted p-value)', fontsize=10)
ax.axvline(-np.log10(0.05), color='grey', linestyle='--', lw=0.8, zorder=0)
ax.set_title('Figure 3 (WGCNA version) — Reactome enrichment: the one axis-concordant\npair (firebrick↔darkgrey) vs. a significant-but-discordant pair (red↔coral)',
             fontsize=11.5, fontweight='bold')

idx = 0
for key, label, color, n in groups:
    n_actual = min(n, len(r[key]))
    y0 = y[idx]
    y1 = y[idx + n_actual - 1]
    ax.axhspan(y1 - 0.5, y0 + 0.5, color=color, alpha=0.08, zorder=0)
    ax.text(-0.3, (y0 + y1) / 2, label, ha='right', va='center', fontsize=8, color=color, fontweight='bold',
            transform=ax.get_yaxis_transform())
    idx += n_actual

ax.set_xlim(0, max(r2[1] for r2 in rows) * 1.15)
plt.tight_layout()
plt.savefig(f'{OUT}/assets/figure3_wgcna_dotplot.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved')
