import pickle
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
r = pickle.load(open(f'{OUT}/figure3_enrichment_results.pkl', 'rb'))

groups = [
    ('mouseMP1__humanMP6', 'Reactome_Pathways_2024', 'mouseMP1↔humanMP6\nCell cycle/DNA replication', '#4a6fa5', 4),
    ('mouseMP3__humanMP2', 'Reactome_Pathways_2024', 'mouseMP3↔humanMP2\nLT-HSC quiescence-linked translation', '#c0562f', 4),
    ('mouseMP7__humanMP5', 'Reactome_Pathways_2024', 'mouseMP7↔humanMP5\nNeutrophil/ER stress/innate immune', '#5a9e6f', 4),
    ('mouseMP2__humanMP1', 'GO_Biological_Process_2025', 'mouseMP2↔humanMP1\nChromatin organization', '#8a6dab', 4),
]

rows = []  # (term, -log10p, overlap_n, color, group_label)
for key, gs, label, color, n in groups:
    df = r[f'{key}__{gs}'].head(n)
    for _, row in df.iterrows():
        term = row['Term']
        if len(term) > 42:
            term = term[:41] + '…'
        overlap_n = int(row['Overlap'].split('/')[0])
        rows.append((term, -np.log10(row['Adjusted P-value']), overlap_n, color, label))

fig, ax = plt.subplots(figsize=(8.5, 0.5 * len(rows) + 2.2))
y = np.arange(len(rows))[::-1]

for yi, (term, nlp, overlap_n, color, label) in zip(y, rows):
    ax.scatter(nlp, yi, s=overlap_n * 22, color=color, alpha=0.85, edgecolor='black', linewidth=0.5, zorder=3)

ax.set_yticks(y)
ax.set_yticklabels([row[0] for row in rows], fontsize=9)
ax.set_xlabel('-log10(adjusted p-value)', fontsize=10)
ax.axvline(-np.log10(0.05), color='grey', linestyle='--', lw=0.8, zorder=0)
ax.set_title('Figure 3 — Biological interpretation of conserved meta-programs\n(dot size = gene overlap count; color = matched mouse↔human MP pair)',
             fontsize=11.5, fontweight='bold')

# group divider lines + side labels
idx = 0
for key, gs, label, color, n in groups:
    n_actual = min(n, len(r[f'{key}__{gs}']))
    y0 = y[idx]
    y1 = y[idx + n_actual - 1]
    ax.axhspan(y1 - 0.5, y0 + 0.5, color=color, alpha=0.07, zorder=0)
    ax.text(-0.3, (y0 + y1) / 2, label, ha='right', va='center', fontsize=8, color=color, fontweight='bold',
            transform=ax.get_yaxis_transform())
    idx += n_actual

ax.set_xlim(0, max(r2[1] for r2 in rows) * 1.15)
plt.tight_layout()
plt.savefig(f'{OUT}/assets/figure3_conserved_mp_dotplot.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved figure 3')
