import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species/nmf_metaprogram'

jmat = pd.read_csv(f'{OUT}/jaccard_matrix.csv', index_col=0)
go_terms = json.load(open(f'{OUT}/mp_go_top_terms.json'))

MOUSE_LABEL = {
    'mouseMP1': 'S-phase/DNA replication', 'mouseMP2': 'Lymphoid-primed progenitor',
    'mouseMP3': 'LT-HSC quiescence/stemness', 'mouseMP4': 'Megakaryocyte (Cd41+)',
    'mouseMP5': 'Type I/II IFN-stimulated genes', 'mouseMP6': 'G2M/mitosis',
    'mouseMP7': 'Promyelocyte/granule biogenesis', 'mouseMP8': 'Mature neutrophil/inflammatory',
}
HUMAN_LABEL = {
    'humanMP1': 'Chromatin/RNA processing', 'humanMP2': 'Cytoplasmic translation (ribosomal)',
    'humanMP3': 'Cytoplasmic translation (ribosomal)', 'humanMP4': 'Immediate-early/stress response',
    'humanMP5': 'ER stress / MHC-I assembly (myeloid)', 'humanMP6': 'DNA replication/repair (cycling)',
}

TOP_USER_MOUSE = {'mouseMP1': 'Cd41+', 'mouseMP2': 'LT-HSC', 'mouseMP3': 'LT-HSC', 'mouseMP4': 'Cd41+',
                  'mouseMP5': 'LT-HSC', 'mouseMP6': 'IA-HSC', 'mouseMP7': 'ST-HSC', 'mouseMP8': 'IA-HSC'}
TOP_USER_HUMAN = {'humanMP1': 'c8', 'humanMP2': 'c0', 'humanMP3': 'c6', 'humanMP4': 'c7',
                  'humanMP5': 'c2', 'humanMP6': 'c5'}

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 6), gridspec_kw={'width_ratios': [1, 1.1]})

im = ax1.imshow(jmat.values, cmap='YlOrRd', vmin=0, vmax=jmat.values.max())
ax1.set_xticks(range(len(jmat.columns)))
ax1.set_xticklabels([f'{c}\n({HUMAN_LABEL[c]})\ntop:{TOP_USER_HUMAN[c]}' for c in jmat.columns], fontsize=7, rotation=90)
ax1.set_yticks(range(len(jmat.index)))
ax1.set_yticklabels([f'{r} ({MOUSE_LABEL[r]}) top:{TOP_USER_MOUSE[r]}' for r in jmat.index], fontsize=7.5)
for i in range(jmat.shape[0]):
    for j in range(jmat.shape[1]):
        v = jmat.values[i, j]
        if v > 0.03:
            ax1.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=7,
                      color='white' if v > jmat.values.max() * 0.6 else 'black')
plt.colorbar(im, ax=ax1, label='Jaccard similarity\n(top-50 genes)', fraction=0.046, pad=0.04)
ax1.set_title('Mouse NMF programs (rows) vs.\nhuman NMF programs (columns)', fontsize=11, fontweight='bold')

matches = pd.read_csv(f'{OUT}/mp_matches.csv')
matches = matches[matches['jaccard'] > 0.03].sort_values('jaccard', ascending=True)
ax2.barh(range(len(matches)), matches['jaccard'], color='#c0562f')
labels = []
for _, r in matches.iterrows():
    m, h = r['mouse_mp'], r['human_mp']
    labels.append(f"{m} ↔ {h}\n{MOUSE_LABEL[m]} / {HUMAN_LABEL[h]}")
ax2.set_yticks(range(len(matches)))
ax2.set_yticklabels(labels, fontsize=8)
ax2.set_xlabel('Jaccard similarity (best-match pairs)')
ax2.set_title('Cross-species meta-program matches', fontsize=11, fontweight='bold')
for i, v in enumerate(matches['jaccard']):
    ax2.text(v + 0.003, i, f'{v:.3f}', va='center', fontsize=8)

plt.tight_layout()
plt.savefig(f'{OUT}/../assets/nmf_metaprogram_matching.png', dpi=200, facecolor='white')
print('saved matching figure')
