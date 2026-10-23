import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
NMF_DIR = f'{OUT}/nmf_metaprogram'

jac = pd.read_csv(f'{NMF_DIR}/jaccard_matrix.csv', index_col=0)
loading = pd.read_csv(f'{NMF_DIR}/mp_loading_correlation.csv', index_col=0)
loading = loading.loc[jac.index, jac.columns]

MOUSE_LABEL = {
    'mouseMP1': 'S-phase/DNA repl.', 'mouseMP2': 'Lymphoid-primed prog.',
    'mouseMP3': 'LT-HSC quiescence', 'mouseMP4': 'Megakaryocyte',
    'mouseMP5': 'IFN-stimulated genes', 'mouseMP6': 'G2M/mitosis',
    'mouseMP7': 'Promyelocyte/granule', 'mouseMP8': 'Mature neutrophil/infl.',
}
HUMAN_LABEL = {
    'humanMP1': 'Chromatin/RNA proc.', 'humanMP2': 'Translation (ribo)',
    'humanMP3': 'Translation (ribo)', 'humanMP4': 'Immediate-early/stress',
    'humanMP5': 'ER stress/MHC-I (myeloid)', 'humanMP6': 'DNA repl./repair (cycling)',
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, mat, title, cmap, vlim, fmt in [
    (axes[0], jac, 'Top-50-gene Jaccard similarity', 'YlOrRd', (0, jac.values.max()), '{:.2f}'),
    (axes[1], loading, 'Full gene-loading correlation\n(Pearson, 15,739 shared genes)', 'RdBu_r', (-0.4, 0.4), '{:.2f}'),
]:
    im = ax.imshow(mat.values, cmap=cmap, vmin=vlim[0], vmax=vlim[1], aspect='auto')
    ax.set_xticks(range(len(mat.columns)))
    ax.set_xticklabels([f'{c}\n{HUMAN_LABEL[c]}' for c in mat.columns], fontsize=7.5, rotation=45, ha='right')
    ax.set_yticks(range(len(mat.index)))
    ax.set_yticklabels([f'{r} ({MOUSE_LABEL[r]})' for r in mat.index], fontsize=8)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = mat.values[i, j]
            thresh = vlim[1] * 0.6 if cmap == 'YlOrRd' else 0.25
            ax.text(j, i, fmt.format(v), ha='center', va='center', fontsize=7,
                     color='white' if abs(v) > thresh else 'black')
    ax.set_title(title, fontsize=11, fontweight='bold')
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.06)

fig.suptitle('Figure 2 — Conserved transcriptional programs: mouse NMF meta-programs × human NMF meta-programs',
             fontsize=12.5, fontweight='bold', y=1.03)
plt.tight_layout()
plt.savefig(f'{OUT}/assets/figure2_mp_correspondence.png', dpi=200, facecolor='white', bbox_inches='tight')
print('saved figure 2')

# also save a combined "best matches by either metric" table
matches = []
for m in jac.index:
    for h in jac.columns:
        matches.append((m, h, jac.loc[m, h], loading.loc[m, h]))
mdf = pd.DataFrame(matches, columns=['mouse_mp', 'human_mp', 'jaccard', 'loading_corr'])
mdf['combined_rank'] = mdf['jaccard'].rank(ascending=False) + mdf['loading_corr'].rank(ascending=False)
mdf = mdf.sort_values('combined_rank')
mdf.to_csv(f'{OUT}/figure2_mp_matches_combined.csv', index=False)
print(mdf.head(8).to_string())
