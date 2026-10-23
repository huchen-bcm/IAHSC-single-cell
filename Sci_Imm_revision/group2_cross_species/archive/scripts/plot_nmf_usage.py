import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species/nmf_metaprogram'

mouse_usage = pd.read_csv(f'{OUT}/mouse_usage.csv', index_col=0)
human_usage = pd.read_csv(f'{OUT}/human_usage.csv', index_col=0)
mouse_obs = pd.read_csv(f'{OUT}/mouse_obs.csv', index_col=0)
human_obs = pd.read_csv(f'{OUT}/human_obs.csv', index_col=0)

mouse_usage_norm = mouse_usage.div(mouse_usage.sum(axis=1), axis=0)
human_usage_norm = human_usage.div(human_usage.sum(axis=1), axis=0)

mu = mouse_usage_norm.join(mouse_obs['HSC_subclass'])
hu = human_usage_norm.join(human_obs['hsc_clusters_lab'])

mouse_order = ['Cd41+ HSC', 'LT-HSC', 'ST-HSC', 'IA-HSC']
human_order = [f'c{i}' for i in range(10)]

mouse_means = mu.groupby('HSC_subclass').mean().loc[mouse_order]
human_means = hu.groupby('hsc_clusters_lab').mean().loc[human_order]

MOUSE_LABEL = {
    'mouseMP1': 'MP1: S-phase/DNA repl.', 'mouseMP2': 'MP2: Lymphoid-primed prog.',
    'mouseMP3': 'MP3: LT-HSC quiescence', 'mouseMP4': 'MP4: Megakaryocyte',
    'mouseMP5': 'MP5: IFN-stimulated genes', 'mouseMP6': 'MP6: G2M/mitosis',
    'mouseMP7': 'MP7: Promyelocyte/granule', 'mouseMP8': 'MP8: Mature neutrophil/infl.',
}
HUMAN_LABEL = {
    'humanMP1': 'MP1: Chromatin/RNA proc.', 'humanMP2': 'MP2: Translation (ribo)',
    'humanMP3': 'MP3: Translation (ribo)', 'humanMP4': 'MP4: Immediate-early/stress',
    'humanMP5': 'MP5: ER stress/MHC-I (myeloid)', 'humanMP6': 'MP6: DNA repl./repair (cycling)',
}

fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

im0 = axes[0].imshow(mouse_means.values, cmap='viridis', aspect='auto')
axes[0].set_xticks(range(len(mouse_means.columns)))
axes[0].set_xticklabels([MOUSE_LABEL[c] for c in mouse_means.columns], rotation=45, ha='right', fontsize=8)
axes[0].set_yticks(range(len(mouse_order)))
axes[0].set_yticklabels(mouse_order, fontsize=9)
for i in range(mouse_means.shape[0]):
    for j in range(mouse_means.shape[1]):
        v = mouse_means.values[i, j]
        axes[0].text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=7.5,
                      color='white' if v < mouse_means.values.max() * 0.6 else 'black')
axes[0].set_title('Mouse: mean normalized program usage\nby HSC subclass (WT)', fontsize=11, fontweight='bold')
plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

im1 = axes[1].imshow(human_means.values, cmap='viridis', aspect='auto')
axes[1].set_xticks(range(len(human_means.columns)))
axes[1].set_xticklabels([HUMAN_LABEL[c] for c in human_means.columns], rotation=45, ha='right', fontsize=8)
axes[1].set_yticks(range(len(human_order)))
ylabels = [f'{c} *' if c in ('c3', 'c4', 'c7') else c for c in human_order]
axes[1].set_yticklabels(ylabels, fontsize=9)
for i in range(human_means.shape[0]):
    for j in range(human_means.shape[1]):
        v = human_means.values[i, j]
        axes[1].text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=7.5,
                      color='white' if v < human_means.values.max() * 0.6 else 'black')
axes[1].set_title('Human: mean normalized program usage\nby cluster (* = reviewer-flagged c3/c4/c7)', fontsize=11, fontweight='bold')
plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

plt.tight_layout()
plt.savefig(f'{OUT}/../assets/nmf_usage_by_group.png', dpi=200, facecolor='white')
print('saved usage figure')
