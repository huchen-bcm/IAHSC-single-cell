import pandas as pd
import numpy as np
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

mouse_profile = pd.read_csv(f'{OUT}/wgcna_mouse_ownmodule_by_subtype_zscore.csv', index_col=0)
human_profile = pd.read_csv(f'{OUT}/wgcna_human_ownmodule_by_cluster_zscore.csv', index_col=0)
go = json.load(open(f'{OUT}/wgcna_all_modules_Reactome_top_terms.json'))
mouse_mod = pd.read_csv(f'{OUT}/wgcna_mouse_gene_modules.csv', index_col=0)
human_mod = pd.read_csv(f'{OUT}/wgcna_human_gene_modules.csv', index_col=0)
mouse_sizes = mouse_mod['moduleColors'].value_counts()
human_sizes = human_mod['moduleColors'].value_counts()


def short_term(terms, maxlen=38):
    if not terms:
        return '(no sig. Reactome term)'
    t = terms[0].split(' (GO:')[0]
    return t if len(t) <= maxlen else t[:maxlen - 1] + '…'


mouse_labels = []
for c in [x.replace('mouse_', '') for x in mouse_profile.columns]:
    mouse_labels.append(f'{c} (n={mouse_sizes[c]})\n{short_term(go.get("mouse_" + c, []))}')

human_labels = []
for c in [x.replace('human_', '') for x in human_profile.columns]:
    human_labels.append(f'{c} (n={human_sizes[c]})\n{short_term(go.get("human_" + c, []))}')

fig, axes = plt.subplots(1, 2, figsize=(16, 6.5), gridspec_kw={'width_ratios': [1, 1.15]})

vmax = 2.5
im0 = axes[0].imshow(mouse_profile.values.T, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
axes[0].set_xticks(range(len(mouse_profile.index)))
axes[0].set_xticklabels(mouse_profile.index, rotation=30, ha='right')
axes[0].set_yticks(range(len(mouse_labels)))
axes[0].set_yticklabels(mouse_labels, fontsize=8)
axes[0].set_title('Mouse: module score by HSC subclass\n(own modules; labels = top Reactome term)', fontsize=11, fontweight='bold')
for i in range(mouse_profile.shape[1]):
    for j in range(mouse_profile.shape[0]):
        axes[0].text(j, i, f'{mouse_profile.values[j,i]:.1f}', ha='center', va='center', fontsize=7,
                      color='white' if abs(mouse_profile.values[j, i]) > vmax * 0.55 else 'black')
plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.06, label='z-score')

im1 = axes[1].imshow(human_profile.values.T, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
axes[1].set_xticks(range(len(human_profile.index)))
axes[1].set_xticklabels(human_profile.index, rotation=30, ha='right')
axes[1].set_yticks(range(len(human_labels)))
axes[1].set_yticklabels(human_labels, fontsize=8)
axes[1].set_title('Human: module score by cluster\n(own modules; labels = top Reactome term)', fontsize=11, fontweight='bold')
for i in range(human_profile.shape[1]):
    for j in range(human_profile.shape[0]):
        axes[1].text(j, i, f'{human_profile.values[j,i]:.1f}', ha='center', va='center', fontsize=7,
                      color='white' if abs(human_profile.values[j, i]) > vmax * 0.55 else 'black')
for xt, c in zip(axes[1].get_xticklabels(), human_profile.index):
    if c in ('c3', 'c4', 'c7'):
        xt.set_color('#c0562f')
        xt.set_fontweight('bold')
plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.06, label='z-score')

plt.tight_layout()
plt.savefig(f'{OUT}/assets/wgcna_module_by_subtype_heatmap_reactome.png', dpi=200, facecolor='white')
print('saved')
