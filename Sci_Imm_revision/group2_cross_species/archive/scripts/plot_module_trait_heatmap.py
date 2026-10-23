import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'


def build_heatmap(species, group_col_order, ax_fig_path, title, module_sizes):
    me = pd.read_csv(f'{OUT}/wgcna_{species}_fullgene_module_eigengenes.csv', index_col=0)
    fg_meta = pd.read_csv(f'{OUT}/metacells_{species}_fullgene_meta.csv', index_col=0)
    proj_meta = pd.read_csv(f'{OUT}/metacells_{species}_meta.csv', index_col=0)  # same metacell groups, has 'proj'
    assert (fg_meta.index == me.index).all()
    assert (proj_meta.index == me.index).all()

    labels = fg_meta['dominant_label'].values
    proj = proj_meta['proj'].values

    traits = pd.DataFrame(index=me.index)
    for g in group_col_order:
        traits[g] = (labels == g).astype(int)
    traits['axis_proj'] = proj

    modules = me.columns.tolist()
    trait_names = traits.columns.tolist()

    r_mat = pd.DataFrame(index=modules, columns=trait_names, dtype=float)
    p_mat = pd.DataFrame(index=modules, columns=trait_names, dtype=float)
    for m in modules:
        for t in trait_names:
            r, p = pearsonr(me[m].astype(float), traits[t].astype(float))
            r_mat.loc[m, t] = r
            p_mat.loc[m, t] = p

    # order modules by correlation with the axis (descending), matching the main narrative
    order = r_mat['axis_proj'].sort_values(ascending=False).index
    r_mat = r_mat.loc[order]
    p_mat = p_mat.loc[order]

    fig, ax = plt.subplots(figsize=(1.35 * len(trait_names) + 2.2, 0.34 * len(modules) + 1.6))
    im = ax.imshow(r_mat.values, cmap='RdBu_r', vmin=-1, vmax=1, aspect='auto')

    for i in range(r_mat.shape[0]):
        for j in range(r_mat.shape[1]):
            r = r_mat.values[i, j]
            p = p_mat.values[i, j]
            pstr = f'{p:.0e}' if p < 0.001 else f'{p:.3f}'
            ax.text(j, i, f'{r:.2f}\n({pstr})', ha='center', va='center', fontsize=6.2,
                     color='white' if abs(r) > 0.6 else 'black')

    ax.set_xticks(range(len(trait_names)))
    xt = [t if t != 'axis_proj' else 'axis\nprojection' for t in trait_names]
    ax.set_xticklabels(xt, rotation=40, ha='right', fontsize=9)
    ax.set_yticks(range(len(modules)))
    ylabels = [f'{m.replace("ME","")} (n={module_sizes.get(m.replace("ME",""),"?")})' for m in r_mat.index]
    ax.set_yticklabels(ylabels, fontsize=8)
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label('Pearson r\n(module eigengene vs. trait)', fontsize=8.5)

    for xt_lbl, name in zip(ax.get_xticklabels(), trait_names):
        if name in ('c3', 'c4', 'c7', 'IA-HSC'):
            xt_lbl.set_color('#c0562f')
            xt_lbl.set_fontweight('bold')

    plt.tight_layout()
    plt.savefig(ax_fig_path, dpi=200, facecolor='white')
    plt.close(fig)
    r_mat.to_csv(ax_fig_path.replace('.png', '_r.csv'))
    p_mat.to_csv(ax_fig_path.replace('.png', '_p.csv'))
    print('saved', ax_fig_path)


mouse_mod = pd.read_csv(f'{OUT}/wgcna_mouse_fullgene_gene_modules.csv', index_col=0)
human_mod = pd.read_csv(f'{OUT}/wgcna_human_fullgene_gene_modules.csv', index_col=0)
mouse_sizes = mouse_mod['moduleColors'].value_counts().to_dict()
human_sizes = human_mod['moduleColors'].value_counts().to_dict()

build_heatmap('mouse', ['IA-HSC', 'LT-HSC', 'ST-HSC', 'Cd41+ HSC'],
              f'{OUT}/assets/wgcna_fullgene_module_trait_mouse.png',
              'Mouse module-trait relationships (full transcriptome WGCNA)', mouse_sizes)

build_heatmap('human', [f'c{i}' for i in range(10)],
              f'{OUT}/assets/wgcna_fullgene_module_trait_human.png',
              'Human module-trait relationships (full transcriptome WGCNA)', human_sizes)
