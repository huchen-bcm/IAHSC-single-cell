import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis'

scores = pd.read_csv(f'{OUT}/scRNA_pseudobulk_gsva_scores_full.csv', index_col=0)
innate = pd.read_csv(f'{OUT}/scRNA_pseudobulk_gsva_scores_InnateImmuneResponse_GOBP.csv', index_col=0).T
leading_edge = pd.read_csv(f'{OUT}/scRNA_pseudobulk_gsva_scores_GSEA_leading_edge.csv', index_col=0).T
scores = pd.concat([scores, innate, leading_edge])

meta = pd.read_csv(f'{OUT}/scRNA_pseudobulk_meta.csv', index_col=0)
meta['cond'] = meta['HSC_subclass'] + ' ' + meta['Treatment']
group_means = scores.T.join(meta['cond']).groupby('cond').mean(numeric_only=True).T

cols_map = {'LT-HSC\nnaive': 'LT-HSC Naive', 'LT-HSC\ninfected': 'LT-HSC Infected',
            'IA-HSC\nnaive': 'IA-HSC Naive', 'IA-HSC\ninfected': 'IA-HSC Infected'}
group_means = group_means[list(cols_map.values())]
group_means.columns = list(cols_map.keys())

# GOBP full term rows -- GO IDs confirmed against MSigDB "Exact source" (gene counts cross-checked
# against the m5.go.bp.v2023.2.Mm.symbols.gmt entries used to build these gene sets, except
# InnateImmuneResponse_GOBP which was newly scored -- see scripts/run_gsva_innate.R history)
label_map = {
    'TypeII_IFN_GOBP': 'Type II IFN response (GO:0034341)',
    'LeukocyteChemotaxis_GOBP': 'Leukocyte chemotaxis (GO:0030595)',
    'MyeloidActivation_GOBP': 'Myeloid leukocyte activation (GO:0002274)',
    'InnateImmuneResponse_GOBP': 'Innate immune response (GO:0045087)',
    'CellCycle_GOBP_MitoticProcess': 'Mitotic cell cycle process (GO:1903047)',
    'CellCycle_GOBP_G2M_transition': 'Cell cycle G2/M phase transition (GO:0044839)',
    # leading-edge genes extracted directly from bulk_IAHSC_vs_LTHSC_gsea/IAHSC_vs_LTHSC.bulk.gsea.GO_Biological_Process_2025.csv
    # (gseapy prerank "Lead_genes" column), best-matching positive-NES term per category -- NOT read off
    # any manuscript figure. See scripts/run_gsva_gsea_leading_edge.R.
    'TypeII_IFN_LeadingEdge': 'Cellular response to type II interferon (GO:0071346)',
    'LeukocyteChemotaxis_LeadingEdge': 'Positive regulation of leukocyte chemotaxis (GO:0002690)',
    'MacrophageActivation_LeadingEdge': 'Regulation of macrophage activation (GO:0043030)',
    'ResponseToBacterium_LeadingEdge': 'Defense response to bacterium (GO:0042742)',
    'CytokineProduction_LeadingEdge': 'Regulation of cytokine production (GO:0001817)',
}

cc_rows = ['CellCycle_GOBP_MitoticProcess', 'CellCycle_GOBP_G2M_transition']
pathway_rows = ['TypeII_IFN_GOBP', 'LeukocyteChemotaxis_GOBP', 'MyeloidActivation_GOBP', 'InnateImmuneResponse_GOBP']
curated_rows = ['TypeII_IFN_LeadingEdge', 'LeukocyteChemotaxis_LeadingEdge', 'MacrophageActivation_LeadingEdge',
                'ResponseToBacterium_LeadingEdge', 'CytokineProduction_LeadingEdge']


def plot_heatmap(rows, out_name, title):
    df = group_means.loc[rows].rename(index=label_map)
    vmax = np.abs(df.values).max()
    fig, ax = plt.subplots(figsize=(5.5, max(1.6, 0.4 * len(df) + 1)))
    ax.imshow(df.values, cmap='RdBu_r', vmin=-vmax, vmax=vmax, aspect='auto')
    ax.set_xticks(range(len(df.columns)))
    ax.set_xticklabels(df.columns, fontsize=8)
    ax.set_yticks(range(len(df.index)))
    ax.set_yticklabels(df.index, fontsize=7)
    for i in range(df.shape[0]):
        for j in range(df.shape[1]):
            v = df.values[i, j]
            text_color = 'white' if abs(v) / vmax > 0.5 else 'black'
            ax.text(j, i, f'{v:.2f}', ha='center', va='center', fontsize=7, color=text_color)
    ax.set_title(title, fontsize=9)
    for spine in ax.spines.values():
        spine.set_visible(True)
    plt.tight_layout()
    plt.savefig(f'{OUT}/assets/{out_name}.pdf', facecolor='white')
    plt.savefig(f'{OUT}/assets/{out_name}.png', dpi=200, facecolor='white')
    print('saved', out_name)


plot_heatmap(pathway_rows + cc_rows, 'four_condition_heatmap_pseudobulk_gsva_GOBPfull',
             'Four-condition pathway/cell-cycle GSVA scores (GOBP full gene sets)\n(scRNA pseudobulk, same algorithm as bulk)')

plot_heatmap(pathway_rows, 'four_condition_heatmap_pseudobulk_gsva_GOBPfull_noCC',
             'Four-condition pathway GSVA scores (GOBP full gene sets, cell cycle removed)\n(scRNA pseudobulk, same algorithm as bulk)')

plot_heatmap(curated_rows, 'four_condition_heatmap_pseudobulk_gsva_leading_edge',
             'Four-condition pathway GSVA scores (GSEA leading-edge gene sets)\n(scRNA pseudobulk, same algorithm as bulk)')
