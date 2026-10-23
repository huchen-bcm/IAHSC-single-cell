import scanpy as sc
import numpy as np
import pandas as pd
from scipy import stats
import sys
sys.path.insert(0, '/home/huc/projects/scAuto/scAnalytics')
from scAnalytics_helpers import score_gene_sets_ucell

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis'

a = sc.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/IAHSC.scRNA.hsc_subclasses.no_erythroid.Dec16.h5ad')

# WT only, LT-HSC / IA-HSC only (same subsetting as the original Group 1 analysis)
wt_mask = a.obs['sample'].astype(str).str.startswith(('17422', '17423'))
cls_mask = a.obs['HSC_subclass'].isin(['LT-HSC', 'IA-HSC'])
adata = a[wt_mask & cls_mask].copy()
adata.obs['Treatment'] = np.where(adata.obs['sample'].astype(str).str.startswith('17422'), 'Naive', 'Infected')
print(adata.obs.groupby(['HSC_subclass', 'Treatment'])['sample'].nunique())
print(adata)

# Build a raw AnnData view for UCell scoring (log-normalized, full transcriptome, per score_gene_sets_ucell's use_raw path)
def load_gmt_term(path, term):
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0] == term:
                return parts[2:]
    raise KeyError(term)

gmt = '/home/huc/projects/Rachana.metabolome/new_data/m5.go.bp.v2023.2.Mm.symbols.gmt'
raw_genes = set(adata.raw.var_names)

gene_sets = {
    'ALT_Type_II_IFN_response_GOBP_full': [g for g in load_gmt_term(gmt, 'GOBP_RESPONSE_TO_TYPE_II_INTERFERON') if g in raw_genes],
    'ALT_Leukocyte_chemotaxis_GOBP_full': [g for g in load_gmt_term(gmt, 'GOBP_LEUKOCYTE_CHEMOTAXIS') if g in raw_genes],
    'ALT_Myeloid_activation_GOBP_full': [g for g in load_gmt_term(gmt, 'GOBP_MYELOID_LEUKOCYTE_ACTIVATION') if g in raw_genes],
    'ALT_Mitotic_cell_cycle_GOBP_full': [g for g in load_gmt_term(gmt, 'GOBP_MITOTIC_CELL_CYCLE_PROCESS') if g in raw_genes],
}
for k, v in gene_sets.items():
    print(k, len(v), 'genes')

score_gene_sets_ucell(adata, gene_sets, use_raw=True)

score_cols = [k + '_gene_signature_ucell_score' for k in gene_sets.keys()]
per_sample = adata.obs.groupby(['HSC_subclass', 'Treatment', 'sample'], observed=True)[score_cols].mean().reset_index()
per_sample.to_csv(f'{OUT}/per_sample_scores_ALT_genesets.csv', index=False)
print(per_sample)

group_summary = per_sample.groupby(['HSC_subclass', 'Treatment'], observed=True)[score_cols].agg(['mean', 'std']).reset_index()
group_summary.to_csv(f'{OUT}/group_summary_scores_ALT_genesets.csv', index=False)

def get_vals(df, cls, treat, col):
    return df[(df.HSC_subclass == cls) & (df.Treatment == treat)][col].values

rows = []
for col in score_cols:
    lt_naive = get_vals(per_sample, 'LT-HSC', 'Naive', col)
    lt_inf = get_vals(per_sample, 'LT-HSC', 'Infected', col)
    ia_naive = get_vals(per_sample, 'IA-HSC', 'Naive', col)
    ia_inf = get_vals(per_sample, 'IA-HSC', 'Infected', col)
    rows.append({
        'score': col,
        'LT_naive_mean': lt_naive.mean(), 'LT_infected_mean': lt_inf.mean(),
        'IA_naive_mean': ia_naive.mean(), 'IA_infected_mean': ia_inf.mean(),
        'IAvsLT_naive_p': stats.ttest_ind(ia_naive, lt_naive).pvalue,
        'IAvsLT_infected_p': stats.ttest_ind(ia_inf, lt_inf).pvalue,
        'LT_infected_vs_naive_p': stats.ttest_ind(lt_inf, lt_naive).pvalue,
        'IA_infected_vs_naive_p': stats.ttest_ind(ia_inf, ia_naive).pvalue,
    })
res = pd.DataFrame(rows)
res.to_csv(f'{OUT}/sample_level_ttests_ALT_genesets.csv', index=False)
print(res.to_string())
