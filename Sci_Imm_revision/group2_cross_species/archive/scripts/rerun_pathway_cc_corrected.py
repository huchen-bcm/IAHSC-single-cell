import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
import gseapy as gp

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

latent = ad.read_h5ad(f'{OUT}/embed_latent50.h5ad')
train = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/adata_merged_human_mouse_orthologs_hvg_integrated.0128_2026.h5ad')
assert (latent.obs_names == train.obs_names).all()

proj_df = pd.read_csv(f'{OUT}/axis_projection_cc_corrected.csv', index_col=0)
proj_new = proj_df.loc[latent.obs_names, 'proj_new'].values

mouse_mask = (latent.obs['system'] == 'mouse').values
human_mask = (latent.obs['system'] == 'human').values

X = train.X
if hasattr(X, 'toarray'):
    X = X.toarray()
genes = train.var_names.to_numpy()

mouse_proj = proj_new[mouse_mask]
human_proj = proj_new[human_mask]
Xm = X[mouse_mask]
Xh = X[human_mask]

print('computing mouse correlations (CC-corrected axis)...')
mouse_corr = np.zeros(len(genes))
for j in range(len(genes)):
    if Xm[:, j].std() > 0:
        mouse_corr[j] = spearmanr(Xm[:, j], mouse_proj).statistic
print('computing human correlations (CC-corrected axis)...')
human_corr = np.zeros(len(genes))
for j in range(len(genes)):
    if Xh[:, j].std() > 0:
        human_corr[j] = spearmanr(Xh[:, j], human_proj).statistic

res = pd.DataFrame({'gene': genes, 'mouse_corr': mouse_corr, 'human_corr': human_corr})
res.to_csv(f'{OUT}/integration_driver_genes_full_CC_CORRECTED.csv', index=False)

res['score'] = res['mouse_corr'] * res['human_corr']
rnk = res[['gene', 'score']].dropna().sort_values('score', ascending=False).drop_duplicates('gene').reset_index(drop=True)
print(f'{len(rnk)} genes ranked (CC-corrected axis)')

pre_res = gp.prerank(
    rnk=rnk,
    gene_sets='MSigDB_Hallmark_2020',
    organism='mouse',
    min_size=5,
    max_size=1000,
    permutation_num=1000,
    seed=0,
    threads=4,
    outdir=None,
)
gres = pre_res.res2d.copy()
gres['FDR q-val'] = gres['FDR q-val'].astype(float)
gres['NES'] = gres['NES'].astype(float)
gres = gres.sort_values('FDR q-val')
gres.to_csv(f'{OUT}/gsea_prerank_hallmark_CC_CORRECTED.csv', index=False)
print(gres[['Term', 'NES', 'NOM p-val', 'FDR q-val']].to_string())
print(f'\n{(gres["FDR q-val"] < 0.25).sum()} terms at FDR<0.25')
print(f'{(gres["FDR q-val"] < 0.05).sum()} terms at FDR<0.05')
