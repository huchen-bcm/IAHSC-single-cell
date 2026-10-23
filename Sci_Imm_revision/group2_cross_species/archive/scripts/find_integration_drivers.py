import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

latent = ad.read_h5ad(f'{OUT}/embed_latent50.h5ad')
train = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/adata_merged_human_mouse_orthologs_hvg_integrated.0128_2026.h5ad')

assert (latent.obs_names == train.obs_names).all()

mouse_mask = (latent.obs['system'] == 'mouse').values
human_mask = (latent.obs['system'] == 'human').values

# discriminant axis: mouse IA-HSC centroid minus mouse LT-HSC centroid, in the 50-dim latent space
ia_mask = (latent.obs['HSC_subclass'] == 'IA-HSC').values
lt_mask = (latent.obs['HSC_subclass'] == 'LT-HSC').values
d = latent.X[ia_mask].mean(axis=0) - latent.X[lt_mask].mean(axis=0)
d = d / np.linalg.norm(d)

proj = latent.X @ d  # projection of every cell onto the discriminant axis
latent.obs['proj'] = proj

print('Validate axis on human side (mean projection by cluster):')
print(latent.obs.loc[human_mask].groupby('hsc_clusters_lab', observed=True)['proj'].mean().sort_values(ascending=False))

X = train.X
if hasattr(X, 'toarray'):
    X = X.toarray()
genes = train.var_names.to_numpy()

mouse_proj = proj[mouse_mask]
human_proj = proj[human_mask]
Xm = X[mouse_mask]
Xh = X[human_mask]

print('computing mouse correlations...')
mouse_corr = np.zeros(len(genes))
for j in range(len(genes)):
    if Xm[:, j].std() > 0:
        mouse_corr[j] = spearmanr(Xm[:, j], mouse_proj).statistic
print('computing human correlations...')
human_corr = np.zeros(len(genes))
for j in range(len(genes)):
    if Xh[:, j].std() > 0:
        human_corr[j] = spearmanr(Xh[:, j], human_proj).statistic

res = pd.DataFrame({'gene': genes, 'mouse_corr': mouse_corr, 'human_corr': human_corr})
res['concordance'] = np.sign(res['mouse_corr']) == np.sign(res['human_corr'])
res['min_abs_corr'] = res[['mouse_corr', 'human_corr']].abs().min(axis=1)
res = res.sort_values('min_abs_corr', ascending=False)
res.to_csv(f'{OUT}/integration_driver_genes_full.csv', index=False)

concordant_pos = res[(res['concordance']) & (res['mouse_corr'] > 0)].sort_values('min_abs_corr', ascending=False)
concordant_neg = res[(res['concordance']) & (res['mouse_corr'] < 0)].sort_values('min_abs_corr', ascending=False)
print(f'{len(concordant_pos)} genes concordantly POSITIVE (up toward IA-HSC/c3-4-7 in both species)')
print(concordant_pos.head(30))
print()
print(f'{len(concordant_neg)} genes concordantly NEGATIVE (up toward LT-HSC/other-clusters in both species)')
print(concordant_neg.head(30))

concordant_pos.to_csv(f'{OUT}/integration_driver_genes_concordant_positive.csv', index=False)
concordant_neg.to_csv(f'{OUT}/integration_driver_genes_concordant_negative.csv', index=False)
print('saved')
