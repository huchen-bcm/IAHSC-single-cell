import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

XIE = '/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/Xie'

adata = sc.read_h5ad(f'{XIE}/adata_merged_Xie_mouse_orthologs_hvg_integrated.best_model.h5ad')
obs = pd.read_csv('/tmp/xie_merged_obs.csv', index_col=0)
assert (adata.obs_names == obs.index).all()

proj = pd.read_csv('/tmp/xie_mouse_axis_cc.csv', index_col=0)  # mouse cells only, has 'proj'
# also need proj for human cells -- recompute quickly from saved embed + axis is simplest:
embed_all = np.load('/tmp/xie_embed_all100.npy')
embed30 = embed_all[:, :30]
mouse_mask = (obs['system'] == 'mouse').values
human_mask = (obs['system'] == 'human').values
ia_mask = (obs['HSC_subclass'] == 'IA-HSC').values
lt_mask = (obs['HSC_subclass'] == 'LT-HSC').values
d = embed30[ia_mask].mean(axis=0) - embed30[lt_mask].mean(axis=0)
d = d / np.linalg.norm(d)
proj_all = embed30 @ d

X = adata.X
if hasattr(X, 'toarray'):
    X = X.toarray()
genes = adata.var_names.to_numpy()

Xm = X[mouse_mask]
Xh = X[human_mask]
mouse_proj = proj_all[mouse_mask]
human_proj = proj_all[human_mask]

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
res['concordant'] = np.sign(res['mouse_corr']) == np.sign(res['human_corr'])
res['min_abs_corr'] = res[['mouse_corr', 'human_corr']].abs().min(axis=1)
res.to_csv(f'{XIE}/xie_integration_driver_genes_full.csv', index=False)
print(f'{len(res)} genes total')

for thresh in [0, 0.1, 0.2]:
    sub = res[(res['mouse_corr'] > thresh) & (res['human_corr'] > thresh)]
    print(f'r>{thresh}: {len(sub)} genes')
