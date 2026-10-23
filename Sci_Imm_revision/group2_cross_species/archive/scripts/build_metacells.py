import anndata as ad
import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

latent = ad.read_h5ad(f'{OUT}/embed_latent50.h5ad')
train = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/adata_merged_human_mouse_orthologs_hvg_integrated.0128_2026.h5ad')
assert (latent.obs_names == train.obs_names).all()

ia_mask = (latent.obs['HSC_subclass'] == 'IA-HSC').values
lt_mask = (latent.obs['HSC_subclass'] == 'LT-HSC').values
d = latent.X[ia_mask].mean(axis=0) - latent.X[lt_mask].mean(axis=0)
d = d / np.linalg.norm(d)
proj = latent.X @ d

X = train.X
if hasattr(X, 'toarray'):
    X = X.toarray()
genes = train.var_names.to_numpy()

for species, cells_per_mc in [('mouse', 50), ('human', 50)]:
    mask = (latent.obs['system'] == species).values
    n_cells = mask.sum()
    K = max(20, n_cells // cells_per_mc)
    print(species, 'n_cells', n_cells, 'K metacells', K)

    Z = latent.X[mask]
    km = MiniBatchKMeans(n_clusters=K, random_state=0, n_init=3, batch_size=1024).fit(Z)
    labels = km.labels_

    Xs = X[mask]
    projs = proj[mask]
    hsc_sub = latent.obs.loc[mask, 'HSC_subclass'].values if species == 'mouse' else None
    hsc_clu = latent.obs.loc[mask, 'hsc_clusters_lab'].values if species == 'human' else None

    mc_expr = np.zeros((K, Xs.shape[1]), dtype=np.float32)
    mc_proj = np.zeros(K)
    mc_size = np.zeros(K, dtype=int)
    mc_label = []
    for k in range(K):
        sel = labels == k
        mc_size[k] = sel.sum()
        if sel.sum() == 0:
            mc_label.append('NA')
            continue
        mc_expr[k] = Xs[sel].mean(axis=0)
        mc_proj[k] = projs[sel].mean()
        if species == 'mouse':
            vals, counts = np.unique(hsc_sub[sel], return_counts=True)
            mc_label.append(vals[np.argmax(counts)])
        else:
            vals, counts = np.unique(hsc_clu[sel].astype(str), return_counts=True)
            mc_label.append(vals[np.argmax(counts)])

    keep = mc_size > 0
    mc_df = pd.DataFrame(mc_expr[keep], columns=genes, index=[f'{species}_mc{k}' for k in np.where(keep)[0]])
    mc_meta = pd.DataFrame({'proj': mc_proj[keep], 'size': mc_size[keep], 'dominant_label': np.array(mc_label)[keep]},
                            index=mc_df.index)
    mc_df.to_csv(f'{OUT}/metacells_{species}_expr.csv')
    mc_meta.to_csv(f'{OUT}/metacells_{species}_meta.csv')
    print(species, 'metacells saved:', mc_df.shape)
