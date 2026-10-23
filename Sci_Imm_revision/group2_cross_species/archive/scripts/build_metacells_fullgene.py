import anndata as ad
import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

latent = ad.read_h5ad(f'{OUT}/embed_latent50.h5ad')

for species, cells_per_mc in [('mouse', 50), ('human', 50)]:
    full = ad.read_h5ad(f'{OUT}/{species}_fullgene_hvg6000.h5ad')
    lat_sub = latent[full.obs_names]
    assert (lat_sub.obs_names == full.obs_names).all()

    n_cells = full.n_obs
    K = max(20, n_cells // cells_per_mc)
    print(species, 'n_cells', n_cells, 'K metacells', K)

    Z = lat_sub.X
    km = MiniBatchKMeans(n_clusters=K, random_state=0, n_init=3, batch_size=1024).fit(Z)
    labels = km.labels_

    X = full.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    genes = full.var_names.to_numpy()

    if species == 'mouse':
        group_col = lat_sub.obs['HSC_subclass'].astype(str).values
    else:
        group_col = lat_sub.obs['hsc_clusters_lab'].astype(str).values

    mc_expr = np.zeros((K, X.shape[1]), dtype=np.float32)
    mc_size = np.zeros(K, dtype=int)
    mc_label = []
    for k in range(K):
        sel = labels == k
        mc_size[k] = sel.sum()
        if sel.sum() == 0:
            mc_label.append('NA')
            continue
        mc_expr[k] = X[sel].mean(axis=0)
        vals, counts = np.unique(group_col[sel], return_counts=True)
        mc_label.append(vals[np.argmax(counts)])

    keep = mc_size > 0
    mc_df = pd.DataFrame(mc_expr[keep], columns=genes, index=[f'{species}_mc{k}' for k in np.where(keep)[0]])
    mc_meta = pd.DataFrame({'size': mc_size[keep], 'dominant_label': np.array(mc_label)[keep]}, index=mc_df.index)
    mc_df.to_csv(f'{OUT}/metacells_{species}_fullgene_expr.csv')
    mc_meta.to_csv(f'{OUT}/metacells_{species}_fullgene_meta.csv')
    print(species, 'metacells saved:', mc_df.shape)
