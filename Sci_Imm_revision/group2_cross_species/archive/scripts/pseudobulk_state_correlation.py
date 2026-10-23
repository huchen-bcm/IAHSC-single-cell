import anndata as ad
import numpy as np
import pandas as pd

NMF_DIR = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species/nmf_metaprogram'
OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

with open(f'{NMF_DIR}/shared_ortholog_genes.txt') as f:
    shared_genes = [l.strip() for l in f if l.strip()]
print('shared ortholog genes:', len(shared_genes))

mouse = ad.read_h5ad(f'{NMF_DIR}/mouse_ortholog_raw.h5ad')
human = ad.read_h5ad(f'{NMF_DIR}/human_ortholog_raw.h5ad')
human.X = human.layers['raw_counts']

mouse = mouse[:, shared_genes].copy()
human = human[:, shared_genes].copy()


def pseudobulk(adata, group_col, order):
    X = adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    groups = adata.obs[group_col].astype(str).values
    pb = np.zeros((len(order), X.shape[1]))
    for i, g in enumerate(order):
        pb[i] = X[groups == g].sum(axis=0)
    cpm = pb / pb.sum(axis=1, keepdims=True) * 1e6
    logcpm = np.log2(cpm + 1)
    return pd.DataFrame(logcpm, index=order, columns=adata.var_names)


mouse_order = ['IA-HSC', 'LT-HSC', 'ST-HSC', 'Cd41+ HSC']
human_order = [f'c{i}' for i in range(10)]

mouse_pb = pseudobulk(mouse, 'HSC_subclass', mouse_order)
human_pb = pseudobulk(human, 'hsc_clusters_lab', human_order)

mouse_pb.to_csv(f'{OUT}/pseudobulk_mouse_state_logcpm.csv')
human_pb.to_csv(f'{OUT}/pseudobulk_human_cluster_logcpm.csv')

corr = pd.DataFrame(index=mouse_order, columns=human_order, dtype=float)
for m in mouse_order:
    for h in human_order:
        corr.loc[m, h] = np.corrcoef(mouse_pb.loc[m], human_pb.loc[h])[0, 1]

print(corr.round(3))
corr.to_csv(f'{OUT}/pseudobulk_crossspecies_state_correlation.csv')

print()
print('IA-HSC ranking:')
print(corr.loc['IA-HSC'].sort_values(ascending=False))
print()
print('LT-HSC ranking:')
print(corr.loc['LT-HSC'].sort_values(ascending=False))
