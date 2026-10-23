import sys
import anndata as ad
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')

SAMPLE_UID = sys.argv[1]  # e.g. "SNG-LB-SS-1S-RS-S20-CD34neg_S19_R1_001_Td0_BCG"
CAPTURE = sys.argv[2]  # e.g. "capture3"
LABEL = sys.argv[3]  # e.g. "S20_Td0_BCG"

VELOCITY_DIR = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/raw_repro_velocity'
OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results'
META_CSV = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/FINAL_phates_meta.csv'

meta = pd.read_csv(META_CSV, index_col=0)
sample_meta = meta[(meta['UNIQUE_ID'] == SAMPLE_UID) & (meta['batchID'] == CAPTURE)].copy()
print(f'{LABEL}: {len(sample_meta)} cells in metadata for this sample/capture')

import os
os.makedirs(OUT, exist_ok=True)
cache_path = f'{OUT}/{LABEL}_velocity.h5ad'

if os.path.exists(cache_path):
    print(f'{LABEL}: loading cached velocity result from {cache_path}')
    adata = ad.read_h5ad(cache_path)
else:
    adata = ad.read_loom(f'{VELOCITY_DIR}/{CAPTURE}/{CAPTURE}.loom')
    adata.var_names_make_unique()

    # convert loom CellID "capture3:AAACCCAGTCCGACGTx" -> metadata rowname "capture3_AAACCCAGTCCGACGT-1"
    def loom_to_meta_bc(cellid):
        cap, bc = cellid.split(':')
        bc = bc.rstrip('x')
        return f'{cap}_{bc}-1'

    adata.obs['meta_bc'] = [loom_to_meta_bc(c) for c in adata.obs_names]
    keep = adata.obs['meta_bc'].isin(sample_meta.index)
    print(f'{LABEL}: {keep.sum()} / {adata.n_obs} loom cells matched to sample')
    adata = adata[keep].copy()
    adata.obs_names = adata.obs['meta_bc'].values

    sample_meta_aligned = sample_meta.loc[adata.obs_names]
    adata.obs['Cluster_colors'] = sample_meta_aligned['clust_names'].values
    adata.obsm['X_umap'] = sample_meta_aligned[['PHATE1', 'PHATE2']].values

    print(f'{LABEL}: cluster composition:')
    print(adata.obs['Cluster_colors'].value_counts())

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
    scv.tl.recover_dynamics(adata, n_jobs=8)
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)
    adata.write(cache_path)

from cellrank.kernels import VelocityKernel, ConnectivityKernel
from cellrank.estimators import GPCCA

vk = VelocityKernel(adata)
vk.compute_transition_matrix()
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck

g = GPCCA(combined_kernel)
n_states_present = adata.obs['Cluster_colors'].nunique()
n_schur = min(20, n_states_present + 5, adata.n_obs - 2)
g.compute_schur(n_components=n_schur)

present_clusters = [c for c in adata.obs['Cluster_colors'].unique() if pd.notna(c) and c not in ('HSC_a', 'HSC_b')]
terminal_dict = {c: adata[adata.obs['Cluster_colors'] == c].obs_names for c in present_clusters}
g.set_terminal_states(terminal_dict)

g.compute_fate_probabilities(solver='gmres')
trans_probs = pd.DataFrame(adata.obsm['lineages_fwd'], columns=adata.obsm['lineages_fwd'].names.tolist())
trans_probs['Cell ID'] = adata.obs_names
trans_probs.to_csv(f'{OUT}/{LABEL}_to_terminal_states.csv', index=False)

print(f'{LABEL}: DONE. Terminal states used: {present_clusters}')
print(trans_probs.head())
