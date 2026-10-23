import anndata as ad
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import os
import gc
import sys
import warnings
warnings.filterwarnings('ignore')

VELOCITY_DIR = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/raw_repro_velocity'
META_CSV = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/FINAL_phates_meta.csv'
meta = pd.read_csv(META_CSV, index_col=0)


def loom_to_meta_bc(cellid):
    cap, bc = cellid.split(':')
    bc = bc.rstrip('x')
    return f'{cap}_{bc}-1'


def run_once(sample_uid, captures, out_csv):
    sample_meta = meta[meta['UNIQUE_ID'] == sample_uid].copy()
    adatas = []
    for cap in captures:
        a = ad.read_loom(f'{VELOCITY_DIR}/{cap}/{cap}.loom')
        a.var_names_make_unique()
        adatas.append(a)
    adata = ad.concat(adatas, join='inner') if len(adatas) > 1 else adatas[0]
    adata.obs['meta_bc'] = [loom_to_meta_bc(c) for c in adata.obs_names]
    keep = adata.obs['meta_bc'].isin(sample_meta.index)
    adata = adata[keep].copy()
    adata.obs_names = adata.obs['meta_bc'].values
    adata.obs_names_make_unique()
    sample_meta_aligned = sample_meta.loc[adata.obs_names]
    adata.obs['Cluster_colors'] = sample_meta_aligned['clust_names'].values
    adata.obsm['X_umap'] = sample_meta_aligned[['PHATE1', 'PHATE2']].values

    scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
    scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
    scv.tl.recover_dynamics(adata, n_jobs=8)
    scv.tl.velocity(adata, mode='dynamical')
    scv.tl.velocity_graph(adata)

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
    trans_probs.to_csv(out_csv, index=False)

    most_likely = trans_probs.set_index('Cell ID').idxmax(axis=1)
    hsc_mask = meta.loc[meta.index.intersection(trans_probs['Cell ID']), 'clust_names'].isin(['HSC_a', 'HSC_b'])
    hsc_cells = hsc_mask[hsc_mask].index
    ml_hsc = most_likely.loc[most_likely.index.intersection(hsc_cells)]
    pct = (ml_hsc == 'CMP_b').mean() * 100
    print(f'{out_csv}: n_hsc={len(ml_hsc)}, pct_CMP_b={pct:.3f}%')
    del adata
    gc.collect()


if __name__ == '__main__':
    which = sys.argv[1]
    run_once('SNG-LB-SS-1S-RS-S20-CD34neg_S19_R1_001_Tm3_BCG', ['capture4'],
              f'/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/rerun_test/S20_Tm3_BCG_{which}.csv')
