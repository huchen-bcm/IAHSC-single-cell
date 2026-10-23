import anndata
import scvelo as scv
import scanpy as sc
import pandas as pd
import numpy as np
import os
import gc
import warnings
warnings.filterwarnings('ignore')

VELOCITY_DIR = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/raw_repro_velocity'
OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group4_citation_check/cellrank_reanalysis/repro_results_v1'
META_CSV = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/FINAL_phates_meta.csv'
os.makedirs(OUT, exist_ok=True)

meta = pd.read_csv(META_CSV, index_col=0)


def loom_to_meta_bc(cellid):
    cap, bc = cellid.split(':')
    bc = bc.rstrip('x')
    return f'{cap}_{bc}-1'


def run_one(sample_uid, captures, label):
    tp_out = f'{OUT}/{label}_to_terminal_states.csv'
    if os.path.exists(tp_out):
        print(f'{label}: already done, skipping')
        return

    sample_meta = meta[meta['UNIQUE_ID'] == sample_uid].copy()

    adatas = []
    for cap in captures:
        a = anndata.read_loom(f'{VELOCITY_DIR}/{cap}/{cap}.loom')
        a.var_names_make_unique()
        adatas.append(a)
    adata = adatas[0].concatenate(adatas[1:], index_unique=None) if len(adatas) > 1 else adatas[0]

    adata.obs['meta_bc'] = [loom_to_meta_bc(c) for c in adata.obs_names]
    keep = adata.obs['meta_bc'].isin(sample_meta.index)
    print(f'{label}: {keep.sum()} / {adata.n_obs} loom cells matched to sample (captures: {captures})')
    adata = adata[keep].copy()
    adata.obs_names = adata.obs['meta_bc'].values
    adata.obs_names_make_unique()

    sample_meta_aligned = sample_meta.loc[adata.obs_names]
    adata.obs['Cluster_colors'] = sample_meta_aligned['clust_names'].values
    adata.obsm['X_umap'] = sample_meta_aligned[['PHATE1', 'PHATE2']].values

    vc = adata.obs['Cluster_colors'].value_counts()
    print(f'{label}: cluster composition: {vc.to_dict()}')

    try:
        # exactly matching Sun et al.'s own script
        scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
        scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
        scv.tl.recover_dynamics(adata, n_jobs=1)
        scv.tl.velocity(adata, mode='dynamical')
        scv.tl.velocity_graph(adata)

        from cellrank.tl.kernels import VelocityKernel, ConnectivityKernel
        from cellrank.tl.estimators import GPCCA

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

        g.compute_absorption_probabilities(use_petsc=True, n_jobs=1, solver='gmres')
        trans_probs = pd.DataFrame(adata.obsm['to_terminal_states'], columns=present_clusters)
        trans_probs['Cell ID'] = adata.obs_names
        trans_probs.to_csv(tp_out, index=False)
        print(f'{label}: DONE. Terminal states used: {present_clusters}')
    except Exception as e:
        print(f'{label}: FAILED: {e}')
    finally:
        del adata
        gc.collect()


if __name__ == '__main__':
    samples = [
        ("SNG-LB-SS-1S-RS-S1-CD34neg_S1_R1_001_Td0_BCG", ['capture9'], 'S1_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S1-CD34neg_S1_R1_001_Tm3_BCG", ['capture10', 'capture14'], 'S1_Tm3_BCG'),
        ("SNG-LB-SS-1S-RS-S5-CD34neg_S4_R1_001_Td0_BCG", ['capture7'], 'S5_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S5-CD34neg_S4_R1_001_Tm3_BCG", ['capture8'], 'S5_Tm3_BCG'),
        ("SNG-LB-SS-1S-RS-S6-CD34neg_S5_R1_001_Td0_BCG", ['capture7'], 'S6_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S6-CD34neg_S5_R1_001_Tm3_BCG", ['capture8'], 'S6_Tm3_BCG'),
        ("SNG-LB-SS-1S-RS-S7-CD34neg_S6_R1_001_Td0_CTL", ['capture7'], 'S7_Td0_CTL'),
        ("SNG-LB-SS-1S-RS-S7-CD34neg_S6_R1_001_Tm3_CTL", ['capture8'], 'S7_Tm3_CTL'),
        ("SNG-LB-SS-1S-RS-S8-CD34neg_S7_R1_001_Td0_BCG", ['capture1'], 'S8_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S8-CD34neg_S7_R1_001_Tm3_BCG", ['capture2', 'capture14'], 'S8_Tm3_BCG'),
        ("SNG-LB-SS-1S-RS-S10-CD34neg_S9_R1_001_Td0_BCG", ['capture9'], 'S10_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S10-CD34neg_S9_R1_001_Tm3_BCG", ['capture6', 'capture10'], 'S10_Tm3_BCG'),
        ("SNG-LB-SS-1S-RS-S11-CD34neg_S10_R1_001_Td0_CTL", ['capture9'], 'S11_Td0_CTL'),
        ("SNG-LB-SS-1S-RS-S11-CD34neg_S10_R1_001_Tm3_CTL", ['capture10'], 'S11_Tm3_CTL'),
        ("SNG-LB-SS-1S-RS-S13-CD34neg_S12_R1_001_Td0_CTL", ['capture5'], 'S13_Td0_CTL'),
        ("SNG-LB-SS-1S-RS-S13-CD34neg_S12_R1_001_Tm3_CTL", ['capture6'], 'S13_Tm3_CTL'),
        ("SNG-LB-SS-1S-RS-S14-CD34neg_S13_R1_001_Td0_BCG", ['capture13'], 'S14_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S14-CD34neg_S13_R1_001_Tm3_BCG", ['capture14'], 'S14_Tm3_BCG'),
        ("SNG-LB-SS-1S-RS-S15-CD34neg_S14_R1_001_Td0_BCG", ['capture3'], 'S15_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S15-CD34neg_S14_R1_001_Tm3_BCG", ['capture4', 'capture6'], 'S15_Tm3_BCG'),
        ("SNG-LB-SS-1S-RS-S16-CD34neg_S15_R1_001_Td0_BCG", ['capture5', 'capture7'], 'S16_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S16-CD34neg_S15_R1_001_Tm3_BCG", ['capture6'], 'S16_Tm3_BCG'),
        ("SNG-LB-SS-1S-RS-S17-CD34neg_S16_R1_001_Td0_BCG", ['capture13'], 'S17_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S17-CD34neg_S16_R1_001_Tm3_BCG", ['capture14'], 'S17_Tm3_BCG'),
        ("SNG-LB-SS-1S-RS-S18-CD34neg_S17_R1_001_Td0_CTL", ['capture13'], 'S18_Td0_CTL'),
        ("SNG-LB-SS-1S-RS-S18-CD34neg_S17_R1_001_Tm3_CTL", ['capture14'], 'S18_Tm3_CTL'),
        ("SNG-LB-SS-1S-RS-S20-CD34neg_S19_R1_001_Td0_BCG", ['capture3'], 'S20_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S20-CD34neg_S19_R1_001_Tm3_BCG", ['capture4'], 'S20_Tm3_BCG'),
        ("SNG-LB-SS-1S-RS-S21-CD34neg_S20_R1_001_Td0_BCG", ['capture5'], 'S21_Td0_BCG'),
        ("SNG-LB-SS-1S-RS-S21-CD34neg_S20_R1_001_Tm3_BCG", ['capture6', 'capture8'], 'S21_Tm3_BCG'),
    ]
    for uid, caps, label in samples:
        print(f'=== starting {label} ===')
        run_one(uid, caps, label)
