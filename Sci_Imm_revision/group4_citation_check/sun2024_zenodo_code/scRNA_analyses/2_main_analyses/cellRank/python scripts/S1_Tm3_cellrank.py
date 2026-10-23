
import anndata
import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt

import os
os.chdir('/scRNA_analyses/2_main_analyses/cellRank/wVelocity_input')

adata = anndata.read_loom("/capture10_output/velocyto/capture10_output.loom")
sample_obs = pd.read_csv("cellID_obs_S1_Tm3.csv")
umap_cord = pd.read_csv("cell_embeddings_S1_Tm3.csv")
cell_clusters = pd.read_csv("clusters_S1_Tm3.csv")

adata = adata[np.isin(adata.obs.index,sample_obs["x"])]

adata_index = pd.DataFrame(adata.obs.index)
adata_index = adata_index.rename(columns = {'CellID':'Cell ID'})
umap_cord = umap_cord.rename(columns = {'Unnamed: 0':'Cell ID'})
umap_ordered = adata_index.merge(umap_cord, on = "Cell ID")
umap_ordered = umap_ordered.iloc[:,1:]
adata.obsm['X_umap'] = umap_ordered.values

cell_clusters = cell_clusters.rename(columns = {'Unnamed: 0':'Cell ID'})
clusters_ordered = adata_index.merge(cell_clusters, on = "Cell ID")
clusters_ordered = clusters_ordered.iloc[:,1:]

adata.uns['Cluster_colors'] = clusters_ordered.values
adata.obs['Cluster_colors'] = clusters_ordered.values

import sys
import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
os.chdir('/scRNA_analyses/2_main_analyses/cellRank/velocity_plots')

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis="umap", legend_fontsize=12, title="S1_Tm3 Velocity Plot", smooth=0.8, min_mass=4, color="Cluster_colors", save="S1_Tm3_velocity_plot.png")

adata.obsm['velocity_umap']
vu = pd.DataFrame(adata.obsm['velocity_umap'])
vu['Cell ID'] = adata_index
vu = vu.rename(columns = {0:'coord1', 1:'coord2'})
vu.to_csv('S1_Tm3_velocity_umap_coords.csv')


os.chdir('/scRNA_analyses/2_main_analyses/cellRank/terminal_state_probs')


from cellrank.tl.kernels import VelocityKernel
vk = VelocityKernel(adata)
print(vk)

vk.compute_transition_matrix()
from cellrank.tl.kernels import ConnectivityKernel
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
from cellrank.tl.estimators import GPCCA

g = GPCCA(combined_kernel)
print(g)

g.compute_schur(n_components=20)

g.set_terminal_states({"MLP_a": adata[adata.obs["Cluster_colors"] == "MLP_a"].obs_names, "GMP_b": adata[adata.obs["Cluster_colors"]=="GMP_b"].obs_names, 
"MLP_b": adata[adata.obs["Cluster_colors"] == "MLP_b"].obs_names, "GMP_a": adata[adata.obs["Cluster_colors"] == "GMP_a"].obs_names, "PreBNK": adata[adata.obs["Cluster_colors"] == "PreBNK"].obs_names,
"MEP_a": adata[adata.obs["Cluster_colors"] == "MEP_a"].obs_names, "MEP_b": adata[adata.obs["Cluster_colors"] == "MEP_b"].obs_names,
"MEP_c": adata[adata.obs["Cluster_colors"] == "MEP_c"].obs_names, "CMP_a": adata[adata.obs["Cluster_colors"] == "CMP_a"].obs_names, "CMP_b": adata[adata.obs["Cluster_colors"] == "CMP_b"].obs_names,
"CMP_c": adata[adata.obs["Cluster_colors"] == "CMP_c"].obs_names})
g.plot_terminal_states(discrete="TRUE", title="S1_Tm3 Terminal states", save="S1_Tm3_terminal_states.png")

g.compute_absorption_probabilities(use_petsc=True, n_jobs=5, solver='gmres')
adata.obsm["to_terminal_states"]
trans_probs = pd.DataFrame(adata.obsm["to_terminal_states"])
trans_probs = trans_probs.rename(columns = {0:'CMP_a', 1:'CMP_b', 2:'CMP_c', 3:'GMP_a', 4:'GMP_b', 5:'MEP_a', 6:'MEP_b', 7:'MEP_c', 8:'MLP_a', 9:'MLP_b', 10:'PreBNK'})  #rename all columns like this manually according to cluster
trans_probs['Cell ID'] = adata_index
trans_probs.to_csv('S1_Tm3_to_terminal_states.csv')
