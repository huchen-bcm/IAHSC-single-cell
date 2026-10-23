
# This script describes the overall cellRank workflow and gives an example of how each sample was processed 

#################
##  overview   ##
#################

#First run velocyto to generate loom files for each capture
#Then run the Rscript "01_make_input_cellrank_wVelocyto" in this cellrank folder to generate inputs for cellrank
#run module load python -> conda activate cellrank -> python



#################
##  overview   ##
#################

#The following is an example of what was run for each capture (capture 2 used as an example)
#The actual python scripts run for each capture are located in the "python scripts" folder 


#CELL RANK estimators and kernels (low-level mode)

import scvelo as scv
import pandas as pd
import numpy as np
import matplotlib as plt
import anndata

import os
os.chdir("")

adata = anndata.read_loom("/capture2_output/velocyto/capture2_output.loom")
sample_obs = pd.read_csv("cellID_obs.csv")
umap_cord = pd.read_csv("cell_embeddings.csv")
cell_clusters = pd.read_csv("clusters.csv")

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

scv.pp.filter_and_normalize(adata, min_shared_counts=20, n_top_genes=2000)
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=30)
scv.pp.moments(adata, n_pcs=None, n_neighbors=None)
scv.tl.recover_dynamics(adata, n_jobs=8)
scv.tl.velocity(adata, mode="dynamical")
scv.tl.velocity_graph(adata)

scv.pl.velocity_embedding_stream(adata, basis="umap", legend_fontsize=12, title="", smooth=0.8, min_mass=4, color="Cluster_colors", save="test_velocity_plot.png")


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
#g.plot_spectrum()
g.compute_macrostates(n_states=11, cluster_key="Cluster_colors")  #use 11 but ok if it defaults to 12
g.plot_macrostates(same_plot=False, save="test_macrostates_separate.png")
g.plot_macrostates(discrete=True, save="test_macrostates_discrete.png")

g.set_terminal_states_from_macrostates(names=["MEP_b", "CMP_c", "CMP_a", "CMP_b_1", "CMP_b_2", "GMP_b", "PreBNK"]) #combine macrostates by original clusters
g.plot_terminal_states()


g.compute_absorption_probabilities(use_petsc=True, n_jobs=5, solver='gmres')
#g.plot_absorption_probabilities(adata, same_plot=False)  #THIS PART NEEDS MORE THAN 50GB MEM - GIVES AN ERROR


##TO SAVE RESULTS

#example of how data frames can be saved
adata.obsm["to_terminal_states"]
trans_probs = pd.DataFrame(adata.obsm["to_terminal_states"])
trans_probs = trans_probs.rename(columns = {0:'MEP_b', 1:'CMPc'})  #rename all columns like this manually according to cluster
trans_probs['Cell ID'] = adata_index
trans_probs.to_csv('test_to_terminal_states.csv')

adata.obsm["terminal_states_memberships"]
memb = pd.DataFrame(adata.obsm["terminal_states_memberships"])
memb = memb.rename(columns = {0:'MEP_b', 1:'CMPc'})  #rename all columns like this manually according to cluster
memb['Cell ID'] = adata_index
memb.to_csv('test_terminal_states_memberships.csv')
