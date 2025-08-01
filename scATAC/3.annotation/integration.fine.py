import snapatac2 as snap
import glob, os
from pathlib import Path
from typing import List, Tuple
import anndata
import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
import scvi
import sys
from plotly.subplots import make_subplots


sys.path.append(os.path.abspath('/home/huc/projects/scAuto/scAnnotation'))
# sys.path.append(os.path.abspath('/home/huc/projects/scAuto/scQC'))
sys.path.append(os.path.abspath('/home/huc/projects/scRNA_resources/markers'))

import scAnnotation_helpers as scAnnotation_helpers
import marker_helpers as marker_helpers
# import helpers  as scQC_helpers

import importlib
importlib.reload(scAnnotation_helpers)
importlib.reload(marker_helpers)

marker_yaml_file = "/home/huc/projects/scRNA_resources/markers/marker_collections/mouse_bone_marrow.yaml"


scRNA = sc.read_h5ad("/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/IAHSC.scRNA.annotation_finalized.June16.h5ad")
scRNA_count_data = sc.read_h5ad("/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/merged_unintegrated.h5ad")

cell_indices = scRNA.obs.index.tolist()
count_adata_ = scRNA_count_data[cell_indices, ]
count_adata_.obs['Cell_types_fine_finalized']= scRNA.obs['Cell_types_fine_finalized'].tolist()

gtf_path = "/home/huc/projects/ref_genome/Mus_musculus/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz"
organism = "mouse"
gencode_data, mito_g, ribo_g, cyc_g = scAnnotation_helpers.process_gene_annotation(gtf_path, organism = organism)

count_adata_.var.index = gencode_data.drop_duplicates('gene').set_index('gene').reindex(index=count_adata_.var.index)['name'].tolist()
count_adata_ = count_adata_[:, ~count_adata_.var.index.duplicated()]
# count_adata_.var.index = gencode_data.drop_duplicates('gene').set_index('gene').reindex(index=count_adata_.var.index)['name'].tolist()

count_adata_.var.index.duplicated().sum()

gene_matrix = snap.read("IAHSC_gene_mat.WT.July18.h5ad", backed = None)

gene_matrix.obs['Cell_types_fine_finalized'] = pd.NA
joint = ad.concat(
    [count_adata_, gene_matrix],
    join='inner',
    label='modality',
    keys=["scRNA", "scATAC"],
    index_unique='_',
)

sc.pp.filter_genes(joint, min_cells=5)
sc.pp.highly_variable_genes(
    joint,
    n_top_genes = 5000,
    flavor="seurat_v3",
    batch_key="modality",
    subset=True
)

scvi.model.SCVI.setup_anndata(joint, batch_key="modality")
# scvi.model.SCVI.setup_anndata(joint, batch_key="modality", layer = "counts")
vae = scvi.model.SCVI(
    joint,
    n_layers=2,
    n_latent=30,
    gene_likelihood="nb",
    dispersion="gene-batch",
)


vae.train(max_epochs=1000, early_stopping=True)


# ax = vae.history['elbo_train'][1:].plot()
# vae.history['elbo_validation'].plot(ax=ax)

joint.obs["celltype_scanvi"] = 'Unknown'
ref_idx = joint.obs['modality'] == "scRNA"
joint.obs["celltype_scanvi"][ref_idx] = joint.obs['Cell_types_fine_finalized'][ref_idx]

lvae = scvi.model.SCANVI.from_scvi_model(
    vae,
    adata=joint,
    labels_key="celltype_scanvi",
    unlabeled_category="Unknown",
)

lvae.train(max_epochs=1000, n_samples_per_label=100)


joint.obs["C_scANVI"] = lvae.predict(joint)
joint.obsm["X_scANVI"] = lvae.get_latent_representation(joint)

sc.pp.neighbors(joint, use_rep="X_scANVI")
sc.tl.umap(joint)

## save the joint anndata object
joint.write("IAHSC_gene_mat.WT.July18.Cell_types_fine_finalized.h5ad", compression='gzip')



# atac.obs['cell_type'] = joint.obs.loc[atac.obs_names + '_query']['C_scANVI'].to_numpy()



# sc.pl.umap(atac, color=['cell_type', "leiden"], wspace=0.45)


# atac.write("pbmc5k_annotated.h5ad", compression="gzip")
