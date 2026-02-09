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

# Load custom helper functions
sys.path.append(os.path.abspath('/home/huc/projects/scAuto/scAnnotation'))
sys.path.append(os.path.abspath('/home/huc/projects/scRNA_resources/markers'))

import scAnnotation_helpers as scAnnotation_helpers
import marker_helpers as marker_helpers
# import helpers  as scQC_helpers

import importlib
importlib.reload(scAnnotation_helpers)
importlib.reload(marker_helpers)

# Load marker list
marker_yaml_file = "/home/huc/projects/scRNA_resources/markers/marker_collections/mouse_bone_marrow.yaml"

# Load gene annotation
gtf_path = "/home/huc/projects/ref_genome/Mus_musculus/refdata-cellranger-arc-mm10-2020-A-2.0.0/genes/genes.gtf.gz"
organism = "mouse"
gencode_data, mito_g, ribo_g, cyc_g = scAnnotation_helpers.process_gene_annotation(gtf_path, organism = organism)


# Load scRNA-seq data
scRNA = sc.read_h5ad("/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/IAHSC.scRNA.annotation_finalized.June16.h5ad")
scRNA_count_data = sc.read_h5ad("/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/merged_unintegrated.h5ad")
cell_indices = scRNA.obs.index.tolist()
count_adata_ = scRNA_count_data[cell_indices, ]
count_adata_.obs['Cell_types_fine_finalized']= scRNA.obs['Cell_types_fine_finalized'].tolist()
count_adata_.var.index = gencode_data.drop_duplicates('gene').set_index('gene').reindex(index=count_adata_.var.index)['name'].tolist()
count_adata_ = count_adata_[:, ~count_adata_.var.index.duplicated()]
count_adata_.var.index.duplicated().sum()
meta_table = pd.read_csv("/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/sample_meta.txt", sep="\t")
meta_table = meta_table.astype(str)
count_adata_.obs = count_adata_.obs.merge(meta_table.set_index('SampleID'), left_on='sample', right_index=True, how='left')
count_adata_ = count_adata_[count_adata_.obs['Genotype'] == 'WT', :]

# only consider HSC
count_adata_ = count_adata_[count_adata_.obs['Cell_types_fine_finalized'].isin(['ST-HSC', "LT-HSC", "Cd41+ HSC", "IA-HSC"]), :]


# Load scATAC-seq gene activity matrix
scATAC_HSC_WT = snap.read_dataset("/home/huc/projects/Ruoqiong_Cao.IAHSC/scATAC/3.annotation/snap_HSC_subset/HSC_dataset_new2.h5ads")
scATAC_HSC_WT_GeneAct = snap.pp.make_gene_matrix(scATAC_HSC_WT, gene_anno=snap.genome.mm10)
scATAC_HSC_WT.close()


gene_matrix = snap.read("/home/huc/projects/Ruoqiong_Cao.IAHSC/scATAC/3.annotation/gene_activity_WT.annotation.finalized.Nov03_copy.h5ad", backed = None)

# Integration
scATAC_HSC_WT_GeneAct.obs['Cell_types_fine_finalized'] = pd.NA

joint = ad.concat(
    [count_adata_, scATAC_HSC_WT_GeneAct],
    join='inner',
    label='modality',
    keys=["scRNA", "scATAC"],
    index_unique='_',
)


sc.pp.filter_genes(joint, min_cells=5)
sc.pp.highly_variable_genes(
    joint,
    n_top_genes = 3000,
    flavor="seurat_v3",
    batch_key="modality",
    subset=True
)

scvi.model.SCVI.setup_anndata(joint, batch_key="modality")
vae = scvi.model.SCVI(
    joint,
    n_layers=2,
    n_latent=50,
    gene_likelihood="nb",
    dispersion="gene-batch",
)

vae.train(max_epochs=1000, early_stopping=True)



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

# sc.pp.neighbors(joint, use_rep="X_scANVI")
# sc.tl.umap(joint)

## save the joint anndata object
joint.write("/home/huc/projects/Ruoqiong_Cao.IAHSC/scATAC/3.annotation/HSC.scVI.scRNA_scATAC_integrated.0204_2026_v2.h5ad", compression='gzip')


## save models
import os
save_path = "./HSC_scVI_scRNA_scATAC_scanvi_model_0204_2026_v2"

# Create the directory if it doesn't exist
if not os.path.exists(save_path):
    os.makedirs(save_path)

lvae.save(save_path, overwrite=True)
