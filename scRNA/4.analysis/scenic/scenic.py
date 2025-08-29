# Import required packages and set up the environment
import scanpy as sc
import omicverse as ov
import numpy as np
import pandas as pd
import sys
import os
import matplotlib.pyplot as plt

# Set up plotting parameters
#ov.plot_set(font_path='Arial')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'



def process_gene_annotation(gtf_path, organism):
    # Validate organism input
    if organism not in ['human', 'mouse']:
        raise ValueError(f"Organism {organism} is not supported. Please choose either 'human' or 'mouse'.")

    # Process GTF info
    if gtf_path.endswith('.zip'):
        # Open the zip file and read the first file inside it
        with zipfile.ZipFile(gtf_path) as z:
            # Get the list of files in the zip archive
            file_list = z.namelist()
            # Open the first file in the archive
            with z.open(file_list[0]) as f:
                gencode_data = pd.read_table(f, skiprows=5, header=None, dtype = str)
    else:
        # If it's not a zip file, read it normally
        gencode_data = pd.read_table(gtf_path, skiprows=5, header=None, dtype = str)

    # gencode_data = pd.read_table(gtf_path, skiprows=5, header=None)

    gencode_data = gencode_data[gencode_data[2] == 'transcript']

    gencode_data['tx'] = [i.split('transcript_id "')[-1].split('"')[0] for i in gencode_data[8]]
    gencode_data['gene'] = [i.split('gene_id "')[-1].split('"')[0].split('.')[0] for i in gencode_data[8]]
    gencode_data['name'] = [i.split('gene_name "')[-1].split('"')[0] for i in gencode_data[8]]
    gencode_data['gene_type'] = [i.split('gene_type "')[-1].split('"')[0] for i in gencode_data[8]]
    gencode_data['tx_type'] = [i.split('transcript_type "')[-1].split('"')[0] for i in gencode_data[8]]
    gencode_data['status'] = [i.split('transcript_status "')[-1].split('"')[0] for i in gencode_data[8]]
    gencode_data['support'] = [i.split('transcript_support_level "')[-1].split('"')[0] for i in gencode_data[8]]
    gencode_data.drop_duplicates('tx', inplace=True)

    # Return processed gene lists
    return gencode_data



## Load gene annotations
gencode_data = process_gene_annotation('/home/huc/projects/ref_genome/Mus_musculus/refdata-gex-mm10-2020-A/genes/genes.gtf', organism='mouse')

## Load data
count_adata = sc.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/merged_unintegrated.h5ad')
adata_all = sc.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/IAHSC.scRNA.annotation_finalized.with_meta_info.Aug01.h5ad')
adata_HSC_only = sc.read_h5ad("/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/IAHSC.scRNA.HSC_only.Aug15.h5ad")

N_CELLS = 2000

if adata_HSC_only.n_obs > N_CELLS:
    rng = np.random.default_rng(42)  # reproducible
    idx = rng.choice(adata_HSC_only.n_obs, N_CELLS, replace=False)
    adata_HSC_only = adata_HSC_only[idx, :].copy()


## Get the counts
adata_HSC_only_counts = count_adata[adata_HSC_only.obs.index, ]
sc.pp.filter_genes(adata_HSC_only_counts, min_cells=5)
adata_HSC_only_counts.var.index = gencode_data.drop_duplicates('gene').set_index('gene').reindex(index=adata_HSC_only_counts.var.index)['name'].tolist()
adata_HSC_only_counts = adata_HSC_only_counts[:, ~adata_HSC_only_counts.var.index.duplicated()]
adata_HSC_only_counts

adata_HSC_only_counts.layers['raw_counts'] = adata_HSC_only_counts.X.copy()


# ---- 0) Pre-filter genes & cast to float32 (biggest win) ----
adata = adata_HSC_only_counts.copy()

# Keep reasonably expressed & variable genes
sc.pp.filter_genes(adata, min_cells=int(adata.n_obs * 0.01))  # ≥1% cells
sc.pp.highly_variable_genes(adata, flavor="seurat_v3", n_top_genes=6000)
adata = adata[:, adata.var['highly_variable']].copy()

# Optional: drop MT/ribo/pseudo unless needed biologically
# adata = adata[:, ~adata.var_names.str.startswith(('MT-', 'RPL', 'RPS'))]

# Make operations cheaper
if not np.issubdtype(adata.X.dtype, np.floating):
    adata.X = adata.X.astype(np.float32)
else:
    adata.X = adata.X.astype(np.float32)


# Set paths to the required database files
# Update these paths to match your local installation
db_glob = "/home/huc/projects/scRNA_resources/scenic_auxiliaries/mm10/*feather"
motif_path = "/home/huc/projects/scRNA_resources/scenic_auxiliaries/mm10/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

# The db_glob pattern should match both ranking database files
# The motif_path should point to the TF-motif annotation file



# Set paths to the required database files
# Update these paths to match your local installation
db_glob = "/home/huc/projects/scRNA_resources/scenic_auxiliaries/mm10/*feather"
motif_path = "/home/huc/projects/scRNA_resources/scenic_auxiliaries/mm10/motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl"

# The db_glob pattern should match both ranking database files
# The motif_path should point to the TF-motif annotation file



# Initialize the SCENIC object
scenic_obj = ov.single.SCENIC(
    adata=adata_HSC_only_counts,                    # Single-cell expression data
    db_glob=db_glob,                # Pattern for ranking database files
    motif_path=motif_path,          # TF-motif annotation file
    n_jobs=25                       # Number of parallel processes
)

print("SCENIC object initialized successfully!")
print(f"Number of ranking databases loaded: {len(scenic_obj.dbs)}")
print(f"Database names: {[db.name for db in scenic_obj.dbs]}")
print(f"Number of genes in dataset: {adata_HSC_only_counts.n_vars}")
print(f"Number of cells in dataset: {adata_HSC_only_counts.n_obs}")




# Run GRN inference using RegDiffusion
print("Starting GRN inference with RegDiffusion...")
print("This typically takes 5-10 minutes for datasets of this size.")

# Perform GRN inference
edgelist = scenic_obj.cal_grn(layer='raw_counts', device = "cuda")

print(f"\nGRN inference completed!")
print(f"Number of potential regulatory relationships: {len(edgelist):,}")
print(f"Number of unique TFs: {edgelist['TF'].nunique()}")
print(f"Number of unique target genes: {edgelist['target'].nunique()}")
print(f"Importance score range: {edgelist['importance'].min():.3f} - {edgelist['importance'].max():.3f}")

# Show first few regulatory relationships
print("\nTop 10 regulatory relationships:")
edgelist.head(10)



# # 0) Keep a handle to the full object
# adata_full = scenic_obj.adata

# # 1) Choose indices for the downsampled cells
# N_CELLS = 20000
# if adata_full.n_obs > N_CELLS:
#     ridx = np.random.RandomState(42).choice(adata_full.n_obs, N_CELLS, replace=False)
# else:
#     ridx = np.arange(adata_full.n_obs)

# # 2) (Optional) ensure float32 for speed/RAM
# adata_sub = adata_full[ridx, :].copy()
# if adata_sub.X.dtype != np.float32:
#     adata_sub.X = adata_sub.X.astype(np.float32)

# 3) Trim fan-out per TF BEFORE module step (huge speed-up)
K = 50
edgelist_trim = (edgelist
                 .sort_values(["TF", "importance"], ascending=[True, False])
                 .groupby("TF", group_keys=False).head(K))

# 4) Run regulon building on the smaller adata
# try:
# scenic_obj.adata = adata_sub  # <<< makes modules_from_adjacencies use the subset
scenic_obj.edgelist = edgelist_trim  # if the object reads from its own edgelist
# If omicverse exposes a way to skip AUCell here, set it (e.g., calc_auc=False). If not, we’ll overwrite later.
regulon_ad = scenic_obj.cal_regulons(
    rho_mask_dropouts=False,   # faster during module building
    thresholds=(0.75, 0.9),
    top_n_targets=(30,),       # faster first pass; increase later if needed
    top_n_regulators=(5, 10, 30),
    )
# finally:
#     # 5) Restore full data on the object
#     scenic_obj.adata = adata_full


print("\nAUCell scoring completed")

#################### Save results
print("\nSaving results...")
# Save the SCENIC object (contains all analysis results)
ov.utils.save(scenic_obj, 'scenic_obj.pkl')

# Save the regulon activity AnnData object
regulon_ad.write('scenic_regulon_ad.h5ad')
