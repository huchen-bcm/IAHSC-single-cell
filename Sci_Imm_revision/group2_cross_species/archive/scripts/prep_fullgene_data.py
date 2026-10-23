import anndata as ad
import numpy as np
import pandas as pd
import gzip
import re
import scanpy as sc

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

# ---- Ensembl -> symbol map for mouse (GRCm39, matches the scRNA reference) ----
gtf = '/home/huc/projects/ref_genome/Mus_musculus/refdata-gex-GRCm39-2024-A/genes/genes.gtf.gz'
ens2sym = {}
with gzip.open(gtf, 'rt') as f:
    for line in f:
        if line.startswith('#'):
            continue
        fields = line.split('\t')
        if fields[2] != 'gene':
            continue
        m_id = re.search(r'gene_id "([^"]+)"', fields[8])
        m_name = re.search(r'gene_name "([^"]+)"', fields[8])
        if m_id and m_name:
            ens2sym[m_id.group(1)] = m_name.group(1)
print('mouse ensembl->symbol map size:', len(ens2sym))

# ---- mouse full raw counts, SAME cell population as Part 0's WGCNA (all mouse cells in the
# actual trained integration model: both genotypes, non-Erythroid, 57,563 cells) ----
latent_check = ad.read_h5ad(f'{OUT}/embed_latent50.h5ad')
mouse_cells_in_model = latent_check.obs_names[latent_check.obs['system'] == 'mouse']

raw = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/merged_unintegrated.h5ad')
mouse = raw[mouse_cells_in_model].copy()
mouse.obs['HSC_subclass'] = latent_check.obs.loc[mouse_cells_in_model, 'HSC_subclass'].values
print('mouse cells (matches Part 0 exactly):', mouse.n_obs)

# map to symbol, drop genes without a symbol, dedupe (sum duplicate symbols)
mouse.var['symbol'] = [ens2sym.get(g, None) for g in mouse.var_names]
mouse = mouse[:, mouse.var['symbol'].notna()].copy()
mouse.var_names = mouse.var['symbol'].values
mouse.var_names_make_unique()
print('mouse genes with symbol:', mouse.n_vars)

sc.pp.normalize_total(mouse, target_sum=1e4)
sc.pp.log1p(mouse)
sc.pp.highly_variable_genes(mouse, n_top_genes=6000, flavor='seurat')
mouse_hvg = mouse[:, mouse.var['highly_variable']].copy()
print('mouse HVG matrix:', mouse_hvg.shape)
mouse_hvg.write(f'{OUT}/mouse_fullgene_hvg6000.h5ad')

# ---- human full raw counts, same cells as before (HSC_a/HSC_b) ----
human = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/h5ad/BCG_scRNA_HSC_annotated.h5ad')
human = human[human.obs['clust_names'].isin(['HSC_a', 'HSC_b'])].copy()
human.X = human.layers['raw_counts'].copy()
print('human cells:', human.n_obs, 'genes:', human.n_vars)

sc.pp.normalize_total(human, target_sum=1e4)
sc.pp.log1p(human)
sc.pp.highly_variable_genes(human, n_top_genes=6000, flavor='seurat')
human_hvg = human[:, human.var['highly_variable']].copy()
print('human HVG matrix:', human_hvg.shape)
human_hvg.write(f'{OUT}/human_fullgene_hvg6000.h5ad')

print('done')
