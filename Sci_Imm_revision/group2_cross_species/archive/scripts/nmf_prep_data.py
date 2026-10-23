import anndata as ad
import pandas as pd
import numpy as np
import scipy.sparse as sp

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species/nmf_metaprogram'

# ---- 1:1 orthologs ---- (mouse raw-counts file is indexed by Ensembl ID, not symbol)
orth = pd.read_csv('/home/huc/projects/Geri.MDS/CUTRUN/DiffBind/mouse_human_marmoset_macaque_orthologs_20231113.csv')
orth = orth[['human_Symbol', 'mouse_Symbol', 'mouse_EnsemblID']].dropna()
orth = orth[(orth['human_Symbol'] != '') & (orth['mouse_Symbol'] != '') & (orth['mouse_EnsemblID'] != '')]
h_counts = orth['human_Symbol'].value_counts()
m_counts = orth['mouse_Symbol'].value_counts()
me_counts = orth['mouse_EnsemblID'].value_counts()
orth11 = orth[orth['human_Symbol'].isin(h_counts[h_counts == 1].index) &
              orth['mouse_Symbol'].isin(m_counts[m_counts == 1].index) &
              orth['mouse_EnsemblID'].isin(me_counts[me_counts == 1].index)].drop_duplicates()
print('strict 1:1 ortholog pairs:', len(orth11))

# ---- mouse raw counts ----
raw = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/merged_unintegrated.h5ad')
april = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/IAHSC.scRNA.hsc_subclasses.April23.h5ad', backed='r')

mouse = raw[april.obs_names].copy()
mouse.obs['HSC_subclass'] = april.obs['HSC_subclass'].values
mouse.obs['sample'] = april.obs['sample'].values
mouse.obs['Genotype'] = ['Dnmt3a-KO' if 'Dnmt3a-KO' in s else 'WT' for s in mouse.obs['sample']]
mouse.obs['infection_status'] = ['infected' if 'infected' in s else 'naive' for s in mouse.obs['sample']]

# WT only, drop Erythroid
mouse = mouse[(mouse.obs['Genotype'] == 'WT') & (mouse.obs['HSC_subclass'] != 'Erythroid')].copy()
print('mouse WT non-erythroid cells:', mouse.n_obs)
print(mouse.obs['HSC_subclass'].value_counts())

mouse_genes_avail = set(mouse.var_names)
orth_m = orth11[orth11['mouse_EnsemblID'].isin(mouse_genes_avail)]
print('orthologs found in mouse data:', len(orth_m))

mouse_sub = mouse[:, orth_m['mouse_EnsemblID'].values].copy()
# rename var to human symbol (canonical cross-species identifier)
mouse_sub.var_names = orth_m['human_Symbol'].values
mouse_sub.var_names_make_unique()
print('mouse_sub:', mouse_sub.shape, 'X dtype', mouse_sub.X.dtype, 'max', mouse_sub.X.max())
mouse_sub.write(f'{OUT}/mouse_ortholog_raw.h5ad')

# ---- human raw counts ----
human = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/h5ad/BCG_scRNA_HSC_annotated.h5ad')
human = human[human.obs['clust_names'].isin(['HSC_a', 'HSC_b'])].copy()
print('human HSC cells:', human.n_obs)
human.X = human.layers['raw_counts'].copy()

human_genes_avail = set(human.var_names)
orth_h = orth11[orth11['human_Symbol'].isin(human_genes_avail)]
# restrict to orthologs also found in mouse (so the two matrices share identical gene universe)
orth_h = orth_h[orth_h['mouse_EnsemblID'].isin(mouse_genes_avail)]
print('orthologs found in human data (and also in mouse):', len(orth_h))

human_sub = human[:, orth_h['human_Symbol'].values].copy()
human_sub.var_names_make_unique()
print('human_sub:', human_sub.shape, 'X dtype', human_sub.X.dtype, 'max', human_sub.X.max())
human_sub.write(f'{OUT}/human_ortholog_raw.h5ad')

# also restrict mouse_sub to exactly this same shared gene universe order for later comparability
shared_genes = sorted(set(mouse_sub.var_names) & set(human_sub.var_names))
print('final shared ortholog gene universe:', len(shared_genes))
with open(f'{OUT}/shared_ortholog_genes.txt', 'w') as f:
    f.write('\n'.join(shared_genes))

mouse_sub.obs[['HSC_subclass', 'sample', 'infection_status']].to_csv(f'{OUT}/mouse_obs.csv')
human_sub.obs[['hsc_clusters_lab', 'clust_names', 'vaccination', 'timepoint', 'BEST']].to_csv(f'{OUT}/human_obs.csv')
print('done')
