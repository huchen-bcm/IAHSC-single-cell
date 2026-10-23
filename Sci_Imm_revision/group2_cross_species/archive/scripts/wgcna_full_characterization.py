import anndata as ad
import numpy as np
import pandas as pd
import gseapy as gp
import time

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

latent = ad.read_h5ad(f'{OUT}/embed_latent50.h5ad')
train = ad.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/BCG/adata_merged_human_mouse_orthologs_hvg_integrated.0128_2026.h5ad')
assert (latent.obs_names == train.obs_names).all()

X = train.X
if hasattr(X, 'toarray'):
    X = X.toarray()
genes = train.var_names.to_numpy()
gene_idx = {g: i for i, g in enumerate(genes)}

mouse_mask = (latent.obs['system'] == 'mouse').values
human_mask = (latent.obs['system'] == 'human').values
mouse_subclass = latent.obs['HSC_subclass'].astype(str).values
human_cluster = latent.obs['hsc_clusters_lab'].astype(str).values

mouse_mod = pd.read_csv(f'{OUT}/wgcna_mouse_gene_modules.csv', index_col=0)
human_mod = pd.read_csv(f'{OUT}/wgcna_human_gene_modules.csv', index_col=0)

mouse_modules = {c: mouse_mod[mouse_mod['moduleColors'] == c].index.tolist() for c in mouse_mod['moduleColors'].unique()}
human_modules = {c: human_mod[human_mod['moduleColors'] == c].index.tolist() for c in human_mod['moduleColors'].unique()}

all_modules = {}
for c, g in mouse_modules.items():
    all_modules[f'mouse_{c}'] = [gg for gg in g if gg in gene_idx]
for c, g in human_modules.items():
    all_modules[f'human_{c}'] = [gg for gg in g if gg in gene_idx]

print('module sizes:', {k: len(v) for k, v in all_modules.items()})

# ---- 1. module score per cell = mean expression of module genes (both species scored on ALL modules) ----
scores = np.zeros((X.shape[0], len(all_modules)), dtype=np.float32)
mod_names = list(all_modules.keys())
for j, name in enumerate(mod_names):
    idxs = [gene_idx[g] for g in all_modules[name]]
    scores[:, j] = X[:, idxs].mean(axis=1)
score_df = pd.DataFrame(scores, columns=mod_names, index=train.obs_names)
score_df.to_csv(f'{OUT}/wgcna_all_module_scores_percell.csv')
print('saved per-cell module scores', score_df.shape)

# ---- 2. mean module score by mouse subclass and human cluster (own AND cross-species modules) ----
mouse_order = ['IA-HSC', 'LT-HSC', 'ST-HSC', 'Cd41+ HSC']
human_order = [f'c{i}' for i in range(10)]

mouse_profile = pd.DataFrame(index=mouse_order, columns=mod_names, dtype=float)
for g in mouse_order:
    sel = mouse_mask & (mouse_subclass == g)
    mouse_profile.loc[g] = score_df.loc[sel].mean(axis=0)

human_profile = pd.DataFrame(index=human_order, columns=mod_names, dtype=float)
for g in human_order:
    sel = human_mask & (human_cluster == g)
    human_profile.loc[g] = score_df.loc[sel].mean(axis=0)

mouse_profile.to_csv(f'{OUT}/wgcna_mouse_subtype_module_profile.csv')
human_profile.to_csv(f'{OUT}/wgcna_human_cluster_module_profile.csv')
print('saved subtype/cluster module profiles')

# ---- 3. cross-species subtype correlation using the shared (all_modules) score feature space, z-scored per module ----
combined = pd.concat([mouse_profile, human_profile])
z = (combined - combined.mean()) / combined.std()
mouse_z = z.loc[mouse_order]
human_z = z.loc[human_order]

corr = pd.DataFrame(index=mouse_order, columns=human_order, dtype=float)
for m in mouse_order:
    for h in human_order:
        corr.loc[m, h] = np.corrcoef(mouse_z.loc[m].astype(float), human_z.loc[h].astype(float))[0, 1]
corr.to_csv(f'{OUT}/wgcna_crossspecies_subtype_correlation.csv')
print(corr.round(3))

# ---- 4. GO enrichment for every module ----
all_go = {}
for name, glist in all_modules.items():
    if len(glist) < 5:
        continue
    for attempt in range(3):
        try:
            enr = gp.enrichr(gene_list=glist, gene_sets=['GO_Biological_Process_2025'], organism='mouse', outdir=None)
            df = enr.results
            df = df[df['Adjusted P-value'] < 0.05].sort_values('Adjusted P-value')
            all_go[name] = df['Term'].head(6).tolist()
            print(name, len(glist), 'genes ->', all_go[name][:3])
            break
        except Exception as e:
            print(name, 'retry', attempt, e)
            time.sleep(4)
    time.sleep(1)

import json
with open(f'{OUT}/wgcna_all_modules_GO_top_terms.json', 'w') as f:
    json.dump(all_go, f, indent=1)
print('all done')
