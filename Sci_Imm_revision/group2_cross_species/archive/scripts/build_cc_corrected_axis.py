import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import LinearRegression

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
latent = ad.read_h5ad(f'{OUT}/embed_latent50.h5ad')

mouse_mask = (latent.obs['system'] == 'mouse').values
human_mask = (latent.obs['system'] == 'human').values

# --- regress cell-cycle (S_score, G2M_score) out of each latent dim, using MOUSE cells only ---
mouse_X = latent.X[mouse_mask]
cc = latent.obs.loc[mouse_mask, ['S_score', 'G2M_score']].values.astype(float)

lr = LinearRegression()
lr.fit(cc, mouse_X)
mouse_X_pred = lr.predict(cc)
mouse_X_resid = mouse_X - mouse_X_pred

mouse_obs = latent.obs.loc[mouse_mask].reset_index(drop=True)
ia_mask = (mouse_obs['HSC_subclass'] == 'IA-HSC').values
lt_mask = (mouse_obs['HSC_subclass'] == 'LT-HSC').values

d_new = mouse_X_resid[ia_mask].mean(axis=0) - mouse_X_resid[lt_mask].mean(axis=0)
d_new = d_new / np.linalg.norm(d_new)

d_old = latent.X[mouse_mask][ia_mask].mean(axis=0) - latent.X[mouse_mask][lt_mask].mean(axis=0)
d_old = d_old / np.linalg.norm(d_old)

print('cosine similarity between old and new (CC-corrected) axis:', np.dot(d_old, d_new))

# project ALL cells (original, uncorrected latent coords) onto the new axis direction
proj_new = latent.X @ d_new
latent.obs['proj_new'] = proj_new
latent.obs['proj_old'] = latent.X @ d_old

# --- validate: does the new axis still separate mouse subclasses sensibly? ---
print('\n=== mouse: new axis projection by HSC_subclass ===')
print(latent.obs.loc[mouse_mask].groupby('HSC_subclass', observed=True)['proj_new'].mean().sort_values(ascending=False))

# --- validate: is the new axis no longer confounded with cell cycle? ---
sub = latent.obs.loc[mouse_mask]
r_g2m_new, p_g2m_new = spearmanr(sub['proj_new'], sub['G2M_score'])
r_s_new, p_s_new = spearmanr(sub['proj_new'], sub['S_score'])
r_g2m_old, p_g2m_old = spearmanr(sub['proj_old'], sub['G2M_score'])
r_s_old, p_s_old = spearmanr(sub['proj_old'], sub['S_score'])
print(f'\nOLD axis vs G2M_score: r={r_g2m_old:.3f} (p={p_g2m_old:.1e})   vs S_score: r={r_s_old:.3f} (p={p_s_old:.1e})')
print(f'NEW axis vs G2M_score: r={r_g2m_new:.3f} (p={p_g2m_new:.1e})   vs S_score: r={r_s_new:.3f} (p={p_s_new:.1e})')

# --- validate: does the new axis still recover c3/c4/c7 in human data? ---
print('\n=== human: new axis projection by hsc_clusters_lab ===')
print(latent.obs.loc[human_mask].groupby('hsc_clusters_lab', observed=True)['proj_new'].mean().sort_values(ascending=False))
print('\n=== human: OLD axis projection by hsc_clusters_lab (for comparison) ===')
print(latent.obs.loc[human_mask].groupby('hsc_clusters_lab', observed=True)['proj_old'].mean().sort_values(ascending=False))

latent.obs[['proj_new', 'proj_old']].to_csv(f'{OUT}/axis_projection_cc_corrected.csv')
np.save(f'{OUT}/axis_vector_cc_corrected.npy', d_new)
print('\nsaved axis_projection_cc_corrected.csv, axis_vector_cc_corrected.npy')
