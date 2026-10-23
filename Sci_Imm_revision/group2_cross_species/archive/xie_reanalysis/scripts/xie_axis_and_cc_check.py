import scanpy as sc
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

XIE = '/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/4.analysis/human_data/Xie'

obs = pd.read_csv('/tmp/xie_merged_obs.csv', index_col=0)
embed_all = np.load('/tmp/xie_embed_all100.npy')
embed30 = embed_all[:, :30]  # matches the notebook's own convention

mouse_mask = (obs['system'] == 'mouse').values
human_mask = (obs['system'] == 'human').values

# --- build axis: mouse IA-HSC centroid - mouse LT-HSC centroid, using first 30 latent dims ---
ia_mask = (obs['HSC_subclass'] == 'IA-HSC').values
lt_mask = (obs['HSC_subclass'] == 'LT-HSC').values
d = embed30[ia_mask].mean(axis=0) - embed30[lt_mask].mean(axis=0)
d = d / np.linalg.norm(d)
proj = embed30 @ d
obs['proj'] = proj

print('=== human (Xie): axis projection by Cluster ===')
print(obs.loc[human_mask].groupby('Cluster', observed=True)['proj'].mean().sort_values(ascending=False))

print('\n=== mouse: axis projection by HSC_subclass ===')
print(obs.loc[mouse_mask].groupby('HSC_subclass', observed=True)['proj'].mean().sort_values(ascending=False))

# --- immediately check cell-cycle confound ---
mouse_cc = sc.read_h5ad('/home/huc/projects/Ruoqiong_Cao.IAHSC/scRNA/3.annotation/processed_anndata/IAHSC.scRNA.hsc_subclasses.April23.h5ad',
                         backed='r').obs[['S_score', 'G2M_score', 'phase']]

mouse_obs_names = obs.index[mouse_mask]
common = mouse_obs_names.intersection(mouse_cc.index)
print(f'\nmouse cells in merged obj: {mouse_mask.sum()}, matched to CC-scored object: {len(common)}')

sub = obs.loc[common].copy()
sub['S_score'] = mouse_cc.loc[common, 'S_score'].values
sub['G2M_score'] = mouse_cc.loc[common, 'G2M_score'].values
sub['phase'] = mouse_cc.loc[common, 'phase'].values

print('\n=== mouse: phase composition by HSC_subclass ===')
print(sub.groupby('HSC_subclass', observed=True)['phase'].value_counts(normalize=True).unstack())

r_g2m, p_g2m = spearmanr(sub['proj'], sub['G2M_score'])
r_s, p_s = spearmanr(sub['proj'], sub['S_score'])
print(f'\naxis proj vs G2M_score: r={r_g2m:.3f} (p={p_g2m:.1e})')
print(f'axis proj vs S_score:   r={r_s:.3f} (p={p_s:.1e})')

sub.to_csv('/tmp/xie_mouse_axis_cc.csv')
