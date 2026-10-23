import pandas as pd
import numpy as np
import gseapy as gp

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'

df = pd.read_csv(f'{OUT}/integration_driver_genes_full.csv')
df = df.dropna(subset=['mouse_corr', 'human_corr'])
df['score'] = df['mouse_corr'] * df['human_corr']
rnk = df[['gene', 'score']].sort_values('score', ascending=False).reset_index(drop=True)
rnk = rnk.drop_duplicates('gene')
print(f'{len(rnk)} genes ranked. score range: {rnk["score"].min():.3f} to {rnk["score"].max():.3f}')
print(rnk.head(10))
print(rnk.tail(10))

pre_res = gp.prerank(
    rnk=rnk,
    gene_sets='MSigDB_Hallmark_2020',
    organism='mouse',
    min_size=5,
    max_size=1000,
    permutation_num=1000,
    seed=0,
    threads=4,
    outdir=None,
)

res = pre_res.res2d.copy()
res['FDR q-val'] = res['FDR q-val'].astype(float)
res['NES'] = res['NES'].astype(float)
res = res.sort_values('FDR q-val')
res.to_csv(f'{OUT}/gsea_prerank_hallmark_crossspecies.csv', index=False)
print(res[['Term', 'NES', 'NOM p-val', 'FDR q-val']].to_string())
print(f'\n{(res["FDR q-val"] < 0.25).sum()} terms at FDR<0.25 (standard GSEA threshold)')
print(f'{(res["FDR q-val"] < 0.05).sum()} terms at FDR<0.05')
