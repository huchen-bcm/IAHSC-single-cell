import pandas as pd
import gseapy as gp

df = pd.read_csv('/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/concordant_pos_r02.csv')
genes = df['gene'].tolist()
print(f'{len(genes)} genes (mouse_corr>0.2 AND human_corr>0.2):', genes)

res = gp.enrichr(gene_list=genes, gene_sets=['Reactome_Pathways_2024'], organism='mouse', outdir=None)
out = res.results.sort_values('Adjusted P-value')
OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
out.to_csv(f'{OUT}/integration_driver_Reactome_positive_r02.csv', index=False)
print(out.head(20)[['Term', 'Overlap', 'P-value', 'Adjusted P-value']].to_string())
print(f"\n{(out['Adjusted P-value'] < 0.05).sum()} terms at FDR<0.05")
