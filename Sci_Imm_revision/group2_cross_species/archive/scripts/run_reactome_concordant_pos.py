import pandas as pd
import gseapy as gp

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group2_cross_species'
df = pd.read_csv(f'{OUT}/integration_driver_genes_concordant_positive.csv')
genes = df['gene'].tolist()
print(f'{len(genes)} genes (Quadrant I: positively concordant, IA-HSC/c3-4-7-high in both species)')

res = gp.enrichr(gene_list=genes, gene_sets=['Reactome_Pathways_2024'], organism='mouse', outdir=None)
out = res.results.sort_values('Adjusted P-value')
out.to_csv(f'{OUT}/integration_driver_Reactome_positive_IA_c347.csv', index=False)
print(out.head(20)[['Term', 'Overlap', 'Adjusted P-value']].to_string())
print(f"\n{(out['Adjusted P-value'] < 0.05).sum()} terms at FDR<0.05")
