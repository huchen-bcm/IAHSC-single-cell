import pandas as pd
import numpy as np
from scipy import stats

OUT = '/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis'
MEND = '/tmp/claude-186759/-home-huc-projects/b754521d-bab1-40e7-9d22-16b36231c512/scratchpad/mendeley_LTHSC'

# ---------- Task A: summarize Mendeley WT infected vs naive (their own Limma-Voom result) ----------
mend = pd.read_excel(f'{MEND}/WT_infect_vs_WT_naive_mendel.xls')
mend['ensembl_noversion'] = mend['Unnamed: 0'].str.split('.').str[0]
print("Mendeley table shape:", mend.shape)
print(mend[['Gene.name','logFC','P.Value','adj.P.Val']].describe())

sig = mend[mend['adj.P.Val'] < 0.05]
print(f"\nSignificant genes (adj.P.Val<0.05): {len(sig)} / {len(mend)}")
print(f"  Up in infected: {(sig['logFC']>0).sum()}, Down in infected: {(sig['logFC']<0).sum()}")

# named genes of interest from the manuscript / our pathway panels
genes_of_interest = ['Irf1','Irf7','Stat1','Gbp6','Gbp7','Ifngr1','Cx3cr1','Ccr1','Ccr2','Ccr5','Ccl9','Mpo','Cebpa']
print("\nManuscript-named genes in Mendeley LT-HSC infected-vs-naive:")
voi = mend[mend['Gene.name'].isin(genes_of_interest)][['Gene.name','logFC','P.Value','adj.P.Val']]
print(voi.to_string())

mend.to_csv(f'{OUT}/mendeley_LTHSC_WTinfected_vs_WTnaive.csv', index=False)

# ---------- Task B: cross-study concordance: our SLAM1-3 (infected LT-HSC, this study) vs Mendeley WT-infected (external) ----------
ours = pd.read_csv('/home/zhijiany/workdir/2025/040825_RNAseq_Ruoqiong/output/star_salmon/salmon.merged.gene_tpm.tsv', sep='\t')
ours['log2tpm_SLAM_mean'] = np.log2(ours[['SLAM1','SLAM2','SLAM3']].mean(axis=1) + 1)
ours['log2tpm_IA_mean'] = np.log2(ours[['IA1','IA2','IA3']].mean(axis=1) + 1)

mend['mend_infected_mean'] = mend[['WT infection_1','WT infection_2']].mean(axis=1)  # already log2-ish (voom logCPM)
mend['mend_naive'] = mend['WT naïve ']

merged = pd.merge(
    ours[['gene_name','log2tpm_SLAM_mean','log2tpm_IA_mean']],
    mend[['Gene.name','mend_infected_mean','mend_naive','logFC']],
    left_on='gene_name', right_on='Gene.name', how='inner'
).drop_duplicates(subset='gene_name')
print(f"\nShared genes for concordance check: {len(merged)}")

rho_infected, p_infected = stats.spearmanr(merged['log2tpm_SLAM_mean'], merged['mend_infected_mean'])
print(f"Concordance -- our SLAM1-3 (infected LT-HSC, this study) vs Mendeley WT-infected (external study): Spearman rho={rho_infected:.3f}, p={p_infected:.2e}, n={len(merged)} genes")

merged.to_csv(f'{OUT}/crossstudy_LTHSC_infected_concordance.csv', index=False)

with open(f'{OUT}/mendeley_LTHSC_summary.txt', 'w') as f:
    f.write(f"Mendeley WT LT-HSC infected-vs-naive: {len(sig)}/{len(mend)} genes significant (adj.P.Val<0.05); {(sig['logFC']>0).sum()} up, {(sig['logFC']<0).sum()} down in infected.\n\n")
    f.write("Manuscript-named genes:\n")
    f.write(voi.to_string() + "\n\n")
    f.write(f"Cross-study concordance (our SLAM1-3 vs Mendeley WT-infected), n={len(merged)} shared genes: Spearman rho={rho_infected:.3f}, p={p_infected:.2e}\n")

print("\nDone.")
