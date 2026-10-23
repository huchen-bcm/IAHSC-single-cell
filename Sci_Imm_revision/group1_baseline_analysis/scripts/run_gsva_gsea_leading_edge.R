suppressMessages(library(GSVA))

OUT <- "/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis"

counts <- read.csv(file.path(OUT, "scRNA_pseudobulk_raw_counts_symbols.csv"), row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
cpm <- t(t(counts) / colSums(counts)) * 1e6
log2cpm <- log2(cpm + 1)

# case-insensitive lookup table: uppercase symbol -> actual matrix rowname
sym_upper <- toupper(rownames(log2cpm))
names(sym_upper) <- rownames(log2cpm)
lookup <- setNames(rownames(log2cpm), sym_upper)

# leading-edge gene sets, extracted directly from the actual GSEA prerank output
# (bulk_IAHSC_vs_LTHSC_gsea/IAHSC_vs_LTHSC.bulk.gsea.GO_Biological_Process_2025.csv, Lead_genes column)
# for the best-matching positive-NES term per category
lead_edge_upper <- list(
  TypeII_IFN_LeadingEdge = c('IFNGR1','LGALS9','IRF8','GAPDH','SYNCRIP','CD47','HCK','DAPK3','LONP1','PDE12','CCL3','CASP1'),
  LeukocyteChemotaxis_LeadingEdge = c('CCR2','C1QBP','LGALS9','CALR','CCR1','HMGB1'),
  MacrophageActivation_LeadingEdge = c('CEBPA','SPHK1','FCGR2B','HSPD1','IL31RA','MCUB','TLR6','KARS1'),
  ResponseToBacterium_LeadingEdge = c('MPEG1','H2BC11','MPO','OPTN','CTSG','PGLYRP1','HMGB2','ELANE','TLR9','JCHAIN','NLRC4','H2BC12','FCER1G','H2BC4','CHAMP1','FAU','TLR6','SLPI','ANXA3','CASP1','TNFRSF14','RPS19','DROSHA','LYST','PRKCD','RASGRP4'),
  CytokineProduction_LeadingEdge = c('CCR2','ERMAP','DTX4','FCGR2B','LGALS9','IRF8','RNF26','CASP3','C5AR2','GAPDH','MAPKAPK2','CD47','CARD9','RGCC','ZBTB1','SOD1','CLEC9A','SRGN','UFC1','NFAM1')
)

gene_sets <- list()
for (nm in names(lead_edge_upper)) {
  mapped <- lookup[lead_edge_upper[[nm]]]
  mapped <- mapped[!is.na(mapped)]
  gene_sets[[nm]] <- unname(mapped)
  cat(nm, ':', length(mapped), '/', length(lead_edge_upper[[nm]]), 'genes matched in expression matrix\n')
}

scores <- gsva(log2cpm, gene_sets, method = "gsva", verbose = FALSE)
out_df <- as.data.frame(t(scores))
write.csv(out_df, file.path(OUT, "scRNA_pseudobulk_gsva_scores_GSEA_leading_edge.csv"))
cat("saved\n")
print(out_df)
