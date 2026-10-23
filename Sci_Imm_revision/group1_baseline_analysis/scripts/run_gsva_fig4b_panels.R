suppressMessages(library(GSVA))

OUT <- "/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis"

counts <- read.csv(file.path(OUT, "scRNA_pseudobulk_raw_counts_symbols.csv"), row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
cat("counts dim:", dim(counts), "\n")

cpm <- t(t(counts) / colSums(counts)) * 1e6
log2cpm <- log2(cpm + 1)

# manuscript's own Fig. 4B leading-edge gene panels, read directly off the published figure
gene_sets <- list(
  ResponseToBacterium_Fig4B    = c("Mpeg1", "Plac8", "Tent5a", "Mt2", "Tifab"),
  MyeloidCellActivation_Fig4B  = c("Myo1f", "Ifngr1", "Sphk1", "Cd48", "Cebpa", "Cst7", "Tyrobp", "Rab44", "Ms4a2", "Fcgr3"),
  MonocyteChemotaxis_Fig4B     = c("Flt1", "Ccl9", "Ccr1", "Ctsg", "Cx3cr1"),
  RegulationCytokineResponse_Fig4B = c("Alox5", "Mefv", "Elane", "Il17rb", "Cd6"),
  TypeIIIFNProduction_Fig4B    = c("Irf8", "Lgals9", "Pglyrp1", "Klrk1", "Ccr2")
)

for (nm in names(gene_sets)) {
  present <- intersect(gene_sets[[nm]], rownames(log2cpm))
  cat(nm, ":", length(present), "/", length(gene_sets[[nm]]), "genes present\n")
}

scores <- gsva(log2cpm, gene_sets, method = "gsva", verbose = FALSE)

out_df <- as.data.frame(t(scores))
write.csv(out_df, file.path(OUT, "scRNA_pseudobulk_gsva_scores_Fig4B_panels.csv"))
cat("saved\n")
print(out_df)
