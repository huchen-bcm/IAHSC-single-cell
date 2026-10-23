suppressMessages(library(GSVA))

OUT <- "/home/huc/projects/Ruoqiong_Cao.IAHSC/Sci_Imm_revision/group1_baseline_analysis"

counts <- read.csv(file.path(OUT, "scRNA_pseudobulk_raw_counts_symbols.csv"), row.names = 1, check.names = FALSE)
counts <- as.matrix(counts)
cat("counts dim:", dim(counts), "\n")

cpm <- t(t(counts) / colSums(counts)) * 1e6
log2cpm <- log2(cpm + 1)

gmt_path <- "/home/huc/projects/Rachana.metabolome/new_data/m5.go.bp.v2023.2.Mm.symbols.gmt"
lines <- readLines(gmt_path)
target <- NULL
for (l in lines) {
  parts <- strsplit(l, "\t")[[1]]
  if (parts[1] == "GOBP_INNATE_IMMUNE_RESPONSE") {
    target <- parts[3:length(parts)]
    break
  }
}
cat("genes in GOBP_INNATE_IMMUNE_RESPONSE:", length(target), "\n")
present <- intersect(target, rownames(log2cpm))
cat("genes present in expression matrix:", length(present), "\n")

gene_sets <- list(InnateImmuneResponse_GOBP = target)

scores <- gsva(log2cpm, gene_sets, method = "gsva", verbose = FALSE)

out_df <- as.data.frame(t(scores))
write.csv(out_df, file.path(OUT, "scRNA_pseudobulk_gsva_scores_InnateImmuneResponse_GOBP.csv"))
cat("saved\n")
print(out_df)
