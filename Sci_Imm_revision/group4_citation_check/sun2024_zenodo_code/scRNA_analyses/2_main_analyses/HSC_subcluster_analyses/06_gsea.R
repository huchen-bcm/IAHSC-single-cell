
library("ggplot2")
library ("DESeq2")
library(statmod)
library(goseq)
library("RColorBrewer")
library(edgeR)
library(devtools)
library(cluster)
library(fgsea)
library(data.table)
library(ggplot2)
library(msigdbr)

# We perform GSEA for HSC subclusters 


# setup ------------------------------------------------------
set.seed(831)
setwd("/scRNA_analyses/2_main_analyses/HSC_subcluster_analyses")
OUT_DIR <- c("9_gsea/")
dir.create(OUT_DIR)


# Read in posterior means without correlations and lfsr -----------------------------------------

PM_correlations <- read.table(file=paste0("5_mash/mash_results/posteriorMeans.txt"), header=TRUE)
lfsr_correlations <- read.table(file=paste0("5_mash/mash_results/lfsr_output.txt"), header=TRUE)

clusters <- colnames(PM_correlations)
for(i in 1:length(clusters)){
  name <- clusters[i]
  betas <- PM_correlations[,paste0(clusters[i])]
  lfsr_vals <- lfsr_correlations[,paste0(clusters[i])]
  rank_vals <- -log10(lfsr_vals)*betas
  names(rank_vals) <- row.names(PM_correlations)
  sorted_ranked_list <- sort(rank_vals, decreasing=FALSE)
  assign(paste0(name,"_table"), sorted_ranked_list)
  print(i)
}

df_list <- list(c0_table, c1_table, c2_table, c3_table, c4_table, c5_table, c6_table, c7_table, c8_table)



# Gene set enrichment of Hallmark pathways --------------------------------

for(i in 1:length(df_list)){
  name <- clusters[i]
  sorted_ranked_list <- df_list[[i]]
  msigdbr_show_species()
  h_df = msigdbr(species = "Homo sapiens", category = "H")  #Hallmark categories n=50
  h_list = h_df %>% split(x = .$gene_symbol, f = .$gs_name)
  
  fgseaRes <- fgsea(pathways = h_list, stats= sorted_ranked_list,
                    maxSize  = 500, nperm=100000)
  
  padjusted <- p.adjust(fgseaRes$pval, method='BH')
  
  gsea_result_df<- data.frame(
    pathways <- fgseaRes$pathway,
    pval <- fgseaRes$pval,
    NES <- fgseaRes$NES,
    padj <- padjusted
  )
  colnames(gsea_result_df) <- c("pathways","pval","NES", "padj")
  write.csv(gsea_result_df, file=paste0(OUT_DIR,"NOT_SORTED_gsea_result_",name,"_no_correlations"))
  gsea_result_df <- gsea_result_df[order(gsea_result_df$pval), ]
  assign(paste0(name,"_gsea_result"), gsea_result_df)
  write.csv(gsea_result_df, file=paste0(OUT_DIR,"SORTED_gsea_result_",name, "_no_correlations"))
  print(i)
}
