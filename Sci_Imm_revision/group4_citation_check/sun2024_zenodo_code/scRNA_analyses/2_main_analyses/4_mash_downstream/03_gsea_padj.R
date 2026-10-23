
#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

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

# In this code a GSEA is performed using betas and log2FC values for genes within each HSPC cluster --------------------

# setup -------------------------------------------------
set.seed(831)
setwd("/scRNA_analyses/2_main_analyses/4_mash_downstream")
OUT_DIR <- c("gsea_padj/")
dir.create(OUT_DIR)



# Read in posterior means and lfsr ----------
PM <- read.table(file=paste0("../MASH/mash_results/posteriorMeans_wcustom.txt"), header=TRUE)
lfsr <- read.table(file=paste0("../MASH/mash_results/lfsr_wcustom_output.txt"), header=TRUE)

clusters <- colnames(PM)
for(i in 1:length(clusters)){
  name <- clusters[i]
  betas <- PM[,paste0(clusters[i])]
  lfsr_vals <- lfsr[,paste0(clusters[i])]
  rank_vals <- -log10(lfsr_vals)*betas
  names(rank_vals) <- row.names(PM)
  sorted_ranked_list <- sort(rank_vals, decreasing=FALSE)
  assign(paste0(name,"_table"), sorted_ranked_list)
  print(i)
}

df_list <- list(HSC_a_table, HSC_b_table, CMP_a_table, CMP_b_table, CMP_c_table, GMP_b1_table, GMP_b2_table, MEP_a1_table, MEP_a2_table, MEP_a3_table, mixed_a_table, mixed_b_table, PreBNK_table)



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
    padj <- fgseaRes$padj
  )
  colnames(gsea_result_df) <- c("pathways","pval","NES", "padj")
  #write.csv(gsea_result_df, file=paste0(OUT_DIR,"NOT_SORTED_gsea_result_",name,"_NO_correlations_with_custom"))
  gsea_result_df <- gsea_result_df[order(gsea_result_df$pval), ]
  assign(paste0(name,"_gsea_result"), gsea_result_df)
  #write.csv(gsea_result_df, file=paste0(OUT_DIR,"SORTED_gsea_result_",name, "_NO_correlations_with_custom"))
  print(i)
  
  
  ##save for paper Table 3
  upd_names <- c("HSC c1", "HSC c2", "CMP c1", "CMP c2", "CMP c3", "GMP c1", "GMP c2", "MEP c1", "MEP c2",
                 "MEP c3", "MLP c1", "MLP c2", "PreBNK")
  write.csv(gsea_result_df, file=paste0(OUT_DIR,"copy_for_paper/SORTED_gsea_result_",upd_names[i],".csv"))
  
}
