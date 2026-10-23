

#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

.libPaths("/project/lbarreiro/USERS/sarah/Rlibs_new")

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

set.seed(831)

setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/")
IN_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/scatterplot_data/"
OUT_DIR <- "gsea/"
dir.create(OUT_DIR)


clusters <- c("CMP", "GMP")
for(i in 1:length(clusters)){
  name <- clusters[i]
  
  cluster_data <- read.csv(file=paste0(IN_DIR, "scatterplot_data_", name , ".csv"))
  percent <- cluster_data[,"percent"]
  p_vals <- cluster_data[,"neglogfdr"]
  rank_vals <- p_vals*percent
  names(rank_vals) <- cluster_data[,"TF_num"]
  sorted_ranked_list <- sort(rank_vals, decreasing=FALSE)
  assign(paste0(name,"_table"), sorted_ranked_list)
  print(i)
}

df_list <- list(CMP_table, GMP_table)

##############################################################
## Get list of TFs that are active in HSC p < 0.15 criteria ##
##############################################################

regulon_data <- read.csv(file=paste0(IN_DIR, "scatterplot_data_CMP.csv")) ##using CMP but same for GMP
HSC_active <- regulon_data$TF_num[which(regulon_data$regulon_combined == "HSC active")]
pathways_list <- list()
pathways_list[[1]] <- as.character(HSC_active)
names(pathways_list) <- c("HSC_active")


for(i in 1:length(df_list)){
  name <- clusters[i]
  sorted_ranked_list <- df_list[[i]]

  fgseaRes <- fgsea(pathways = pathways_list, stats= sorted_ranked_list,
                    maxSize  = 500, nperm=100000)
  
  padjusted <- p.adjust(fgseaRes$pval, method='BH')
  
  gsea_result_df<- data.frame(
    pathways <- fgseaRes$pathway,
    pval <- fgseaRes$pval,
    NES <- fgseaRes$NES,
    padj <- fgseaRes$padj
  )
  colnames(gsea_result_df) <- c("pathways","pval","NES", "padj")
  write.csv(gsea_result_df, file=paste0(OUT_DIR,name,"_gsea_result"))
  print(i)
  
  file_name = paste0(OUT_DIR, paste0("GMP_gsea.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plotEnrichment(pathways_list[["HSC_active"]],
                 sorted_ranked_list))
  dev.off()
}



###############################################################
## Get list of TFs that are active in HSC p <= 0.05 criteria ##
###############################################################

regulon_data <- read.csv(file=paste0(IN_DIR, "0.05_scatterplot_data_CMP.csv")) ##using CMP but same for GMP
HSC_active <- regulon_data$TF_num[which(regulon_data$regulon_combined == "HSC active")]
pathways_list <- list()
pathways_list[[1]] <- as.character(HSC_active)
names(pathways_list) <- c("HSC_active")


for(i in 1:length(df_list)){
  name <- clusters[i]
  sorted_ranked_list <- df_list[[i]]
  
  fgseaRes <- fgsea(pathways = pathways_list, stats= sorted_ranked_list,
                    maxSize  = 500, nperm=100000)
  
  padjusted <- p.adjust(fgseaRes$pval, method='BH')
  
  gsea_result_df<- data.frame(
    pathways <- fgseaRes$pathway,
    pval <- fgseaRes$pval,
    NES <- fgseaRes$NES,
    padj <- fgseaRes$padj
  )
  colnames(gsea_result_df) <- c("pathways","pval","NES", "padj")
  write.csv(gsea_result_df, file=paste0(OUT_DIR,name,"_0.05_gsea_result"))
  print(i)
  
  file_name = paste0(OUT_DIR, paste0("GMP_gsea_0.05.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plotEnrichment(pathways_list[["HSC_active"]],
                       sorted_ranked_list))
  dev.off()

}




###############################################################
## Get list of TFs that are active in HSC p <= 0.1 criteria ##
###############################################################

regulon_data <- read.csv(file=paste0(IN_DIR, "0.1_scatterplot_data_CMP.csv")) ##using CMP but same for GMP
HSC_active <- regulon_data$TF_num[which(regulon_data$regulon_combined == "HSC active")]
pathways_list <- list()
pathways_list[[1]] <- as.character(HSC_active)
names(pathways_list) <- c("HSC_active")


for(i in 1:length(df_list)){
  name <- clusters[i]
  sorted_ranked_list <- df_list[[i]]
  
  fgseaRes <- fgsea(pathways = pathways_list, stats= sorted_ranked_list,
                    maxSize  = 500, nperm=100000)
  
  padjusted <- p.adjust(fgseaRes$pval, method='BH')
  
  gsea_result_df<- data.frame(
    pathways <- fgseaRes$pathway,
    pval <- fgseaRes$pval,
    NES <- fgseaRes$NES,
    padj <- fgseaRes$padj
  )
  colnames(gsea_result_df) <- c("pathways","pval","NES", "padj")
  write.csv(gsea_result_df, file=paste0(OUT_DIR,name,"_0.1_gsea_result"))
  print(i)
  
  file_name = paste0(OUT_DIR, paste0("GMP_gsea_0.1.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plotEnrichment(pathways_list[["HSC_active"]],
                       sorted_ranked_list))
  dev.off()
  
}



