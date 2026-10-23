#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

.libPaths("/project/lbarreiro/USERS/sarah/Rlibs_new")

install.packages("ggbeeswarm")
library(stringr)
library(textTinyR)
library(pbapply)
library(readr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(plyr)
library(ggbeeswarm)

setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/")
OUT_DIR <- "emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/"
IN_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/scatterplot_data_individual/"

############################
## Combo of HSCa and HSCb ##
############################

## read in regulon data
REGULON_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic/get_diff_regulons_outs/"
regulon_HSCa <- read.csv(file=paste0(REGULON_DIR,"HSC_a_wilcox_data_no_outliers.csv"))
regulon_HSCb <- read.csv(file=paste0(REGULON_DIR,"HSC_b_wilcox_data_no_outliers.csv"))
regulons_final <- unique(c(regulon_HSCa$X, regulon_HSCb$X))

# ## read in TF expression data 
# TF_expression_dir <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/MASH_emmreml_downstream/DE_overlapping_TFs/"
# FULL_TF_data <- read.csv(file=paste0(TF_expression_dir, "TF_exp_FULL_data.csv"))
# HSCa_exp_data <- subset(FULL_TF_data, FULL_TF_data$celltype == "HSC_a")
# 
# 
# ## Create a meta pval for each regulon 
# 
# 
# xi <- -2*(log(pval_1) + log(pval_2))
# multivariate_p <- 1 - pchisq(xi, 4)


clusters <- c("CMP", "GMP")

for(i in 1:length(clusters)){
  ## read in motif enrichment data 
  enrichment_data <- read.table(file=paste0("emmreml_downstream/HINT_motif_enrichment_analysis/hint_Enrichment/",clusters[i],"_homer_DE_peaks/fulltest_statistics.txt"), header=TRUE, fill=TRUE)
  
  str_sub(enrichment_data$FACTOR, start = 1, end = 9, omit_na = FALSE) <- ""
  
  
  pathway_TFs <- intersect(enrichment_data$FACTOR, regulons_final)
  motif_data_final <- subset(enrichment_data, enrichment_data$FACTOR %in% pathway_TFs)
  
  final_pathway <- motif_data_final$FACTOR[which(motif_data_final$CORR.P.VALUE < 0.05)]
  
  
  ## order regulon TFs by p-value 
  regulon_TFs <- regulons_final[which(regulons_final %in% pathway_TFs)]
  regulon_value <- vector(length = length(regulon_TFs))
  for(k in 1:length(regulon_TFs)){
    if(regulon_TFs[k] %in% regulon_HSCa$X){
      if(regulon_TFs[k] %in% regulon_HSCb$X){
        regulon_value[k] <- min(regulon_HSCa$pval_wilcox[which(regulon_HSCa$X == regulon_TFs[k])],
                                regulon_HSCb$pval_wilcox[which(regulon_HSCb$X == regulon_TFs[k])])
      }
      else{
        regulon_value[k] <-regulon_HSCa$pval_wilcox[which(regulon_HSCa$X == regulon_TFs[k])]
      }
    }
    else{
      regulon_value[k] <-regulon_HSCb$pval_wilcox[which(regulon_HSCb$X == regulon_TFs[k])]
    }
  }
  
  
  rank_vals <- -log10(regulon_value)
  names(rank_vals) <- regulon_TFs
  sorted_ranked_list <- sort(rank_vals, decreasing=FALSE)
  assign(paste0(name,"_table"), sorted_ranked_list)
  
  pathways_list <- list()
  pathways_list[[1]] <- as.character(final_pathway)
  names(pathways_list) <- c("motif_enriched")
  
  name <- clusters[i]
  
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
  
  write.csv(gsea_result_df, file=paste0(OUT_DIR,"gsea_regulon/",name,"_gsea_result"))
  print(i)
  
  file_name = paste0(OUT_DIR, paste0("gsea_regulon/CMP_gsea.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plotEnrichment(pathways_list[["motif_enriched"]],
                       sorted_ranked_list))
  dev.off()
}


###############
## HSCa only ##
###############

## read in regulon data
REGULON_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic/get_diff_regulons_outs/"
regulon_HSCa <- read.csv(file=paste0(REGULON_DIR,"HSC_a_wilcox_data_no_outliers.csv"))
regulons_final <- unique(c(regulon_HSCa$X))



clusters <- c("CMP", "GMP")

for(i in 1:length(clusters)){
  name <- clusters[i]
  ## read in motif enrichment data 
  enrichment_data <- read.table(file=paste0("emmreml_downstream/HINT_motif_enrichment_analysis/hint_Enrichment/",clusters[i],"_homer_DE_peaks/fulltest_statistics.txt"), header=TRUE, fill=TRUE)
  
  str_sub(enrichment_data$FACTOR, start = 1, end = 9, omit_na = FALSE) <- ""
  
  
  pathway_TFs <- intersect(enrichment_data$FACTOR, regulons_final)
  motif_data_final <- subset(enrichment_data, enrichment_data$FACTOR %in% pathway_TFs)
  
  final_pathway <- motif_data_final$FACTOR[which(motif_data_final$CORR.P.VALUE < 0.01)]
  
  
  ## order regulon TFs by p-value 
  regulon_TFs <- regulons_final[which(regulons_final %in% pathway_TFs)]
  regulon_data_final <- regulon_HSCa[which(regulon_HSCa$X %in% regulon_TFs),]
  
  rank_vals <- -log10(regulon_data_final$pval_wilcox)
  names(rank_vals) <- regulon_data_final$X
  sorted_ranked_list <- sort(rank_vals, decreasing=FALSE)
  assign(paste0(name,"_table"), sorted_ranked_list)
  
  pathways_list <- list()
  pathways_list[[1]] <- as.character(final_pathway)
  names(pathways_list) <- c("motif_enriched")
  
  
  
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
  
  write.csv(gsea_result_df, file=paste0(OUT_DIR,"gsea_regulon/",name,"_gsea_result_HSCa_0.01"))
  print(i)
  
  file_name = paste0(OUT_DIR, paste0("gsea_regulon/CMP_gsea_HSCa_0.01.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plotEnrichment(pathways_list[["motif_enriched"]],
                       sorted_ranked_list))
  dev.off()
}




################################################
## HSCa only regulon p.value/percentage score ##
################################################

## read in regulon data
REGULON_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic/get_diff_regulons_outs/"
regulon_HSCa <- read.csv(file=paste0(REGULON_DIR,"HSC_a_wilcox_data_no_outliers.csv"))
regulons_final <- unique(c(regulon_HSCa$X))



clusters <- c("CMP", "GMP")

for(i in 1:length(clusters)){
  ## read in motif enrichment data 
  name <- clusters[i]
  enrichment_data <- read.csv(file=paste0(IN_DIR, "0.1_scatterplot_data_", name , ".csv"))
  
  str_sub(enrichment_data$name, start = 1, end = 9, omit_na = FALSE) <- ""
  
  
  pathway_TFs <- intersect(enrichment_data$name, regulons_final)
  motif_data_final <- subset(enrichment_data, enrichment_data$name %in% pathway_TFs)
  
  final_pathway1 <- motif_data_final[(which(motif_data_final$peak_fdr > -log10(0.1))),]
  final_pathway <- final_pathway1[(which(final_pathway1$perc > 0.1)),]$name
  
  
  ## order regulon TFs by p-value 
  regulon_TFs <- regulons_final[which(regulons_final %in% pathway_TFs)]
  regulon_data_final <- regulon_HSCa[which(regulon_HSCa$X %in% regulon_TFs),]
  
  rank_vals <- -log10(regulon_data_final$pval_wilcox)
  names(rank_vals) <- regulon_data_final$X
  sorted_ranked_list <- sort(rank_vals, decreasing=FALSE)
  assign(paste0(name,"_table"), sorted_ranked_list)
  
  pathways_list <- list()
  pathways_list[[1]] <- as.character(final_pathway)
  names(pathways_list) <- c("motif_enriched")
  
  
  
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
  
  write.csv(gsea_result_df, file=paste0(OUT_DIR,"gsea_regulon/",name,"_gsea_result"))
  print(i)
  
  file_name = paste0(OUT_DIR, paste0("gsea_regulon/GMP_gsea.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plotEnrichment(pathways_list[["motif_enriched"]],
                       sorted_ranked_list))
  dev.off()
}
