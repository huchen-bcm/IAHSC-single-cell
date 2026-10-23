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
IN_DIR <- "emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/"


# Cluster and TF info -----------------------------------------------
clusters <- c("CMP", "GMP", "MEP", "CLP", "PreBNK")
myColorScale <- c('#FF7F0E', '#33CC00', '#8C5782', '#46B8DA', '#1A0099')
plot_list1 <- list()
plot_list2 <- list()
plot_list3 <- list()
plot_list4 <- list()
plot_list5 <- list()
plot_list6 <- list()
TF_info <- read.csv(file=paste0("emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/TF_cluster_information_with_names.csv"))

TF_clusters <- unique(TF_info$clusters)


#Set sig regulons for each cluster group
REGULON_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic/get_diff_regulons_outs/"
regulon_CMPa <- read.csv(file=paste0(REGULON_DIR,"CMP_a_wilcox_data_no_outliers.csv"))
regulon_CMPb <- read.csv(file=paste0(REGULON_DIR,"CMP_b_wilcox_data_no_outliers.csv"))
regulon_CMPc <- read.csv(file=paste0(REGULON_DIR,"CMP_c_wilcox_data_no_outliers.csv"))
sig_regulons_cmp <- c(regulon_CMPa$X[which(regulon_CMPa$pval_wilcox < 0.15)], regulon_CMPb$X[which(regulon_CMPb$pval_wilcox < 0.15)], regulon_CMPc$X[which(regulon_CMPc$pval_wilcox < 0.15)])
sig_regulons_cmp_final <- unique(sig_regulons_cmp)

regulon_GMPa <- read.csv(file=paste0(REGULON_DIR,"GMP_a_wilcox_data_no_outliers.csv"))
regulon_GMPb <- read.csv(file=paste0(REGULON_DIR,"GMP_b_wilcox_data_no_outliers.csv"))
sig_regulons_gmp <- c(regulon_GMPa$X[which(regulon_GMPa$pval_wilcox < 0.15)], regulon_GMPb$X[which(regulon_GMPb$pval_wilcox < 0.15)])
sig_regulons_gmp_final <- unique(sig_regulons_gmp)

regulon_MLPa <- read.csv(file=paste0(REGULON_DIR,"MLP_a_wilcox_data_no_outliers.csv"))
regulon_MLPb <- read.csv(file=paste0(REGULON_DIR,"MLP_b_wilcox_data_no_outliers.csv"))
sig_regulons_mlp <- c(regulon_MLPa$X[which(regulon_MLPa$pval_wilcox < 0.15)], regulon_MLPb$X[which(regulon_MLPb$pval_wilcox < 0.15)])
sig_regulons_mlp_final <- unique(sig_regulons_mlp)

regulon_PreBNK <- read.csv(file=paste0(REGULON_DIR,"PreBNK_wilcox_data_no_outliers.csv"))
sig_regulons_PreBNK_final <- regulon_PreBNK$X[which(regulon_PreBNK$pval_wilcox < 0.15)]

regulon_MEPa <- read.csv(file=paste0(REGULON_DIR,"MEP_a_wilcox_data_no_outliers.csv"))
regulon_MEPb <- read.csv(file=paste0(REGULON_DIR,"MEP_b_wilcox_data_no_outliers.csv"))
regulon_MEPc <- read.csv(file=paste0(REGULON_DIR,"MEP_c_wilcox_data_no_outliers.csv"))
sig_regulons_mep <- c(regulon_MEPa$X[which(regulon_MEPa$pval_wilcox < 0.15)], regulon_MEPb$X[which(regulon_MEPb$pval_wilcox < 0.15)], regulon_MEPc$X[which(regulon_MEPc$pval_wilcox < 0.15)])
sig_regulons_mep_final <- unique(sig_regulons_mep)

sig_regulons_list <- list(sig_regulons_cmp_final, sig_regulons_gmp_final, sig_regulons_mep_final, sig_regulons_mlp_final, sig_regulons_PreBNK_final)

plot_list <- list()
for(i in 1:length(clusters)){
  #Read in statistical enrichment data
  enrichment_data <- read.table(file=paste0("emmreml_downstream/HINT_motif_enrichment_analysis/hint_Enrichment/",clusters[i],"_homer_DE_peaks/fulltest_statistics.txt"), header=TRUE, fill=TRUE)
  
  #Read in bedtools overlap data
  bedtools_overlap <- read.table(file=paste0("emmreml_downstream/HINT_motif_enrichment_analysis/bedtools_output/peaks_with_motifs/", clusters[i],"_HINT_out.bed"))
  
  #Get the percent DA peaks bound by each enriched TF
  tfs <- enrichment_data$FACTOR
  percentage <- vector(length=length(tfs))
  for(j in 1:length(tfs)){
    num_DA_peaks <- unique(bedtools_overlap$V4[which(bedtools_overlap$V8 %in% tfs[j])])
    total_peaks <- unique(bedtools_overlap$V4)
    percentage[j] <- length(num_DA_peaks)/length(total_peaks)
  }
  
  #Read in HSC regulon data
  REGULON_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic/get_diff_regulons_outs/"
  regulon_HSCa <- read.csv(file=paste0(REGULON_DIR,"HSC_a_wilcox_data_no_outliers.csv"))
  regulon_HSCb <- read.csv(file=paste0(REGULON_DIR,"HSC_b_wilcox_data_no_outliers.csv"))
  sig_regulons <- c(regulon_HSCa$X[which(round(regulon_HSCa$pval_wilcox,2) <= 0.1)], regulon_HSCb$X[which(round(regulon_HSCb$pval_wilcox,2) <= 0.1)])
  sig_regulons_final <- unique(sig_regulons)
  
  #HSC regulon data
  regulon_lab <- rep("no", length(tfs))
  TF_info_copy <- TF_info
  str_sub(TF_info_copy$TF.motif, start = 1, end = 9, omit_na = FALSE) <- ""
  for(j in 1:length(sig_regulons_final)){
    if(sig_regulons_final[j] %in% TF_info_copy$TF.motif)
    {
      TF_matched_group <- unique(TF_info_copy$clusters[which(TF_info_copy$TF.motif == sig_regulons_final[j] )]) #cluster number
      group_motifs <- TF_info$TF.motif[which(TF_info$clusters %in% TF_matched_group)] #motif 
      regulon_lab[which(tfs %in% group_motifs)] <- "yes" # regulon lab = + HSC regulon activity 
    }
  }
  
  
  #matched cluster regulon data
  regulon_lab_clust <- rep("no", length(tfs))
  for(j in 1:length(sig_regulons_list[[i]])){
    if(sig_regulons_list[[i]][j] %in% TF_info_copy$TF.motif)
    {
      TF_matched_group <- unique(TF_info_copy$clusters[which(TF_info_copy$TF.motif == sig_regulons_list[[i]][j])])
      group_motifs <- TF_info$TF.motif[which(TF_info$clusters %in% TF_matched_group)]
      regulon_lab_clust[which(tfs %in% group_motifs)] <- "yes"
    }
  }
  
  #read in and format TF EXPRESSION data
  TF_expression_dir <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/MASH_emmreml_downstream/DE_overlapping_TFs/"
  FULL_TF_data <- read.csv(file=paste0(TF_expression_dir, "TF_exp_FULL_data.csv"))
  HSCa_exp_data <- subset(FULL_TF_data, FULL_TF_data$celltype == "HSC_a")
  HSCb_exp_data <- subset(FULL_TF_data, FULL_TF_data$celltype == "HSC_b")
  HSC_sig <- unique(c(as.character(HSCa_exp_data$TF_name[which(HSCa_exp_data$val == 1)]), as.character(HSCb_exp_data$TF_name[which(HSCb_exp_data$val == 1)])))
  
  hsc_expression_lab <- rep("no", length(tfs))
  for(j in 1:length(HSC_sig)){
    if(HSC_sig[j] %in% TF_info_copy$TF.motif)
    {
      TF_matched_group <- unique(TF_info_copy$clusters[which(TF_info_copy$TF.motif == HSC_sig[j])])
      group_motifs <- TF_info$TF.motif[which(TF_info$clusters %in% TF_matched_group)]
      hsc_expression_lab[which(tfs %in% group_motifs)] <- "yes"
    }
  }
  
  
  ##combine HSC_regulon and HSC_exp into a single ACTIVATION label 
  HSC_combined_activity <- paste0(regulon_lab, "_", hsc_expression_lab)
  HSC_combined_activity[which(grepl("yes",HSC_combined_activity)==TRUE)] <- "HSC active"
  HSC_combined_activity[which(regulon_lab_clust=="yes")] <- "current_active"
  
  
  lab_fdr <- rep("no", length(tfs))
  lab_fdr[which(enrichment_data$CORR.P.VALUE < 0.05)] <- "yes"
  
  regulon_lab_edited <- regulon_lab
  regulon_lab_edited[which(lab_fdr=="no")] <- "no"
  
  regulon_lab_clust_edited <- regulon_lab_clust
  regulon_lab_clust_edited[which(lab_fdr=="no")] <- "no"
  
  hsc_expression_lab_edited <- hsc_expression_lab
  hsc_expression_lab_edited[which(lab_fdr=="no")] <- "no"
  
  
  df_scatter <- data.frame(
    name <- tfs,
    peak_fdr <- -log10(enrichment_data$CORR.P.VALUE),
    perc <- percentage,
    hsc_reg <- regulon_lab_edited,
    clust_reg <- regulon_lab_clust_edited,
    hsc_exp <- hsc_expression_lab_edited,
    l <- HSC_combined_activity
  )
  colnames(df_scatter) <- c("name", "peak_fdr", "perc", "hsc_reg", "clust_reg","hsc_exp", "l")
  write.csv(df_scatter, file=paste0(OUT_DIR,"scatterplot_data_individual/0.1_scatterplot_data_",clusters[i],".csv"))
  
  
   

  # p <- ggplot(df_scatter, aes(x=perc, y=peak_fdr)) + geom_point(aes(size=peak_fdr, fill=l), shape=21, color="black") + theme_test() +
  #   labs(y="-log10(FDR)", x="Percent sites bound", size="-log10(FDR)") +
  #   scale_fill_manual(values=c("#0072B5","#EE4C97","grey")) + geom_vline(xintercept = 0.2, color="#EF5350", linetype="dashed") + geom_hline(yintercept=15, color="#EF5350", linetype="dashed")
  # plot_list[[i]] <- p

}

# for(k in 1:length(plot_list)) {
#   file_name = paste0(OUT_DIR, paste0(clusters[k],"_scatter_SIMPLE_individual.tiff"))
#   tiff(file_name, units="in", width=4, height=3, res=400)
#   print(plot_list[[k]])
#   dev.off()
#   print(k)
#   
# }