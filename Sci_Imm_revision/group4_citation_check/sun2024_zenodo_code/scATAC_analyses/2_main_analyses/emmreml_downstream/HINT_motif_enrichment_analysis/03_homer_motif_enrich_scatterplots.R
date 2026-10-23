#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

.libPaths("/project/lbarreiro/USERS/sarah/Rlibs_new")

library(stringr)
library(textTinyR)
library(pbapply)
library(readr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)
library(ggbeeswarm)

setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/")
OUT_DIR <- "emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/"

#Read in data for scatter plot
plot_data <- read.csv(file=paste0(OUT_DIR, "FDR_percentage_FULL_table.csv")) # from hint_motif_enrich_results_summary.R
colnames(plot_data) <- c("row_num", "TF_name", "celltype", "neglogfdr", "percent", "TF_num")


# from hint_motif_enrich_results_summary.R
TF_info <- read.csv(file=paste0("emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/TF_cluster_information_with_names.csv"))

#Read in regulon data from each cluster
#Read in all regulon data
REGULON_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic/get_diff_regulons_outs/"
regulon_CMPa <- read.csv(file=paste0(REGULON_DIR,"CMP_a_wilcox_data_no_outliers.csv"))
regulon_CMPb <- read.csv(file=paste0(REGULON_DIR,"CMP_b_wilcox_data_no_outliers.csv"))
regulon_CMPc <- read.csv(file=paste0(REGULON_DIR,"CMP_c_wilcox_data_no_outliers.csv"))
sig_regulons_cmp <- c(as.character(regulon_CMPa$X[which(regulon_CMPa$pval_wilcox < 0.15)]), as.character(regulon_CMPb$X[which(regulon_CMPb$pval_wilcox < 0.15)]), as.character(regulon_CMPc$X[which(regulon_CMPc$pval_wilcox < 0.15)]))
sig_regulons_cmp_final <- unique(sig_regulons_cmp)

regulon_GMPa <- read.csv(file=paste0(REGULON_DIR,"GMP_a_wilcox_data_no_outliers.csv"))
regulon_GMPb <- read.csv(file=paste0(REGULON_DIR,"GMP_b_wilcox_data_no_outliers.csv"))
sig_regulons_gmp <- c(as.character(regulon_GMPa$X[which(regulon_GMPa$pval_wilcox < 0.15)]), as.character(regulon_GMPb$X[which(regulon_GMPb$pval_wilcox < 0.15)]))
sig_regulons_gmp_final <- unique(sig_regulons_gmp)

regulon_MLPa <- read.csv(file=paste0(REGULON_DIR,"MLP_a_wilcox_data_no_outliers.csv"))
regulon_MLPb <- read.csv(file=paste0(REGULON_DIR,"MLP_b_wilcox_data_no_outliers.csv"))
sig_regulons_mlp <- c(as.character(regulon_MLPa$X[which(regulon_MLPa$pval_wilcox < 0.15)]), as.character(regulon_MLPb$X[which(regulon_MLPb$pval_wilcox < 0.15)]))
sig_regulons_mlp_final <- unique(sig_regulons_mlp)

regulon_PreBNK <- read.csv(file=paste0(REGULON_DIR,"PreBNK_wilcox_data_no_outliers.csv"))
sig_regulons_PreBNK_final <- as.character(regulon_PreBNK$X[which(regulon_PreBNK$pval_wilcox < 0.15)])

regulon_MEPa <- read.csv(file=paste0(REGULON_DIR,"MEP_a_wilcox_data_no_outliers.csv"))
regulon_MEPb <- read.csv(file=paste0(REGULON_DIR,"MEP_b_wilcox_data_no_outliers.csv"))
regulon_MEPc <- read.csv(file=paste0(REGULON_DIR,"MEP_c_wilcox_data_no_outliers.csv"))
sig_regulons_mep <- c(as.character(regulon_MEPa$X[which(regulon_MEPa$pval_wilcox < 0.15)]), as.character(regulon_MEPb$X[which(regulon_MEPb$pval_wilcox < 0.15)]), as.character(regulon_MEPc$X[which(regulon_MEPc$pval_wilcox < 0.15)]))
sig_regulons_mep_final <- unique(sig_regulons_mep)

regulon_HSCa <- read.csv(file=paste0(REGULON_DIR,"HSC_a_wilcox_data_no_outliers.csv"))
regulon_HSCb <- read.csv(file=paste0(REGULON_DIR,"HSC_b_wilcox_data_no_outliers.csv"))
sig_regulons_hsc <- c(as.character(regulon_HSCa$X[which(regulon_HSCa$pval_wilcox < 0.15)]), as.character(regulon_HSCb$X[which(regulon_HSCb$pval_wilcox < 0.15)]))
sig_regulons_hsc_final <- unique(sig_regulons_hsc)

sig_regulons_list <- list(sig_regulons_hsc_final, sig_regulons_cmp_final, sig_regulons_gmp_final, sig_regulons_mep_final, sig_regulons_mlp_final, sig_regulons_PreBNK_final)


#read in and format TF expression data
TF_expression_dir <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis_Raul/MASH_emmreml_downstream/DE_overlapping_TFs/"
FULL_TF_data <- read.csv(file=paste0(TF_expression_dir, "TF_exp_FULL_data.csv"))
HSCa_exp_data <- subset(FULL_TF_data, FULL_TF_data$celltype == "HSC_a")
HSCb_exp_data <- subset(FULL_TF_data, FULL_TF_data$celltype == "HSC_b")
HSC_sig <- unique(c(as.character(HSCa_exp_data$TF_name[which(HSCa_exp_data$val == 1)]), as.character(HSCb_exp_data$TF_name[which(HSCb_exp_data$val == 1)])))


# Scatterplots simplified version -----------------------------------------


plot_list <- list()
plot_list2 <- list()
plot_list3 <- list()
plot_list4 <- list()
clusters <- c("HSC", "CMP", "GMP", "MEP", "MLP", "PreBNK")
for(i in 1:length(clusters)){
  
  plot_data_final <- subset(plot_data, plot_data$celltype == clusters[i])
  
  if (clusters[i] == "GMP")
  {
    plot_lab <- plot_data_final$TF_name
    plot_lab[which(plot_data_final$percent < 0.15)] <- ""
    plot_lab[which(plot_data_final$neglogfdr < 12)] <- ""
    plot_data_final$plot_lab <- plot_lab
  }
  else{
    plot_lab <- plot_data_final$TF_name
    plot_lab[which(plot_data_final$percent < 0.15)] <- ""
    plot_lab[which(plot_data_final$neglogfdr < 20)] <- ""
    plot_data_final$plot_lab <- plot_lab
  }
  
  
  #Make another lab to indicate TFs with sig regulons in MATCHED CLUSTER
  regulon_group <- sig_regulons_list[[i]]
  matched_regulon <- rep("no", length(row.names(plot_data_final)))
  TF_info_copy <- TF_info
  str_sub(TF_info_copy$TF.motif, start = 1, end = 9, omit_na = FALSE) <- ""
  for(j in 1:length(regulon_group)){
    if(regulon_group[j] %in% TF_info_copy$TF.motif)
    {
      TF_matched <- unique(TF_info_copy$clusters[which(TF_info_copy$TF.motif == regulon_group[j] )])
      matched_regulon[which(plot_data_final$TF_num %in% TF_matched)] <- "yes"
    }
  }
  matched_regulon[which(matched_regulon == "")] <- "no"
  
  
  #Make another lab to indicate TFs with sig regulons in HSC
  regulon_group <- sig_regulons_hsc_final
  HSC_regulon <- rep("no", length(row.names(plot_data_final)))
  TF_info_copy <- TF_info
  str_sub(TF_info_copy$TF.motif, start = 1, end = 9, omit_na = FALSE) <- ""
  for(j in 1:length(regulon_group)){
    if(regulon_group[j] %in% TF_info_copy$TF.motif)
    {
      TF_matched <- unique(TF_info_copy$clusters[which(TF_info_copy$TF.motif == regulon_group[j] )])
      HSC_regulon[which(plot_data_final$TF_num %in% TF_matched)] <- "yes"
    }
  }
  HSC_regulon[which(HSC_regulon == "")] <- "no"
  
  #Make another lab to indicate TFs with diff expression in HSC
  DE_exp_group <- HSC_sig
  HSC_exp <- rep("no", length(row.names(plot_data_final)))
  TF_info_copy <- TF_info
  str_sub(TF_info_copy$TF.motif, start = 1, end = 9, omit_na = FALSE) <- ""
  for(j in 1:length(DE_exp_group)){
    if(DE_exp_group[j] %in% TF_info_copy$TF.motif)
    {
      TF_matched <- unique(TF_info_copy$clusters[which(TF_info_copy$TF.motif == DE_exp_group[j] )])
      HSC_exp[which(plot_data_final$TF_num %in% TF_matched)] <- "yes"
    }
  }
  HSC_exp[which(HSC_exp == "")] <- "no"
  
  
  #Edit plot data final to include regulon labels and TF DE exp labels 
  
  ##combine HSC_regulon and HSC_exp into a single ACTIVATION label 
  table(paste0(HSC_regulon, "_", HSC_exp))
  HSC_combined_activity <- paste0(HSC_regulon, "_", HSC_exp)
  HSC_combined_activity[which(grepl("yes",HSC_combined_activity)==TRUE)] <- "HSC active"
  HSC_combined_activity[which(matched_regulon=="yes")] <- "current_active"
  plot_data_final$regulon_combined <- HSC_combined_activity
  #plot_data_final$regulon_combined <- factor(plot_data_final$regulon_combined, c("no_no", "HSC_active", "current_active"))
  #write.csv(plot_data_final, file=paste0(OUT_DIR,"scatterplot_data/scatterplot_data_",clusters[i],".csv"))
  
  plot_data_final_violin <- subset(plot_data_final, plot_data_final$regulon_combined %in% c("no_no", "HSC active"))
  plot_data_final_violin$combined_score <- plot_data_final_violin$neglogfdr*plot_data_final_violin$percent
  
  
  #PLOT
  
  if(i %in% c(2)){
    p <- ggplot(plot_data_final, aes(x=percent, y=neglogfdr)) + geom_point(aes(size=neglogfdr, fill=regulon_combined), shape=21, color="black") + theme_test() +
      geom_text_repel(aes(label=plot_lab), size=3) + labs(y="-log10(FDR)", x="Percent sites bound", size="-log10(FDR)") + ylim(c(0,80)) + xlim(c(0,0.5)) +
      scale_fill_manual(values=c("#0072B5","#EE4C97","grey")) + geom_vline(xintercept = 0.15, color="#EF5350", linetype="dashed") + geom_hline(yintercept=20, color="#EF5350", linetype="dashed")
    plot_list[[i]] <- p
    
    ##Violin plot
    wilcox_clust_fdr <- wilcox.test(plot_data_final_violin$neglogfdr[which(plot_data_final_violin$regulon_combined == "HSC active")], plot_data_final_violin$neglogfdr[which(plot_data_final_violin$regulon_combined == "no_no")])$p.value
    
    p2 <- ggplot(plot_data_final_violin, aes(x=regulon_combined, y=neglogfdr)) + geom_violin(aes(fill=regulon_combined))+
      theme_classic() + labs(x="", y="-log10(FDR)", title=paste0(clusters[i], ": P = ",round(wilcox_clust_fdr,3))) +  theme(legend.position="none") +
      scale_fill_manual(values=c("#EE4C97","grey" )) + theme(axis.text=element_text(size=12)) + coord_cartesian(ylim = c(0,40))
    plot_list2[[i]] <- p2
    
    wilcox_clust_percent <- wilcox.test(plot_data_final_violin$percent[which(plot_data_final_violin$regulon_combined == "HSC active")], plot_data_final_violin$percent[which(plot_data_final_violin$regulon_combined == "no_no")])$p.value
    
    p3 <- ggplot(plot_data_final_violin, aes(x=regulon_combined, y=percent)) + geom_violin(aes(fill=regulon_combined))+
      theme_classic() + labs(x="", y="Percent DR peaks with Motif", title=paste0(clusters[i], ": P = ",round(wilcox_clust_percent,3))) +  theme(legend.position="none") +
      scale_fill_manual(values=c("#EE4C97","grey" )) + theme(axis.text=element_text(size=12)) + coord_cartesian(ylim = c(0,0.8))
    plot_list3[[i]] <- p3
    
    wilcox_clust_combined <- wilcox.test(plot_data_final_violin$combined_score[which(plot_data_final_violin$regulon_combined == "HSC active")], plot_data_final_violin$combined_score[which(plot_data_final_violin$regulon_combined == "no_no")])$p.value
    
    p4 <- ggplot(plot_data_final_violin, aes(x=regulon_combined, y=combined_score)) + geom_violin(aes(fill=regulon_combined))+
      theme_classic() + labs(x="", y="-log10(FDR)*%sites", title=paste0(clusters[i], ": P = ",round(wilcox_clust_combined,3))) +  theme(legend.position="none") +
      scale_fill_manual(values=c("#EE4C97","grey" )) + theme(axis.text=element_text(size=12)) + coord_cartesian(ylim = c(0,15))
    plot_list4[[i]] <- p4
    
  }
  else if(i %in% c(3)) {
    p <- ggplot(plot_data_final, aes(x=percent, y=neglogfdr)) + geom_point(aes(size=neglogfdr, fill=regulon_combined), shape=21, color="black") + theme_test() +
      geom_text_repel(aes(label=plot_lab), size=3) + labs(y="-log10(FDR)", x="Percent sites bound", size="-log10(FDR)") + ylim(c(0,45)) + xlim(c(0,0.5)) +
      scale_fill_manual(values=c("#0072B5","#EE4C97","grey")) + geom_vline(xintercept = 0.15, color="#EF5350", linetype="dashed") + geom_hline(yintercept=12, color="#EF5350", linetype="dashed")
    plot_list[[i]] <- p
    
    ##Violin plot
    wilcox_clust_fdr <- wilcox.test(plot_data_final_violin$neglogfdr[which(plot_data_final_violin$regulon_combined == "HSC active")], plot_data_final_violin$neglogfdr[which(plot_data_final_violin$regulon_combined == "no_no")])$p.value
    
    p2 <- ggplot(plot_data_final_violin, aes(x=regulon_combined, y=neglogfdr)) + geom_violin(aes(fill=regulon_combined))+
      theme_classic() + labs(x="", y="-log10(FDR)", title=paste0(clusters[i], ": P = ",round(wilcox_clust_fdr,3))) +  theme(legend.position="none") +
      scale_fill_manual(values=c("#EE4C97","grey" )) + theme(axis.text=element_text(size=12)) + coord_cartesian(ylim = c(0,40))
    plot_list2[[i]] <- p2
    
    wilcox_clust_percent <- wilcox.test(plot_data_final_violin$percent[which(plot_data_final_violin$regulon_combined == "HSC active")], plot_data_final_violin$percent[which(plot_data_final_violin$regulon_combined == "no_no")])$p.value
    
    p3 <- ggplot(plot_data_final_violin, aes(x=regulon_combined, y=percent)) + geom_violin(aes(fill=regulon_combined))+
      theme_classic() + labs(x="", y="Percent DR peaks with Motif", title=paste0(clusters[i], ": P = ",round(wilcox_clust_percent,3))) +  theme(legend.position="none") +
      scale_fill_manual(values=c("#EE4C97","grey" )) + theme(axis.text=element_text(size=12)) + coord_cartesian(ylim = c(0,0.8))
    plot_list3[[i]] <- p3
    
    wilcox_clust_combined <- wilcox.test(plot_data_final_violin$combined_score[which(plot_data_final_violin$regulon_combined == "HSC active")], plot_data_final_violin$combined_score[which(plot_data_final_violin$regulon_combined == "no_no")])$p.value
    
    p4 <- ggplot(plot_data_final_violin, aes(x=regulon_combined, y=combined_score)) + geom_violin(aes(fill=regulon_combined))+
      theme_classic() + labs(x="", y="-log10(FDR)*%sites", title=paste0(clusters[i], ": P = ",round(wilcox_clust_combined,3))) +  theme(legend.position="none") +
      scale_fill_manual(values=c("#EE4C97","grey" )) + theme(axis.text=element_text(size=12)) + coord_cartesian(ylim = c(0,15))
    plot_list4[[i]] <- p4
  }
  else {
    p <- ggplot(plot_data_final, aes(x=percent, y=neglogfdr)) + geom_point(aes(size=neglogfdr, fill=regulon_combined), shape=21, color="black") + theme_test() +
      geom_text_repel(aes(label=plot_lab), size=3) + labs(y="-log10(FDR)", x="Percent sites bound", size="-log10(FDR)") + 
      scale_fill_manual(values=c("#0072B5","#EE4C97","grey"))
    plot_list[[i]] <- p
    
    ##Violin plot
    wilcox_clust_fdr <- wilcox.test(plot_data_final_violin$neglogfdr[which(plot_data_final_violin$regulon_combined == "HSC active")], plot_data_final_violin$neglogfdr[which(plot_data_final_violin$regulon_combined == "no_no")])$p.value
    
    p2 <- ggplot(plot_data_final_violin, aes(x=regulon_combined, y=neglogfdr)) + geom_violin(aes(fill=regulon_combined))+
      theme_classic() + labs(x="", y="-log10(FDR)", title=paste0(clusters[i], ": P = ",round(wilcox_clust_fdr,3))) +  theme(legend.position="none") +
      scale_fill_manual(values=c("#EE4C97","grey" )) + theme(axis.text=element_text(size=12)) + coord_cartesian(ylim = c(0,40))
    plot_list2[[i]] <- p2
    
    wilcox_clust_percent <- wilcox.test(plot_data_final_violin$percent[which(plot_data_final_violin$regulon_combined == "HSC active")], plot_data_final_violin$percent[which(plot_data_final_violin$regulon_combined == "no_no")])$p.value
    
    p3 <- ggplot(plot_data_final_violin, aes(x=regulon_combined, y=percent)) + geom_violin(aes(fill=regulon_combined))+
      theme_classic() + labs(x="", y="Percent DR peaks with Motif", title=paste0(clusters[i], ": P = ",round(wilcox_clust_percent,3))) +  theme(legend.position="none") +
      scale_fill_manual(values=c("#EE4C97","grey" )) + theme(axis.text=element_text(size=12)) + coord_cartesian(ylim = c(0,0.8))
    plot_list3[[i]] <- p3
    
    wilcox_clust_combined <- wilcox.test(plot_data_final_violin$combined_score[which(plot_data_final_violin$regulon_combined == "HSC active")], plot_data_final_violin$combined_score[which(plot_data_final_violin$regulon_combined == "no_no")])$p.value
    
    p4 <- ggplot(plot_data_final_violin, aes(x=regulon_combined, y=combined_score)) + geom_violin(aes(fill=regulon_combined))+
      theme_classic() + labs(x="", y="-log10(FDR)*%sites", title=paste0(clusters[i], ": P = ",round(wilcox_clust_combined,3))) +  theme(legend.position="none") +
      scale_fill_manual(values=c("#EE4C97","grey" )) + theme(axis.text=element_text(size=12)) + coord_cartesian(ylim = c(0,15))
    plot_list4[[i]] <- p4
  }
  
  
  #write.csv(plot_data_final, file=paste0(OUT_DIR, clusters[i], "_final_scatterplot_data_SIMPLE.csv"))
  
}


for(k in 1:length(plot_list)) {
  file_name = paste0(OUT_DIR, paste0(clusters[k],"_scatter_SIMPLE.tiff"))
  tiff(file_name, units="in", width=4, height=3, res=400)
  print(plot_list[[k]])
  dev.off()
  print(k)
  
}

for(k in 1:length(plot_list2)) {
  file_name = paste0(OUT_DIR, paste0(clusters[k],"_violin_FDR.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plot_list2[[k]])
  dev.off()
  print(k)
  
}

for(k in 1:length(plot_list3)) {
  file_name = paste0(OUT_DIR, paste0(clusters[k],"_violin_percent.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plot_list3[[k]])
  dev.off()
  print(k)
  
}

for(k in 1:length(plot_list4)) {
  file_name = paste0(OUT_DIR, paste0(clusters[k],"_violin_combined.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plot_list4[[k]])
  dev.off()
  print(k)
  
}