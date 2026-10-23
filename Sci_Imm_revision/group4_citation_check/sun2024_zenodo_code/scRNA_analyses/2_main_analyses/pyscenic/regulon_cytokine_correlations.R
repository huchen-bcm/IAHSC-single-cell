library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(reticulate)
library(tidyr)
library(reshape2)
library(plyr)
library(grid)
library(gridExtra)
library(data.table)
library(SCENIC)
library(SCopeLoomR)
library(doRNG)
library(ggrepel)

setwd("/scRNA_analyses/2_main_analyses/pyscenic")
OUT_DIR <- "regulon_cytokine_correlations/"
dir.create(OUT_DIR)

## Get data from sce object:
seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"




#Analysis excluding S2, S3, S12 (filtered out due DR gene analysis due to fewer than 20 cells)
auc_mat <- read.csv(file=paste0("scenic_cluster_HSc_a_auc_mtx.csv"))
df <- seurat_obj@meta.data
samples <- unique(df$UNIQUE_ID)
for(i in 1:length(samples)){
  temp <- row.names(subset(df, df$UNIQUE_ID %in% samples[i]))
  auc_subset <- subset(auc_mat, auc_mat$Cell %in% temp)
  row.names(auc_subset) <- auc_subset$Cell
  auc_subset <- auc_subset[,2:length(colnames(auc_subset))]
  means <- colMeans(auc_subset)
  if(i ==1)
  {
    donor_data <- means
  }
  else
  {
    donor_data <- cbind(donor_data, means)
  }
}
colnames(donor_data) <- samples

meta_data <- read.csv(file="../all_samples_meta_data.csv")



Tm3_all <- c("SNG-LB-SS-1S-RS-S1-CD34neg_S1_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S5-CD34neg_S4_R1_001_Tm3_BCG",
             "SNG-LB-SS-1S-RS-S6-CD34neg_S5_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S8-CD34neg_S7_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S9-CD34neg_S8_R1_001_Tm3_BCG",
             "SNG-LB-SS-1S-RS-S10-CD34neg_S9_R1_001_Tm3_BCG",  "SNG-LB-SS-1S-RS-S14-CD34neg_S13_R1_001_Tm3_BCG",
             "SNG-LB-SS-1S-RS-S15-CD34neg_S14_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S16-CD34neg_S15_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S17-CD34neg_S16_R1_001_Tm3_BCG",
             "SNG-LB-SS-1S-RS-S19-CD34neg_S18_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S20-CD34neg_S19_R1_001_Tm3_BCG", "SNG-LB-SS-1S-RS-S21-CD34neg_S20_R1_001_Tm3_BCG", 
             "SNG-LB-SS-1S-RS-S7-CD34neg_S6_R1_001_Tm3_CTL", "SNG-LB-SS-1S-RS-S11-CD34neg_S10_R1_001_Tm3_CTL",
             "SNG-LB-SS-1S-RS-S13-CD34neg_S12_R1_001_Tm3_CTL", "SNG-LB-SS-1S-RS-S18-CD34neg_S17_R1_001_Tm3_CTL")


Td0_all <- c("SNG-LB-SS-1S-RS-S1-CD34neg_S1_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S5-CD34neg_S4_R1_001_Td0_BCG",
             "SNG-LB-SS-1S-RS-S6-CD34neg_S5_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S8-CD34neg_S7_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S9-CD34neg_S8_R1_001_Td0_BCG",
             "SNG-LB-SS-1S-RS-S10-CD34neg_S9_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S14-CD34neg_S13_R1_001_Td0_BCG",
             "SNG-LB-SS-1S-RS-S15-CD34neg_S14_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S16-CD34neg_S15_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S17-CD34neg_S16_R1_001_Td0_BCG",
             "SNG-LB-SS-1S-RS-S19-CD34neg_S18_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S20-CD34neg_S19_R1_001_Td0_BCG", "SNG-LB-SS-1S-RS-S21-CD34neg_S20_R1_001_Td0_BCG", 
             "SNG-LB-SS-1S-RS-S7-CD34neg_S6_R1_001_Td0_CTL", "SNG-LB-SS-1S-RS-S11-CD34neg_S10_R1_001_Td0_CTL",
             "SNG-LB-SS-1S-RS-S13-CD34neg_S12_R1_001_Td0_CTL", "SNG-LB-SS-1S-RS-S18-CD34neg_S17_R1_001_Td0_CTL")


logfc_mat <- data.frame()
count <- 0
name_short <- vector()
name_count <- 1
for(i in 1:length(Tm3_all)){
  name <- Tm3_all[i]
  if(name %in% colnames(donor_data))
  {
    d <- name
    Td0_pair_name <- Td0_all[which(Tm3_all == d)]
    if(Td0_pair_name %in% colnames(donor_data) && is.na(donor_data[,Td0_pair_name])==FALSE)
    {
      Td0_col <- donor_data[,Td0_pair_name]
      Tm3_col <- donor_data[,d]
      log2fc <- Tm3_col - Td0_col
      name_short[name_count] <- subset(meta_data, meta_data$Sample %in% as.character(d))$donor
      name_count <- name_count+1
      if(count==0)
      {
        logfc_mat <- log2fc
      }
      else
      {
        logfc_mat <- cbind(logfc_mat, log2fc)
      }
      
      count <- count + 1
    }
    
  }
  
}

colnames(logfc_mat) <- name_short
row.names(logfc_mat) <- substr(row.names(logfc_mat),1,nchar(row.names(logfc_mat))-3)


wilcox_data <- read.csv(file=paste0("get_diff_regulons_outs/HSC_a_wilcox_data_no_outliers.csv"), header=TRUE)
sig_regulons <- wilcox_data$X[which(wilcox_data$pval_wilcox < 0.15)]
logfc_mat_sub <- subset(logfc_mat, row.names(logfc_mat) %in% sig_regulons)
logfc_mat_t <- t(logfc_mat_sub)
write.table(logfc_mat_t, file=paste0(OUT_DIR, "HSCa_sig_regulon_logfc.txt"))
logfc_mat_final <- scale(logfc_mat_t)

#IL1B
regulon_name <- vector(length=length(colnames(logfc_mat_final)))
pval_il1b <- vector(length=length(colnames(logfc_mat_final)))
rho_il1b <- vector(length=length(colnames(logfc_mat_final)))
plot_list1 <- list()
for(i in 1:length(colnames(logfc_mat_final))){
  
    regulon_name[i] <- colnames(logfc_mat_final)[i]
    reg_vals <- logfc_mat_final[,i]
    
    cytokine_data <- data.frame(fread("../PBMC_cytokine_data_CHM.txt"))
    meta_data <- read.csv(file="../all_samples_meta_data.csv")
    il1b <- vector(length=length(reg_vals))
    treatment <- vector(length=length(reg_vals))
    
    for(j in 1:length(reg_vals))
    {
      cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == names(reg_vals)[j])
      cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
      il1b[j] <- as.numeric(cyt_subset2$IL1B_FC)
      treatment[j] <- unique(meta_data$vaccination[which(meta_data$donor == names(reg_vals)[j])])
      print(j)
    }
    
    val <- cor.test(reg_vals, il1b, method="spearman")
    pval_il1b[i] <- val$p.value
    rho_il1b[i] <- val$estimate
    
    df_reg <- data.frame(
      score <- reg_vals,
      c <- il1b,
      t <- treatment
    )
    colnames(df_reg) <- c("score", "c", "Group")
    
    ann_text<-data.frame(
      x = -0.5, y = 2,
      label = paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_il1b[i],4))
    )
    
    p1 <- ggplot(df_reg, aes(x=score, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_linedraw() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
      geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + 
      labs(x = paste0(colnames(logfc_mat_final)[i]," score"), y ="FC IL1B by PBMCs", title= paste0(colnames(logfc_mat_final)[i]," score vs. IL1B")) +
      scale_fill_manual(values=c("#42B540", "grey"))
    plot_list1[[i]] <- p1

    
    # ## i == 25
    # if(colnames(logfc_mat_final)[i] == "KLF5"){
    #   file_name = paste0(OUT_DIR,"KLF5_hsca_corr_IL1B.tiff")
    #   tiff(file_name, units="in", width=4, height=3, res=300)
    #   ggplot(df_reg, aes(x=score, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_classic() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
    #     geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab(paste0(colnames(logfc_mat_final)[i]," score")) +
    #     ylab("FC IL1B by PBMCs") + ggtitle(paste0(colnames(logfc_mat_final)[i]," score vs. IL1B"))+
    #     geom_text(data=ann_text, aes(x=-0.4, y=2, label=label)) + scale_fill_manual(values=c("#F1788D", "grey")) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
    #   dev.off()
    # }
    # 
    # ## i == 26
    # if(colnames(logfc_mat_final)[i] == "KLF6"){
    #   file_name = paste0(OUT_DIR,"KLF6_hsca_corr_IL1B.tiff")
    #   tiff(file_name, units="in", width=4, height=3, res=300)
    #   ggplot(df_reg, aes(x=score, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_classic() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
    #     geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab(paste0(colnames(logfc_mat_final)[i]," score")) +
    #     ylab("FC IL1B by PBMCs") + ggtitle(paste0(colnames(logfc_mat_final)[i]," score vs. IL1B"))+
    #     geom_text(data=ann_text, aes(x=-0.4, y=2, label=label)) + scale_fill_manual(values=c("#F1788D", "grey")) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
    #   dev.off()
    # }
    # 
    # ## i == 9
    # if(colnames(logfc_mat_final)[i] == "EGR1"){
    #   file_name = paste0(OUT_DIR,"EGR1_hsca_corr_IL1B.tiff")
    #   tiff(file_name, units="in", width=4, height=3, res=300)
    #   ggplot(df_reg, aes(x=score, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_classic() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
    #     geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab(paste0(colnames(logfc_mat_final)[i]," score")) +
    #     ylab("FC IL1B by PBMCs") + ggtitle(paste0(colnames(logfc_mat_final)[i]," score vs. IL1B"))+
    #     geom_text(data=ann_text, aes(x=-0.4, y=2, label=label)) + scale_fill_manual(values=c("#F1788D", "grey")) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
    #   dev.off()
    # }
    # 
    # ## i == 5
    # if(colnames(logfc_mat_final)[i] == "CEBPB"){
    #   file_name = paste0(OUT_DIR,"CEBPB_hsca_corr_IL1B.tiff")
    #   tiff(file_name, units="in", width=4, height=3, res=300)
    #   ggplot(df_reg, aes(x=score, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_classic() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
    #     geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab(paste0(colnames(logfc_mat_final)[i]," score")) +
    #     ylab("FC IL1B by PBMCs") + ggtitle(paste0(colnames(logfc_mat_final)[i]," score vs. IL1B"))+
    #     geom_text(data=ann_text, aes(x=-0.4, y=2, label=label)) + scale_fill_manual(values=c("#F1788D", "grey")) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
    #   dev.off()
    # }
    # 
}




#IL6
regulon_name <- vector(length=length(colnames(logfc_mat_final)))
pval_il6 <- vector(length=length(colnames(logfc_mat_final)))
rho_il6 <- vector(length=length(colnames(logfc_mat_final)))
plot_list2 <- list()
for(i in 1:length(colnames(logfc_mat_final))){
  
  regulon_name[i] <- colnames(logfc_mat_final)[i]
  reg_vals <- logfc_mat_final[,i]
  
  cytokine_data <- data.frame(fread("../PBMC_cytokine_data_CHM.txt"))
  meta_data <- read.csv(file="../all_samples_meta_data.csv")
  il6 <- vector(length=length(reg_vals))
  treatment <- vector(length=length(reg_vals))
  
  for(j in 1:length(reg_vals))
  {
    cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == names(reg_vals)[j])
    cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
    il6[j] <- as.numeric(cyt_subset2$IL6_FC)
    treatment[j] <- unique(meta_data$vaccination[which(meta_data$donor == names(reg_vals)[j])])
    print(j)
  }
  
  val <- cor.test(reg_vals, il6, method="spearman")
  pval_il6[i] <- val$p.value
  rho_il6[i] <- val$estimate
  
  df_reg <- data.frame(
    score <- reg_vals,
    c <- il6,
    t <- treatment
  )
  colnames(df_reg) <- c("score", "c", "Group")
  
  p2 <- ggplot(df_reg, aes(x=score, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_classic() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
    geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab(paste0(colnames(logfc_mat_final)[i]," score")) +
    ylab("FC IL6 by PBMCs") + ggtitle(paste0(colnames(logfc_mat_final)[i]," score vs. IL6"))+
    annotate("text", x=-0.6, y=2, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_il6[i],6))) + scale_fill_manual(values=c("#F3A546", "#00A2B3"))
  plot_list2[[i]] <- p2
  
  # ## i == 25
  # if(colnames(logfc_mat_final)[i] == "KLF5"){
  #   file_name = paste0(OUT_DIR,"KLF5_hsca_corr_IL6.tiff")
  #   tiff(file_name, units="in", width=4, height=3, res=300)
  #   ggplot(df_reg, aes(x=score, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_classic() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  #     geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab(paste0(colnames(logfc_mat_final)[i]," score")) +
  #     ylab("FC IL6 by PBMCs") + ggtitle(paste0(colnames(logfc_mat_final)[i]," score vs. IL6"))+
  #     annotate("text", x=-0.4, y=2, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_il6[i],4))) +  scale_fill_manual(values=c("#F3A546", "grey")) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
  #   dev.off()
  # }
  # 
  # 
  # ## i == 9
  # if(colnames(logfc_mat_final)[i] == "EGR1"){
  #   file_name = paste0(OUT_DIR,"EGR1_hsca_corr_IL6.tiff")
  #   tiff(file_name, units="in", width=4, height=3, res=300)
  #   ggplot(df_reg, aes(x=score, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_classic() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  #     geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab(paste0(colnames(logfc_mat_final)[i]," score")) +
  #     ylab("FC IL6 by PBMCs") + ggtitle(paste0(colnames(logfc_mat_final)[i]," score vs. IL6"))+
  #     annotate("text", x=-0.4, y=2, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval_il6[i],4))) +  scale_fill_manual(values=c("#F3A546", "grey")) + theme(axis.text=element_text(size=12),axis.title=element_text(size=14))
  #   dev.off()
  # }
}


for(k in 1:length(plot_list1)) 
{
  #file_name = paste0(OUT_DIR, colnames(logfc_mat_final)[k], "_hsca_corr_IL1B.tiff")
  #tiff(file_name, units="in", width=3.5, height=3, res=300)
  #print(plot_list1[[k]])
  #dev.off()
  
  p <- plot_list1[[k]]
  final_file <- file.path(paste0(OUT_DIR, colnames(logfc_mat_final)[k], "_hsca_corr_IL1B.pdf"))
  ggsave(filename = final_file, plot = p, height = 3.2, width = 4)
  
  print(k)
  
}

for(k in 1:length(plot_list2)) 
{
  p <- plot_list2[[k]]
  final_file <- file.path(paste0(OUT_DIR, colnames(logfc_mat_final)[k], "_hsca_corr_IL6.pdf"))
  ggsave(filename = final_file, plot = p, height = 3.2, width = 4)
  
  print(k)
  
}

#IL1B and IL6 rho comparison
  
df_comparison <- data.frame(
  val <- c(abs(rho_il1b), abs(rho_il6)),
  g <- c(rep("IL1B", length(rho_il1b)), rep("IL6", length(rho_il6)))
)

tiff(paste0(OUT_DIR, "IL6_IL1B_corr_violin.tiff"), units="in", width=4, height=4, res=300)
ggplot(df_comparison, aes(x=g, y=val, fill=g)) + 
  geom_violin(trim=FALSE)  + geom_jitter(shape=16,size=2,alpha=0.6, position=position_jitter(0.2)) + theme_classic() + labs(x="",y="Rho") + scale_fill_manual(values=c("#F1788D","#F3A546"))+
  theme(legend.position="none")
dev.off()

tiff(paste0(OUT_DIR, "IL6_IL1B_corr_boxplot.tiff"), units="in", width=3.2, height=3.2, res=300)
ggplot(df_comparison, aes(x=g, y=val, fill=g)) + 
  geom_boxplot()  + geom_jitter(shape=21,size=2,alpha=0.5, position=position_jitter(0.2), color="black", fill="grey") + theme_classic() + labs(x="",y="Rho") + scale_fill_manual(values=c("#F1788D","#F3A546"))+
  theme(legend.position="none") + theme(axis.text=element_text(size=12),
                                        axis.title=element_text(size=14))
dev.off()

wilcox.test(df_comparison$val[which(df_comparison$g == "IL1B")], df_comparison$val[which(df_comparison$g == "IL6")])


