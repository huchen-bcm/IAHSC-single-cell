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
OUT_DIR <- "get_diff_regulons_outs/"
dir.create(OUT_DIR)

## Get data from sce object:
seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"



#Analysis excluding S2, S3, S12
auc_mat <- read.csv(file=paste0("scenic_cluster_HSC_b_auc_mtx.csv"))
df <- seurat_obj@meta.data
samples <- unique(df$UNIQUE_ID)
for(i in 1:length(samples))
{
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
for(i in 1:length(Tm3_all))
{
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

#To make all plots with p < 0.06
pval <- vector(length=length(row.names(logfc_mat)))
pval_wilcox <- vector(length=length(row.names(logfc_mat)))
tf_name <- vector()
for(i in 1:length(row.names(logfc_mat)))
{
  values <- logfc_mat[i, ]
  BCG_donors <- unique(meta_data$donor[which(meta_data$vaccination == "BCG")])
  CTL_donors <- unique(meta_data$donor[which(meta_data$vaccination == "CTL")])
  CTL_vals <- values[which(names(values) %in% CTL_donors)]
  BCG_vals <- values[which(names(values) %in% BCG_donors)]
  
  res <- t.test(BCG_vals, CTL_vals)
  pval[i] <- res$p.value
  
  res_wilcox <- wilcox.test(BCG_vals, CTL_vals)
  pval_wilcox[i] <- res_wilcox$p.value
  }
  
}
padj_wilcox <- p.adjust(pval_wilcox, method='BH')
wilcox_data <- data.frame(cbind(pval_wilcox, padj_wilcox))
row.names(wilcox_data) <- row.names(logfc_mat)
write.csv(wilcox_data, file=paste0(OUT_DIR, "HSC_b_wilcox_data_no_outliers.csv")) ##RENAME HSC c2

