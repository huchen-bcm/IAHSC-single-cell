#R version 4.1.0
library(ggplot2)
library(cowplot)
library(reticulate)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(Seurat)
library(sctransform)
library(data.table)

# In this script we determine all the sample-cluster groups that have fewer than 15 cells
# These are not included in future analyses 


# setup -------------------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/HSC_subcluster_analyses/") 
OUT_DIR <- "3_filter15/"
dir.create(OUT_DIR)


# Determine the number of cells from each sample in each term fate group ----------
HSC <- readRDS(file=paste0("1_UMAP_label_clusters/HSC_final_UNIQUE_IDS.rds"))
DefaultAssay(HSC) <- "integrated"

meta_data <- read.csv(file="all_samples_meta_data.csv")
allCells_by_cluster_list <- SplitObject(HSC, split.by = "seurat_clusters")

for(i in 1:length(allCells_by_cluster_list)) 
{
  name <- names(allCells_by_cluster_list)[i]
  obj <- allCells_by_cluster_list[[i]]
  md <- obj@meta.data
  
  obj_split <- SplitObject(obj, split.by = "UNIQUE_ID")
  sample_nums <- vector(length=length(obj_split))
  sample_names <- vector(length=length(obj_split))
  timepoint <- vector(length=length(obj_split))
  vaccination <- vector(length=length(obj_split))
  donor <- vector(length=length(obj_split))
  for(j in 1:length(obj_split))
  {
    obj_temp <- obj_split[[j]]
    sample_nums[j] <- dim(obj_temp)[2]
    sample_names[j] <- names(obj_split)[j]
    timepoint[j] <- subset(meta_data, meta_data$Sample %in% sample_names[j])$timepoint
    vaccination[j] <- subset(meta_data, meta_data$Sample %in% sample_names[j])$vaccination
    donor[j] <- subset(meta_data, meta_data$Sample %in% sample_names[j])$donor
  }
  df <- data.frame(cbind(sample_names, sample_nums, timepoint, vaccination, donor))
  assign(paste0("c",name,"_df"),df)
  print(i)
}

cell_df_list <- list(c0_df, c1_df, c2_df, c3_df, c4_df, c5_df, c6_df, c7_df, c8_df, c9_df)
clusters <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9")


#filer15 list
for(i in 1:length(clusters))
{
  df <- cell_df_list[[i]]
  df_subset <- subset(df, as.numeric(df$sample_nums) <15)
  names <- df_subset$sample_names
  assign(paste0(clusters[i],"_filter15"), names)
}
filter_15_list <- list(c0_filter15, c1_filter15, c2_filter15, c3_filter15, c4_filter15, c5_filter15, c6_filter15, c7_filter15, c8_filter15, c9_filter15)
saveRDS(filter_15_list, file=paste0(OUT_DIR,"filter_15_list.rds"))
filter_15_list <- readRDS(file=paste0(OUT_DIR,"filter_15_list.rds"))
