#R version 4.1.0
install.packages("textTinyR")

library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(Seurat)
library(textTinyR)
library(pbapply)

# setup --------------------------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/HSC_subcluster_analyses")
IN_DIR <- "1_UMAP_label_clusters/"
out_dir <- "2_pseudobulk/"
dir.create(out_dir)

# Read in processed object with UNIQUE ID labels ------------------------
all <- readRDS(file=paste0(IN_DIR,"HSC_final_UNIQUE_IDS.rds"))
DefaultAssay(all) <- "RNA"


# Split the object by celltype ---------------------------------------------
all_list <- SplitObject(all, split.by = "seurat_clusters")



# Save a separate object for each cluster ---------------------------------
for (i in 1:length(all_list)){
  name <- names(all_list)[i]
  saveRDS(all_list[[i]], file = paste0(out_dir,name,"_singlets.rds"))
  print(i)
}


# Calculate pseudobulk sums -----------------------------------------------
clusters <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9")
for (i in 1:length(clusters)){
  name <- clusters[i]
  
  ## read in single cell data (raw UMI counts stored in RNA slot) for each cluster
  dat <- readRDS(paste0(out_dir,name,"_singlets.rds"))
  ## get raw UMI counts
  raw_data_sparse <- GetAssayData(dat, assay = "RNA", slot = "counts")
  meta_data <- dat@meta.data
  sample_colname <- "UNIQUE_ID"
  
  IDs <- as.data.frame(meta_data)[, sample_colname]
  unique_ID_list <- as.list(unique(IDs))
  
  ## calculate pseudobulk by summing across all cells identified in each individual,cdt pair
  pseudobulk <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Sums(raw_data_sparse[,IDs == x, drop = FALSE], rowSums = TRUE)}))
  cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(raw_data_sparse[,IDs == x, drop = FALSE])})
  colnames(pseudobulk) <- names(cellcount) <- unique_ID_list
  rownames(pseudobulk) <- rownames(x = dat)
  
  saveRDS(pseudobulk, file = paste0(out_dir,name,"_pseudobulk.rds"))
  
  print(i)
}
