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

# In this script we select all cells that were labelled as HSC_a (new lab HSC c1) or HSC_b (new lab HSC c2) and we re-cluster the cells to form high-resolution HSC subclusters

# setup ------------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/HSC_subcluster_analyses")
IN_DIR <- "../label_transfer/"
OUT_DIR <- "0_SCT_cluster/"
dir.create(OUT_DIR)


# Read in final labeled obj ----------------------------------------------
#S9 Td0 needs to be filtered out from this object!
all <- readRDS(file=paste0(IN_DIR,"allCells_integrated_label_transfer_UNIQUE_IDS.rds"))
DefaultAssay(all) <- "RNA"

#Filter out S9 Td0 
unique_ids <- unique(all$UNIQUE_ID)
unique_ids = unique_ids[!(unique_ids %in% "SNG-LB-SS-1S-RS-S9-CD34neg_S8_R1_001_Td0_BCG" )]
all<- subset(
  x = all,
  subset = UNIQUE_ID %in% unique_ids
)

# Get cell barcodes of all HSCs -------------------------------------------
all_meta <- all@meta.data
HSC_barcodes <- row.names(all_meta)[which(all_meta$clust_names %in% c('HSC_a', 'HSC_b'))]


# read in merged_filtered seurat object from initial processing --------------------------
# We will perform cluster on the HSCs from this object ------------------------------------
mergeAllCells <- readRDS(file = "../../1_data_processing_clustering/merged_filtered_seuratObject.rds")
mergeAllCells$batchID <- factor(mergeAllCells$batchID)

# Subset to keep only HSC barcodes ------------------------------------------------
mergeAllCells$temp <- row.names(mergeAllCells@meta.data)
HSC_obj <- subset(
  x = mergeAllCells,
  subset = temp %in% HSC_barcodes
)

HSC_obj$temp <- NULL


# Split the merged object by timepoint ------------------------------------
HSC_obj_list <- SplitObject(HSC_obj, split.by = "timepoint")
#   Td0     Tm3
#  8,480   12,110

Tm3 <- HSC_obj_list[[1]]
Td0 <- HSC_obj_list[[2]]


# Split Tm3 into CTL and BCG ----------------------------------------------
Tm3_list <- SplitObject(Tm3, split.by = "vaccination")
#   BCG       CTL
#  9,144     2,966

Tm3_BCG <- Tm3_list[[1]]
Tm3_CTL <- Tm3_list[[2]]


# Create final seurat object list split 3 ways ----------------------------
mergeAllCells_list <- c(Td0, Tm3_BCG, Tm3_CTL)


# Run scTransform separately  ---------------------------------------------
#scTransform replaces NormalizeData(), ScaleData(), and FindVariableFeatures(); FindvariableFeatures() identifies features that are outliers on a mean-variance plot
mergeAllCells_list <- lapply(X = mergeAllCells_list, FUN = function(x) {
  x <- SCTransform(x, vars.to.regress = c("percent.mt","batchID"))
})
saveRDS(mergeAllCells_list, file = paste0(OUT_DIR,"HSC_obj_afterSCTransform.rds"))


# Integration -------------------------------------------------------------
features <- SelectIntegrationFeatures(object.list = mergeAllCells_list, nfeatures = 3000)
mergeAllCells_list <- PrepSCTIntegration(object.list = mergeAllCells_list, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = mergeAllCells_list, normalization.method = "SCT", 
                                  anchor.features = features)

allCells_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT", verbose=TRUE)
DefaultAssay(allCells_integrated) <- "integrated"

saveRDS(allCells_integrated, file = paste0(OUT_DIR,"HSC_obj_SCT_withIntegration.rds"))



# Run standard clustering workflow ----------------------------------------

allCells_integrated <- RunPCA(allCells_integrated, npcs = 30)
allCells_integrated <- RunUMAP(allCells_integrated, dims = 1:30, seed.use = 2020)
allCells_integrated <- FindNeighbors(allCells_integrated, dims = 1:20)
allCells_integrated <- FindClusters(allCells_integrated, resolution = 0.5)

saveRDS(allCells_integrated, file = paste0(OUT_DIR,"HSC_obj_SCT_integrated_withClusters.rds"))


