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

setwd("/scRNA_analyses/2_main_analyses/pyscenic")

## Get data from sce object:
seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"

#each individual cluster
all_list <- SplitObject(seurat_obj, split.by = "clust_names")
names(all_list)

for(i in 1:length(all_list))
{
  obj <- all_list[[i]]
  name <- names(all_list)[i]
  sce <- as.SingleCellExperiment(obj)
  exprMat <- counts(sce)
  
  
  cellInfo <- data.frame(obj@meta.data)
  
  
  loom <- build_loom(paste0(name, "_humanBM.loom"), dgem=exprMat)
  loom <- add_cell_annotation(loom, cellInfo)
  close_loom(loom)
  
}

#HSCa and HSCb combined
HSC_obj<- subset(
  x = seurat_obj,
  subset = clust_names %in% c("HSC_a", "HSC_b")
)

sce <- as.SingleCellExperiment(HSC_obj)
exprMat <- counts(sce)

cellInfo <- data.frame(HSC_obj@meta.data)


loom <- build_loom(paste0("HSC_a_b_humanBM.loom"), dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)



#HSCa, HSC_b, CMP_b combined
HSC_obj<- subset(
  x = seurat_obj,
  subset = clust_names %in% c("HSC_a", "HSC_b", "CMP_b")
)

final_obj<- subset(
  x = HSC_obj,
  subset = timepoint %in% c("Td0")
)  ##15k cells

sce <- as.SingleCellExperiment(final_obj)
exprMat <- counts(sce)

cellInfo <- data.frame(final_obj@meta.data)


loom <- build_loom(paste0("HSCab_CMPb_humanBM.loom"), dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)


#CMPb_GMP combined 

CMP_GMP_obj<- subset(
  x = seurat_obj,
  subset = clust_names %in% c("CMP_b", "GMP_a", "GMP_b")
)

final_obj<- subset(
  x = CMP_GMP_obj,
  subset = timepoint %in% c("Td0")
)  ##13k cells

sce <- as.SingleCellExperiment(final_obj)
exprMat <- counts(sce)

cellInfo <- data.frame(final_obj@meta.data)


loom <- build_loom(paste0("CMPb_GMPab_humanBM.loom"), dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)


#CMPb_GMPa combined 

CMP_GMP_obj<- subset(
  x = seurat_obj,
  subset = clust_names %in% c("CMP_b", "GMP_a")
)

final_obj<- subset(
  x = CMP_GMP_obj,
  subset = timepoint %in% c("Td0")
)  ##10k cells

sce <- as.SingleCellExperiment(final_obj)
exprMat <- counts(sce)

cellInfo <- data.frame(final_obj@meta.data)


loom <- build_loom(paste0("CMPb_GMPa_humanBM.loom"), dgem=exprMat)
loom <- add_cell_annotation(loom, cellInfo)
close_loom(loom)
