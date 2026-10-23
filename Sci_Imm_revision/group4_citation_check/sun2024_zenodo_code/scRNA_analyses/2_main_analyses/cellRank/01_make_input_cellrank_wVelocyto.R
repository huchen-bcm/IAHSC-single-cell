if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pcaMethods")
install_github("velocyto-team/velocyto.R")
install.packages("pagoda2")

library(devtools)
library(velocyto.R)
library(pagoda2)
library(Seurat)
library(hdf5r)
library(SeuratDisk)
library(SingleCellExperiment)

# Here we read in loom files for each scRNA capture (14 total captures) that were generated from running velocyto separately for each capture 
# For each donor within each capture: we extract cells and their embeddings + cluster IDs
# This is the input needed to run CellRank downstream 



setwd("/scRNA_analyses/2_main_analyses/cellRank")



# Capture 1 ---------------------------------------------------------------
IN_DIR <- "/capture1_output/velocyto/" 
OUT_DIR <- "wVelocity_input/"
dir.create(OUT_DIR)

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture1_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"

S8_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S8-CD34neg_S7_R1_001_Td0_BCG"
)
cells_S8_Td0 <- Cells(S8_Td0_obj)
cells_S8_Td0 <- gsub("capture1_", "capture1_output:", cells_S8_Td0)
cells_S8_Td0 <- gsub("-1", "x", cells_S8_Td0)
write.csv(cells_S8_Td0, file = paste0(OUT_DIR, "cellID_obs_S8_Td0.csv"), row.names = FALSE)

embed_S8_Td0 <- Embeddings(S8_Td0_obj, reduction = "umap")
embed_names_S8_Td0 <- row.names(embed_S8_Td0)
embed_names_S8_Td0 <- gsub("capture1_", "capture1_output:", embed_names_S8_Td0)
embed_names_S8_Td0 <- gsub("-1", "x", embed_names_S8_Td0)
row.names(embed_S8_Td0) <- embed_names_S8_Td0
write.csv(embed_S8_Td0, file = paste0(OUT_DIR, "cell_embeddings_S8_Td0.csv"))

clusters_S8_Td0 <- data.frame(S8_Td0_obj@meta.data$clust_names)
row.names(clusters_S8_Td0) <- embed_names_S8_Td0
write.csv(clusters_S8_Td0, file = paste0(OUT_DIR, "clusters_S8_Td0.csv"))



# Capture 2 ---------------------------------------------------------------


IN_DIR <- "/capture2_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture2_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S8_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S8-CD34neg_S7_R1_001_Tm3_BCG"
)

cells_S8_Tm3 <- Cells(S8_Tm3_obj)
cells_S8_Tm3 <- gsub("capture2_", "capture2_output:", cells_S8_Tm3)
cells_S8_Tm3 <- gsub("-1", "x", cells_S8_Tm3)
write.csv(cells_S8_Tm3, file = paste0(OUT_DIR, "cellID_obs_S8_Tm3.csv"), row.names = FALSE)

embed_S8_Tm3 <- Embeddings(S8_Tm3_obj, reduction = "umap")
embed_names_S8_Tm3 <- row.names(embed_S8_Tm3)
embed_names_S8_Tm3 <- gsub("capture2_", "capture2_output:", embed_names_S8_Tm3)
embed_names_S8_Tm3 <- gsub("-1", "x", embed_names_S8_Tm3)
row.names(embed_S8_Tm3) <- embed_names_S8_Tm3
write.csv(embed_S8_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S8_Tm3.csv"))

clusters_S8_Tm3 <- data.frame(S8_Tm3_obj@meta.data$clust_names)
row.names(clusters_S8_Tm3) <- embed_names_S8_Tm3
write.csv(clusters_S8_Tm3, file = paste0(OUT_DIR, "clusters_S8_Tm3.csv"))


# Capture 3 ---------------------------------------------------------------

IN_DIR <- "/capture3_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture3_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S20_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S20-CD34neg_S19_R1_001_Td0_BCG"
)

S15_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S15-CD34neg_S14_R1_001_Td0_BCG"
)

cells_S20_Td0 <- Cells(S20_Td0_obj)
cells_S20_Td0 <- gsub("capture3_", "capture3_output:", cells_S20_Td0)
cells_S20_Td0 <- gsub("-1", "x", cells_S20_Td0)
write.csv(cells_S20_Td0, file = paste0(OUT_DIR, "cellID_obs_S20_Td0.csv"), row.names = FALSE)

cells_S15_Td0 <- Cells(S15_Td0_obj)
cells_S15_Td0 <- gsub("capture3_", "capture3_output:", cells_S15_Td0)
cells_S15_Td0 <- gsub("-1", "x", cells_S15_Td0)
write.csv(cells_S15_Td0, file = paste0(OUT_DIR, "cellID_obs_S15_Td0.csv"), row.names = FALSE)


embed_S20_Td0 <- Embeddings(S20_Td0_obj, reduction = "umap")
embed_names_S20_Td0 <- row.names(embed_S20_Td0)
embed_names_S20_Td0 <- gsub("capture3_", "capture3_output:", embed_names_S20_Td0)
embed_names_S20_Td0 <- gsub("-1", "x", embed_names_S20_Td0)
row.names(embed_S20_Td0) <- embed_names_S20_Td0
write.csv(embed_S20_Td0, file = paste0(OUT_DIR, "cell_embeddings_S20_Td0.csv"))

embed_S15_Td0 <- Embeddings(S15_Td0_obj, reduction = "umap")
embed_names_S15_Td0 <- row.names(embed_S15_Td0)
embed_names_S15_Td0 <- gsub("capture3_", "capture3_output:", embed_names_S15_Td0)
embed_names_S15_Td0 <- gsub("-1", "x", embed_names_S15_Td0)
row.names(embed_S15_Td0) <- embed_names_S15_Td0
write.csv(embed_S15_Td0, file = paste0(OUT_DIR, "cell_embeddings_S15_Td0.csv"))


clusters_S20_Td0 <- data.frame(S20_Td0_obj@meta.data$clust_names)
row.names(clusters_S20_Td0) <- embed_names_S20_Td0
write.csv(clusters_S20_Td0, file = paste0(OUT_DIR, "clusters_S20_Td0.csv"))

clusters_S15_Td0 <- data.frame(S15_Td0_obj@meta.data$clust_names)
row.names(clusters_S15_Td0) <- embed_names_S15_Td0
write.csv(clusters_S15_Td0, file = paste0(OUT_DIR, "clusters_S15_Td0.csv"))





# Capture 4 ---------------------------------------------------------------

IN_DIR <- "/capture4_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture4_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S20_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S20-CD34neg_S19_R1_001_Tm3_BCG"
)

S15_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S15-CD34neg_S14_R1_001_Tm3_BCG"
)

cells_S20_Tm3 <- Cells(S20_Tm3_obj)
cells_S20_Tm3 <- gsub("capture4_", "capture4_output:", cells_S20_Tm3)
cells_S20_Tm3 <- gsub("-1", "x", cells_S20_Tm3)
write.csv(cells_S20_Tm3, file = paste0(OUT_DIR, "cellID_obs_S20_Tm3.csv"), row.names = FALSE)

cells_S15_Tm3 <- Cells(S15_Tm3_obj)
cells_S15_Tm3 <- gsub("capture4_", "capture4_output:", cells_S15_Tm3)
cells_S15_Tm3 <- gsub("-1", "x", cells_S15_Tm3)
write.csv(cells_S15_Tm3, file = paste0(OUT_DIR, "cellID_obs_S15_Tm3.csv"), row.names = FALSE)


embed_S20_Tm3 <- Embeddings(S20_Tm3_obj, reduction = "umap")
embed_names_S20_Tm3 <- row.names(embed_S20_Tm3)
embed_names_S20_Tm3 <- gsub("capture4_", "capture4_output:", embed_names_S20_Tm3)
embed_names_S20_Tm3 <- gsub("-1", "x", embed_names_S20_Tm3)
row.names(embed_S20_Tm3) <- embed_names_S20_Tm3
write.csv(embed_S20_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S20_Tm3.csv"))

embed_S15_Tm3 <- Embeddings(S15_Tm3_obj, reduction = "umap")
embed_names_S15_Tm3 <- row.names(embed_S15_Tm3)
embed_names_S15_Tm3 <- gsub("capture4_", "capture4_output:", embed_names_S15_Tm3)
embed_names_S15_Tm3 <- gsub("-1", "x", embed_names_S15_Tm3)
row.names(embed_S15_Tm3) <- embed_names_S15_Tm3
write.csv(embed_S15_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S15_Tm3.csv"))


clusters_S20_Tm3 <- data.frame(S20_Tm3_obj@meta.data$clust_names)
row.names(clusters_S20_Tm3) <- embed_names_S20_Tm3
write.csv(clusters_S20_Tm3, file = paste0(OUT_DIR, "clusters_S20_Tm3.csv"))

clusters_S15_Tm3 <- data.frame(S15_Tm3_obj@meta.data$clust_names)
row.names(clusters_S15_Tm3) <- embed_names_S15_Tm3
write.csv(clusters_S15_Tm3, file = paste0(OUT_DIR, "clusters_S15_Tm3.csv"))




# Capture 5 ---------------------------------------------------------------

IN_DIR <- "/capture5_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture5_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S21_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S21-CD34neg_S20_R1_001_Td0_BCG"
)

S13_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S13-CD34neg_S12_R1_001_Td0_CTL"
)

S16_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S16-CD34neg_S15_R1_001_Td0_BCG"
)


cells_S21_Td0 <- Cells(S21_Td0_obj)
cells_S21_Td0 <- gsub("capture5_", "capture5_output:", cells_S21_Td0)
cells_S21_Td0 <- gsub("-1", "x", cells_S21_Td0)
write.csv(cells_S21_Td0, file = paste0(OUT_DIR, "cellID_obs_S21_Td0.csv"), row.names = FALSE)

cells_S13_Td0 <- Cells(S13_Td0_obj)
cells_S13_Td0 <- gsub("capture5_", "capture5_output:", cells_S13_Td0)
cells_S13_Td0 <- gsub("-1", "x", cells_S13_Td0)
write.csv(cells_S13_Td0, file = paste0(OUT_DIR, "cellID_obs_S13_Td0.csv"), row.names = FALSE)

cells_S16_Td0 <- Cells(S16_Td0_obj)
cells_S16_Td0 <- gsub("capture5_", "capture5_output:", cells_S16_Td0)
cells_S16_Td0 <- gsub("-1", "x", cells_S16_Td0)
write.csv(cells_S16_Td0, file = paste0(OUT_DIR, "cellID_obs_S16_Td0.csv"), row.names = FALSE)

embed_S21_Td0 <- Embeddings(S21_Td0_obj, reduction = "umap")
embed_names_S21_Td0 <- row.names(embed_S21_Td0)
embed_names_S21_Td0 <- gsub("capture5_", "capture5_output:", embed_names_S21_Td0)
embed_names_S21_Td0 <- gsub("-1", "x", embed_names_S21_Td0)
row.names(embed_S21_Td0) <- embed_names_S21_Td0
write.csv(embed_S21_Td0, file = paste0(OUT_DIR, "cell_embeddings_S21_Td0.csv"))

embed_S13_Td0 <- Embeddings(S13_Td0_obj, reduction = "umap")
embed_names_S13_Td0 <- row.names(embed_S13_Td0)
embed_names_S13_Td0 <- gsub("capture5_", "capture5_output:", embed_names_S13_Td0)
embed_names_S13_Td0 <- gsub("-1", "x", embed_names_S13_Td0)
row.names(embed_S13_Td0) <- embed_names_S13_Td0
write.csv(embed_S13_Td0, file = paste0(OUT_DIR, "cell_embeddings_S13_Td0.csv"))

embed_S16_Td0 <- Embeddings(S16_Td0_obj, reduction = "umap")
embed_names_S16_Td0 <- row.names(embed_S16_Td0)
embed_names_S16_Td0 <- gsub("capture5_", "capture5_output:", embed_names_S16_Td0)
embed_names_S16_Td0 <- gsub("-1", "x", embed_names_S16_Td0)
row.names(embed_S16_Td0) <- embed_names_S16_Td0
write.csv(embed_S16_Td0, file = paste0(OUT_DIR, "cell_embeddings_S16_Td0.csv"))


clusters_S21_Td0 <- data.frame(S21_Td0_obj@meta.data$clust_names)
row.names(clusters_S21_Td0) <- embed_names_S21_Td0
write.csv(clusters_S21_Td0, file = paste0(OUT_DIR, "clusters_S21_Td0.csv"))

clusters_S13_Td0 <- data.frame(S13_Td0_obj@meta.data$clust_names)
row.names(clusters_S13_Td0) <- embed_names_S13_Td0
write.csv(clusters_S13_Td0, file = paste0(OUT_DIR, "clusters_S13_Td0.csv"))

clusters_S16_Td0 <- data.frame(S16_Td0_obj@meta.data$clust_names)
row.names(clusters_S16_Td0) <- embed_names_S16_Td0
write.csv(clusters_S16_Td0, file = paste0(OUT_DIR, "clusters_S16_Td0.csv"))


# Capture 6 ---------------------------------------------------------------

IN_DIR <- "/capture6_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture6_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S21_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S21-CD34neg_S20_R1_001_Tm3_BCG"
)

S13_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S13-CD34neg_S12_R1_001_Tm3_CTL"
)

S16_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S16-CD34neg_S15_R1_001_Tm3_BCG"
)


cells_S21_Tm3 <- Cells(S21_Tm3_obj)
cells_S21_Tm3 <- gsub("capture6_", "capture6_output:", cells_S21_Tm3)
cells_S21_Tm3 <- gsub("-1", "x", cells_S21_Tm3)
write.csv(cells_S21_Tm3, file = paste0(OUT_DIR, "cellID_obs_S21_Tm3.csv"), row.names = FALSE)

cells_S13_Tm3 <- Cells(S13_Tm3_obj)
cells_S13_Tm3 <- gsub("capture6_", "capture6_output:", cells_S13_Tm3)
cells_S13_Tm3 <- gsub("-1", "x", cells_S13_Tm3)
write.csv(cells_S13_Tm3, file = paste0(OUT_DIR, "cellID_obs_S13_Tm3.csv"), row.names = FALSE)

cells_S16_Tm3 <- Cells(S16_Tm3_obj)
cells_S16_Tm3 <- gsub("capture6_", "capture6_output:", cells_S16_Tm3)
cells_S16_Tm3 <- gsub("-1", "x", cells_S16_Tm3)
write.csv(cells_S16_Tm3, file = paste0(OUT_DIR, "cellID_obs_S16_Tm3.csv"), row.names = FALSE)

embed_S21_Tm3 <- Embeddings(S21_Tm3_obj, reduction = "umap")
embed_names_S21_Tm3 <- row.names(embed_S21_Tm3)
embed_names_S21_Tm3 <- gsub("capture6_", "capture6_output:", embed_names_S21_Tm3)
embed_names_S21_Tm3 <- gsub("-1", "x", embed_names_S21_Tm3)
row.names(embed_S21_Tm3) <- embed_names_S21_Tm3
write.csv(embed_S21_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S21_Tm3.csv"))

embed_S13_Tm3 <- Embeddings(S13_Tm3_obj, reduction = "umap")
embed_names_S13_Tm3 <- row.names(embed_S13_Tm3)
embed_names_S13_Tm3 <- gsub("capture6_", "capture6_output:", embed_names_S13_Tm3)
embed_names_S13_Tm3 <- gsub("-1", "x", embed_names_S13_Tm3)
row.names(embed_S13_Tm3) <- embed_names_S13_Tm3
write.csv(embed_S13_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S13_Tm3.csv"))

embed_S16_Tm3 <- Embeddings(S16_Tm3_obj, reduction = "umap")
embed_names_S16_Tm3 <- row.names(embed_S16_Tm3)
embed_names_S16_Tm3 <- gsub("capture6_", "capture6_output:", embed_names_S16_Tm3)
embed_names_S16_Tm3 <- gsub("-1", "x", embed_names_S16_Tm3)
row.names(embed_S16_Tm3) <- embed_names_S16_Tm3
write.csv(embed_S16_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S16_Tm3.csv"))


clusters_S21_Tm3 <- data.frame(S21_Tm3_obj@meta.data$clust_names)
row.names(clusters_S21_Tm3) <- embed_names_S21_Tm3
write.csv(clusters_S21_Tm3, file = paste0(OUT_DIR, "clusters_S21_Tm3.csv"))

clusters_S13_Tm3 <- data.frame(S13_Tm3_obj@meta.data$clust_names)
row.names(clusters_S13_Tm3) <- embed_names_S13_Tm3
write.csv(clusters_S13_Tm3, file = paste0(OUT_DIR, "clusters_S13_Tm3.csv"))

clusters_S16_Tm3 <- data.frame(S16_Tm3_obj@meta.data$clust_names)
row.names(clusters_S16_Tm3) <- embed_names_S16_Tm3
write.csv(clusters_S16_Tm3, file = paste0(OUT_DIR, "clusters_S16_Tm3.csv"))



# Capture 7 ---------------------------------------------------------------

IN_DIR <- "/capture7_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture7_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S5_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S5-CD34neg_S4_R1_001_Td0_BCG"
)

S6_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S6-CD34neg_S5_R1_001_Td0_BCG"
)

S7_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S7-CD34neg_S6_R1_001_Td0_CTL"
)


cells_S5_Td0 <- Cells(S5_Td0_obj)
cells_S5_Td0 <- gsub("capture7_", "capture7_output:", cells_S5_Td0)
cells_S5_Td0 <- gsub("-1", "x", cells_S5_Td0)
write.csv(cells_S5_Td0, file = paste0(OUT_DIR, "cellID_obs_S5_Td0.csv"), row.names = FALSE)

cells_S6_Td0 <- Cells(S6_Td0_obj)
cells_S6_Td0 <- gsub("capture7_", "capture7_output:", cells_S6_Td0)
cells_S6_Td0 <- gsub("-1", "x", cells_S6_Td0)
write.csv(cells_S6_Td0, file = paste0(OUT_DIR, "cellID_obs_S6_Td0.csv"), row.names = FALSE)

cells_S7_Td0 <- Cells(S7_Td0_obj)
cells_S7_Td0 <- gsub("capture7_", "capture7_output:", cells_S7_Td0)
cells_S7_Td0 <- gsub("-1", "x", cells_S7_Td0)
write.csv(cells_S7_Td0, file = paste0(OUT_DIR, "cellID_obs_S7_Td0.csv"), row.names = FALSE)

embed_S5_Td0 <- Embeddings(S5_Td0_obj, reduction = "umap")
embed_names_S5_Td0 <- row.names(embed_S5_Td0)
embed_names_S5_Td0 <- gsub("capture7_", "capture7_output:", embed_names_S5_Td0)
embed_names_S5_Td0 <- gsub("-1", "x", embed_names_S5_Td0)
row.names(embed_S5_Td0) <- embed_names_S5_Td0
write.csv(embed_S5_Td0, file = paste0(OUT_DIR, "cell_embeddings_S5_Td0.csv"))

embed_S6_Td0 <- Embeddings(S6_Td0_obj, reduction = "umap")
embed_names_S6_Td0 <- row.names(embed_S6_Td0)
embed_names_S6_Td0 <- gsub("capture7_", "capture7_output:", embed_names_S6_Td0)
embed_names_S6_Td0 <- gsub("-1", "x", embed_names_S6_Td0)
row.names(embed_S6_Td0) <- embed_names_S6_Td0
write.csv(embed_S6_Td0, file = paste0(OUT_DIR, "cell_embeddings_S6_Td0.csv"))
  
embed_S7_Td0 <- Embeddings(S7_Td0_obj, reduction = "umap")
embed_names_S7_Td0 <- row.names(embed_S7_Td0)
embed_names_S7_Td0 <- gsub("capture7_", "capture7_output:", embed_names_S7_Td0)
embed_names_S7_Td0 <- gsub("-1", "x", embed_names_S7_Td0)
row.names(embed_S7_Td0) <- embed_names_S7_Td0
write.csv(embed_S7_Td0, file = paste0(OUT_DIR, "cell_embeddings_S7_Td0.csv"))


clusters_S5_Td0 <- data.frame(S5_Td0_obj@meta.data$clust_names)
row.names(clusters_S5_Td0) <- embed_names_S5_Td0
write.csv(clusters_S5_Td0, file = paste0(OUT_DIR, "clusters_S5_Td0.csv"))

clusters_S6_Td0 <- data.frame(S6_Td0_obj@meta.data$clust_names)
row.names(clusters_S6_Td0) <- embed_names_S6_Td0
write.csv(clusters_S6_Td0, file = paste0(OUT_DIR, "clusters_S6_Td0.csv"))

clusters_S7_Td0 <- data.frame(S7_Td0_obj@meta.data$clust_names)
row.names(clusters_S7_Td0) <- embed_names_S7_Td0
write.csv(clusters_S7_Td0, file = paste0(OUT_DIR, "clusters_S7_Td0.csv"))



# Capture 8 ---------------------------------------------------------------
IN_DIR <- "/capture8_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture8_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S5_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S5-CD34neg_S4_R1_001_Tm3_BCG"
)

S6_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S6-CD34neg_S5_R1_001_Tm3_BCG"
)

S7_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S7-CD34neg_S6_R1_001_Tm3_CTL"
)


cells_S5_Tm3 <- Cells(S5_Tm3_obj)
cells_S5_Tm3 <- gsub("capture8_", "capture8_output:", cells_S5_Tm3)
cells_S5_Tm3 <- gsub("-1", "x", cells_S5_Tm3)
write.csv(cells_S5_Tm3, file = paste0(OUT_DIR, "cellID_obs_S5_Tm3.csv"), row.names = FALSE)

cells_S6_Tm3 <- Cells(S6_Tm3_obj)
cells_S6_Tm3 <- gsub("capture8_", "capture8_output:", cells_S6_Tm3)
cells_S6_Tm3 <- gsub("-1", "x", cells_S6_Tm3)
write.csv(cells_S6_Tm3, file = paste0(OUT_DIR, "cellID_obs_S6_Tm3.csv"), row.names = FALSE)

cells_S7_Tm3 <- Cells(S7_Tm3_obj)
cells_S7_Tm3 <- gsub("capture8_", "capture8_output:", cells_S7_Tm3)
cells_S7_Tm3 <- gsub("-1", "x", cells_S7_Tm3)
write.csv(cells_S7_Tm3, file = paste0(OUT_DIR, "cellID_obs_S7_Tm3.csv"), row.names = FALSE)

embed_S5_Tm3 <- Embeddings(S5_Tm3_obj, reduction = "umap")
embed_names_S5_Tm3 <- row.names(embed_S5_Tm3)
embed_names_S5_Tm3 <- gsub("capture8_", "capture8_output:", embed_names_S5_Tm3)
embed_names_S5_Tm3 <- gsub("-1", "x", embed_names_S5_Tm3)
row.names(embed_S5_Tm3) <- embed_names_S5_Tm3
write.csv(embed_S5_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S5_Tm3.csv"))

embed_S6_Tm3 <- Embeddings(S6_Tm3_obj, reduction = "umap")
embed_names_S6_Tm3 <- row.names(embed_S6_Tm3)
embed_names_S6_Tm3 <- gsub("capture8_", "capture8_output:", embed_names_S6_Tm3)
embed_names_S6_Tm3 <- gsub("-1", "x", embed_names_S6_Tm3)
row.names(embed_S6_Tm3) <- embed_names_S6_Tm3
write.csv(embed_S6_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S6_Tm3.csv"))

embed_S7_Tm3 <- Embeddings(S7_Tm3_obj, reduction = "umap")
embed_names_S7_Tm3 <- row.names(embed_S7_Tm3)
embed_names_S7_Tm3 <- gsub("capture8_", "capture8_output:", embed_names_S7_Tm3)
embed_names_S7_Tm3 <- gsub("-1", "x", embed_names_S7_Tm3)
row.names(embed_S7_Tm3) <- embed_names_S7_Tm3
write.csv(embed_S7_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S7_Tm3.csv"))


clusters_S5_Tm3 <- data.frame(S5_Tm3_obj@meta.data$clust_names)
row.names(clusters_S5_Tm3) <- embed_names_S5_Tm3
write.csv(clusters_S5_Tm3, file = paste0(OUT_DIR, "clusters_S5_Tm3.csv"))

clusters_S6_Tm3 <- data.frame(S6_Tm3_obj@meta.data$clust_names)
row.names(clusters_S6_Tm3) <- embed_names_S6_Tm3
write.csv(clusters_S6_Tm3, file = paste0(OUT_DIR, "clusters_S6_Tm3.csv"))

clusters_S7_Tm3 <- data.frame(S7_Tm3_obj@meta.data$clust_names)
row.names(clusters_S7_Tm3) <- embed_names_S7_Tm3
write.csv(clusters_S7_Tm3, file = paste0(OUT_DIR, "clusters_S7_Tm3.csv"))


# Capture 9 ---------------------------------------------------------------

IN_DIR <- "/capture9_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture9_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S10_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S10-CD34neg_S9_R1_001_Td0_BCG"
)

S1_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S1-CD34neg_S1_R1_001_Td0_BCG"
)

S11_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S11-CD34neg_S10_R1_001_Td0_CTL"
)


cells_S10_Td0 <- Cells(S10_Td0_obj)
cells_S10_Td0 <- gsub("capture9_", "capture9_output:", cells_S10_Td0)
cells_S10_Td0 <- gsub("-1", "x", cells_S10_Td0)
write.csv(cells_S10_Td0, file = paste0(OUT_DIR, "cellID_obs_S10_Td0.csv"), row.names = FALSE)

cells_S1_Td0 <- Cells(S1_Td0_obj)
cells_S1_Td0 <- gsub("capture9_", "capture9_output:", cells_S1_Td0)
cells_S1_Td0 <- gsub("-1", "x", cells_S1_Td0)
write.csv(cells_S1_Td0, file = paste0(OUT_DIR, "cellID_obs_S1_Td0.csv"), row.names = FALSE)

cells_S11_Td0 <- Cells(S11_Td0_obj)
cells_S11_Td0 <- gsub("capture9_", "capture9_output:", cells_S11_Td0)
cells_S11_Td0 <- gsub("-1", "x", cells_S11_Td0)
write.csv(cells_S11_Td0, file = paste0(OUT_DIR, "cellID_obs_S11_Td0.csv"), row.names = FALSE)

embed_S10_Td0 <- Embeddings(S10_Td0_obj, reduction = "umap")
embed_names_S10_Td0 <- row.names(embed_S10_Td0)
embed_names_S10_Td0 <- gsub("capture9_", "capture9_output:", embed_names_S10_Td0)
embed_names_S10_Td0 <- gsub("-1", "x", embed_names_S10_Td0)
row.names(embed_S10_Td0) <- embed_names_S10_Td0
write.csv(embed_S10_Td0, file = paste0(OUT_DIR, "cell_embeddings_S10_Td0.csv"))

embed_S1_Td0 <- Embeddings(S1_Td0_obj, reduction = "umap")
embed_names_S1_Td0 <- row.names(embed_S1_Td0)
embed_names_S1_Td0 <- gsub("capture9_", "capture9_output:", embed_names_S1_Td0)
embed_names_S1_Td0 <- gsub("-1", "x", embed_names_S1_Td0)
row.names(embed_S1_Td0) <- embed_names_S1_Td0
write.csv(embed_S1_Td0, file = paste0(OUT_DIR, "cell_embeddings_S1_Td0.csv"))

embed_S11_Td0 <- Embeddings(S11_Td0_obj, reduction = "umap")
embed_names_S11_Td0 <- row.names(embed_S11_Td0)
embed_names_S11_Td0 <- gsub("capture9_", "capture9_output:", embed_names_S11_Td0)
embed_names_S11_Td0 <- gsub("-1", "x", embed_names_S11_Td0)
row.names(embed_S11_Td0) <- embed_names_S11_Td0
write.csv(embed_S11_Td0, file = paste0(OUT_DIR, "cell_embeddings_S11_Td0.csv"))


clusters_S10_Td0 <- data.frame(S10_Td0_obj@meta.data$clust_names)
row.names(clusters_S10_Td0) <- embed_names_S10_Td0
write.csv(clusters_S10_Td0, file = paste0(OUT_DIR, "clusters_S10_Td0.csv"))

clusters_S1_Td0 <- data.frame(S1_Td0_obj@meta.data$clust_names)
row.names(clusters_S1_Td0) <- embed_names_S1_Td0
write.csv(clusters_S1_Td0, file = paste0(OUT_DIR, "clusters_S1_Td0.csv"))

clusters_S11_Td0 <- data.frame(S11_Td0_obj@meta.data$clust_names)
row.names(clusters_S11_Td0) <- embed_names_S11_Td0
write.csv(clusters_S11_Td0, file = paste0(OUT_DIR, "clusters_S11_Td0.csv"))


# Capture 10 ---------------------------------------------------------------
IN_DIR <- "/capture10_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture10_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S10_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S10-CD34neg_S9_R1_001_Tm3_BCG"
)

S1_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S1-CD34neg_S1_R1_001_Tm3_BCG"
)

S11_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S11-CD34neg_S10_R1_001_Tm3_CTL"
)


cells_S10_Tm3 <- Cells(S10_Tm3_obj)
cells_S10_Tm3 <- gsub("capture10_", "capture10_output:", cells_S10_Tm3)
cells_S10_Tm3 <- gsub("-1", "x", cells_S10_Tm3)
write.csv(cells_S10_Tm3, file = paste0(OUT_DIR, "cellID_obs_S10_Tm3.csv"), row.names = FALSE)

cells_S1_Tm3 <- Cells(S1_Tm3_obj)
cells_S1_Tm3 <- gsub("capture10_", "capture10_output:", cells_S1_Tm3)
cells_S1_Tm3 <- gsub("-1", "x", cells_S1_Tm3)
write.csv(cells_S1_Tm3, file = paste0(OUT_DIR, "cellID_obs_S1_Tm3.csv"), row.names = FALSE)

cells_S11_Tm3 <- Cells(S11_Tm3_obj)
cells_S11_Tm3 <- gsub("capture10_", "capture10_output:", cells_S11_Tm3)
cells_S11_Tm3 <- gsub("-1", "x", cells_S11_Tm3)
write.csv(cells_S11_Tm3, file = paste0(OUT_DIR, "cellID_obs_S11_Tm3.csv"), row.names = FALSE)

embed_S10_Tm3 <- Embeddings(S10_Tm3_obj, reduction = "umap")
embed_names_S10_Tm3 <- row.names(embed_S10_Tm3)
embed_names_S10_Tm3 <- gsub("capture10_", "capture10_output:", embed_names_S10_Tm3)
embed_names_S10_Tm3 <- gsub("-1", "x", embed_names_S10_Tm3)
row.names(embed_S10_Tm3) <- embed_names_S10_Tm3
write.csv(embed_S10_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S10_Tm3.csv"))

embed_S1_Tm3 <- Embeddings(S1_Tm3_obj, reduction = "umap")
embed_names_S1_Tm3 <- row.names(embed_S1_Tm3)
embed_names_S1_Tm3 <- gsub("capture10_", "capture10_output:", embed_names_S1_Tm3)
embed_names_S1_Tm3 <- gsub("-1", "x", embed_names_S1_Tm3)
row.names(embed_S1_Tm3) <- embed_names_S1_Tm3
write.csv(embed_S1_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S1_Tm3.csv"))

embed_S11_Tm3 <- Embeddings(S11_Tm3_obj, reduction = "umap")
embed_names_S11_Tm3 <- row.names(embed_S11_Tm3)
embed_names_S11_Tm3 <- gsub("capture10_", "capture10_output:", embed_names_S11_Tm3)
embed_names_S11_Tm3 <- gsub("-1", "x", embed_names_S11_Tm3)
row.names(embed_S11_Tm3) <- embed_names_S11_Tm3
write.csv(embed_S11_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S11_Tm3.csv"))


clusters_S10_Tm3 <- data.frame(S10_Tm3_obj@meta.data$clust_names)
row.names(clusters_S10_Tm3) <- embed_names_S10_Tm3
write.csv(clusters_S10_Tm3, file = paste0(OUT_DIR, "clusters_S10_Tm3.csv"))

clusters_S1_Tm3 <- data.frame(S1_Tm3_obj@meta.data$clust_names)
row.names(clusters_S1_Tm3) <- embed_names_S1_Tm3
write.csv(clusters_S1_Tm3, file = paste0(OUT_DIR, "clusters_S1_Tm3.csv"))

clusters_S11_Tm3 <- data.frame(S11_Tm3_obj@meta.data$clust_names)
row.names(clusters_S11_Tm3) <- embed_names_S11_Tm3
write.csv(clusters_S11_Tm3, file = paste0(OUT_DIR, "clusters_S11_Tm3.csv"))


# Capture 11 ---------------------------------------------------------------

IN_DIR <- "/capture11_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture11_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S12_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S12-CD34neg_S11_R1_001_Td0_BCG"
)

S2_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S2-CD34neg_S2_R1_001_Td0_BCG"
)

S3_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S3-CD34neg_S3_R1_001_Td0_CTL"
)


cells_S12_Td0 <- Cells(S12_Td0_obj)
cells_S12_Td0 <- gsub("capture11_", "capture11_output:", cells_S12_Td0)
cells_S12_Td0 <- gsub("-1", "x", cells_S12_Td0)
write.csv(cells_S12_Td0, file = paste0(OUT_DIR, "cellID_obs_S12_Td0.csv"), row.names = FALSE)

cells_S2_Td0 <- Cells(S2_Td0_obj)
cells_S2_Td0 <- gsub("capture11_", "capture11_output:", cells_S2_Td0)
cells_S2_Td0 <- gsub("-1", "x", cells_S2_Td0)
write.csv(cells_S2_Td0, file = paste0(OUT_DIR, "cellID_obs_S2_Td0.csv"), row.names = FALSE)

cells_S3_Td0 <- Cells(S3_Td0_obj)
cells_S3_Td0 <- gsub("capture11_", "capture11_output:", cells_S3_Td0)
cells_S3_Td0 <- gsub("-1", "x", cells_S3_Td0)
write.csv(cells_S3_Td0, file = paste0(OUT_DIR, "cellID_obs_S3_Td0.csv"), row.names = FALSE)

embed_S12_Td0 <- Embeddings(S12_Td0_obj, reduction = "umap")
embed_names_S12_Td0 <- row.names(embed_S12_Td0)
embed_names_S12_Td0 <- gsub("capture11_", "capture11_output:", embed_names_S12_Td0)
embed_names_S12_Td0 <- gsub("-1", "x", embed_names_S12_Td0)
row.names(embed_S12_Td0) <- embed_names_S12_Td0
write.csv(embed_S12_Td0, file = paste0(OUT_DIR, "cell_embeddings_S12_Td0.csv"))

embed_S2_Td0 <- Embeddings(S2_Td0_obj, reduction = "umap")
embed_names_S2_Td0 <- row.names(embed_S2_Td0)
embed_names_S2_Td0 <- gsub("capture11_", "capture11_output:", embed_names_S2_Td0)
embed_names_S2_Td0 <- gsub("-1", "x", embed_names_S2_Td0)
row.names(embed_S2_Td0) <- embed_names_S2_Td0
write.csv(embed_S2_Td0, file = paste0(OUT_DIR, "cell_embeddings_S2_Td0.csv"))

embed_S3_Td0 <- Embeddings(S3_Td0_obj, reduction = "umap")
embed_names_S3_Td0 <- row.names(embed_S3_Td0)
embed_names_S3_Td0 <- gsub("capture11_", "capture11_output:", embed_names_S3_Td0)
embed_names_S3_Td0 <- gsub("-1", "x", embed_names_S3_Td0)
row.names(embed_S3_Td0) <- embed_names_S3_Td0
write.csv(embed_S3_Td0, file = paste0(OUT_DIR, "cell_embeddings_S3_Td0.csv"))


clusters_S12_Td0 <- data.frame(S12_Td0_obj@meta.data$clust_names)
row.names(clusters_S12_Td0) <- embed_names_S12_Td0
write.csv(clusters_S12_Td0, file = paste0(OUT_DIR, "clusters_S12_Td0.csv"))

clusters_S2_Td0 <- data.frame(S2_Td0_obj@meta.data$clust_names)
row.names(clusters_S2_Td0) <- embed_names_S2_Td0
write.csv(clusters_S2_Td0, file = paste0(OUT_DIR, "clusters_S2_Td0.csv"))

clusters_S3_Td0 <- data.frame(S3_Td0_obj@meta.data$clust_names)
row.names(clusters_S3_Td0) <- embed_names_S3_Td0
write.csv(clusters_S3_Td0, file = paste0(OUT_DIR, "clusters_S3_Td0.csv"))



# Capture 12 ---------------------------------------------------------------

IN_DIR <- "/capture12_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture12_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S12_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S12-CD34neg_S11_R1_001_Tm3_BCG"
)

S2_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S2-CD34neg_S2_R1_001_Tm3_BCG"
)

S3_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S3-CD34neg_S3_R1_001_Tm3_CTL"
)


cells_S12_Tm3 <- Cells(S12_Tm3_obj)
cells_S12_Tm3 <- gsub("capture12_", "capture12_output:", cells_S12_Tm3)
cells_S12_Tm3 <- gsub("-1", "x", cells_S12_Tm3)
write.csv(cells_S12_Tm3, file = paste0(OUT_DIR, "cellID_obs_S12_Tm3.csv"), row.names = FALSE)

cells_S2_Tm3 <- Cells(S2_Tm3_obj)
cells_S2_Tm3 <- gsub("capture12_", "capture12_output:", cells_S2_Tm3)
cells_S2_Tm3 <- gsub("-1", "x", cells_S2_Tm3)
write.csv(cells_S2_Tm3, file = paste0(OUT_DIR, "cellID_obs_S2_Tm3.csv"), row.names = FALSE)

cells_S3_Tm3 <- Cells(S3_Tm3_obj)
cells_S3_Tm3 <- gsub("capture12_", "capture12_output:", cells_S3_Tm3)
cells_S3_Tm3 <- gsub("-1", "x", cells_S3_Tm3)
write.csv(cells_S3_Tm3, file = paste0(OUT_DIR, "cellID_obs_S3_Tm3.csv"), row.names = FALSE)

embed_S12_Tm3 <- Embeddings(S12_Tm3_obj, reduction = "umap")
embed_names_S12_Tm3 <- row.names(embed_S12_Tm3)
embed_names_S12_Tm3 <- gsub("capture12_", "capture12_output:", embed_names_S12_Tm3)
embed_names_S12_Tm3 <- gsub("-1", "x", embed_names_S12_Tm3)
row.names(embed_S12_Tm3) <- embed_names_S12_Tm3
write.csv(embed_S12_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S12_Tm3.csv"))

embed_S2_Tm3 <- Embeddings(S2_Tm3_obj, reduction = "umap")
embed_names_S2_Tm3 <- row.names(embed_S2_Tm3)
embed_names_S2_Tm3 <- gsub("capture12_", "capture12_output:", embed_names_S2_Tm3)
embed_names_S2_Tm3 <- gsub("-1", "x", embed_names_S2_Tm3)
row.names(embed_S2_Tm3) <- embed_names_S2_Tm3
write.csv(embed_S2_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S2_Tm3.csv"))

embed_S3_Tm3 <- Embeddings(S3_Tm3_obj, reduction = "umap")
embed_names_S3_Tm3 <- row.names(embed_S3_Tm3)
embed_names_S3_Tm3 <- gsub("capture12_", "capture12_output:", embed_names_S3_Tm3)
embed_names_S3_Tm3 <- gsub("-1", "x", embed_names_S3_Tm3)
row.names(embed_S3_Tm3) <- embed_names_S3_Tm3
write.csv(embed_S3_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S3_Tm3.csv"))


clusters_S12_Tm3 <- data.frame(S12_Tm3_obj@meta.data$clust_names)
row.names(clusters_S12_Tm3) <- embed_names_S12_Tm3
write.csv(clusters_S12_Tm3, file = paste0(OUT_DIR, "clusters_S12_Tm3.csv"))

clusters_S2_Tm3 <- data.frame(S2_Tm3_obj@meta.data$clust_names)
row.names(clusters_S2_Tm3) <- embed_names_S2_Tm3
write.csv(clusters_S2_Tm3, file = paste0(OUT_DIR, "clusters_S2_Tm3.csv"))

clusters_S3_Tm3 <- data.frame(S3_Tm3_obj@meta.data$clust_names)
row.names(clusters_S3_Tm3) <- embed_names_S3_Tm3
write.csv(clusters_S3_Tm3, file = paste0(OUT_DIR, "clusters_S3_Tm3.csv"))



# Capture 13 ---------------------------------------------------------------

IN_DIR <- "/capture13_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture13_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S14_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S14-CD34neg_S13_R1_001_Td0_BCG"
)

S17_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S17-CD34neg_S16_R1_001_Td0_BCG"
)

S18_Td0_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S18-CD34neg_S17_R1_001_Td0_CTL"
)


cells_S14_Td0 <- Cells(S14_Td0_obj)
cells_S14_Td0 <- gsub("capture13_", "capture13_output:", cells_S14_Td0)
cells_S14_Td0 <- gsub("-1", "x", cells_S14_Td0)
write.csv(cells_S14_Td0, file = paste0(OUT_DIR, "cellID_obs_S14_Td0.csv"), row.names = FALSE)

cells_S17_Td0 <- Cells(S17_Td0_obj)
cells_S17_Td0 <- gsub("capture13_", "capture13_output:", cells_S17_Td0)
cells_S17_Td0 <- gsub("-1", "x", cells_S17_Td0)
write.csv(cells_S17_Td0, file = paste0(OUT_DIR, "cellID_obs_S17_Td0.csv"), row.names = FALSE)

cells_S18_Td0 <- Cells(S18_Td0_obj)
cells_S18_Td0 <- gsub("capture13_", "capture13_output:", cells_S18_Td0)
cells_S18_Td0 <- gsub("-1", "x", cells_S18_Td0)
write.csv(cells_S18_Td0, file = paste0(OUT_DIR, "cellID_obs_S18_Td0.csv"), row.names = FALSE)

embed_S14_Td0 <- Embeddings(S14_Td0_obj, reduction = "umap")
embed_names_S14_Td0 <- row.names(embed_S14_Td0)
embed_names_S14_Td0 <- gsub("capture13_", "capture13_output:", embed_names_S14_Td0)
embed_names_S14_Td0 <- gsub("-1", "x", embed_names_S14_Td0)
row.names(embed_S14_Td0) <- embed_names_S14_Td0
write.csv(embed_S14_Td0, file = paste0(OUT_DIR, "cell_embeddings_S14_Td0.csv"))

embed_S17_Td0 <- Embeddings(S17_Td0_obj, reduction = "umap")
embed_names_S17_Td0 <- row.names(embed_S17_Td0)
embed_names_S17_Td0 <- gsub("capture13_", "capture13_output:", embed_names_S17_Td0)
embed_names_S17_Td0 <- gsub("-1", "x", embed_names_S17_Td0)
row.names(embed_S17_Td0) <- embed_names_S17_Td0
write.csv(embed_S17_Td0, file = paste0(OUT_DIR, "cell_embeddings_S17_Td0.csv"))

embed_S18_Td0 <- Embeddings(S18_Td0_obj, reduction = "umap")
embed_names_S18_Td0 <- row.names(embed_S18_Td0)
embed_names_S18_Td0 <- gsub("capture13_", "capture13_output:", embed_names_S18_Td0)
embed_names_S18_Td0 <- gsub("-1", "x", embed_names_S18_Td0)
row.names(embed_S18_Td0) <- embed_names_S18_Td0
write.csv(embed_S18_Td0, file = paste0(OUT_DIR, "cell_embeddings_S18_Td0.csv"))


clusters_S14_Td0 <- data.frame(S14_Td0_obj@meta.data$clust_names)
row.names(clusters_S14_Td0) <- embed_names_S14_Td0
write.csv(clusters_S14_Td0, file = paste0(OUT_DIR, "clusters_S14_Td0.csv"))

clusters_S17_Td0 <- data.frame(S17_Td0_obj@meta.data$clust_names)
row.names(clusters_S17_Td0) <- embed_names_S17_Td0
write.csv(clusters_S17_Td0, file = paste0(OUT_DIR, "clusters_S17_Td0.csv"))

clusters_S18_Td0 <- data.frame(S18_Td0_obj@meta.data$clust_names)
row.names(clusters_S18_Td0) <- embed_names_S18_Td0
write.csv(clusters_S18_Td0, file = paste0(OUT_DIR, "clusters_S18_Td0.csv"))


# Capture 14 ---------------------------------------------------------------

IN_DIR <- "/capture14_output/velocyto/"
OUT_DIR <- "wVelocity_input/"

#Read in data
ldat <- read.loom.matrices(paste0(IN_DIR, "capture14_output.loom"))
emat <- ldat$spliced
nmat <- ldat$unspliced

seurat_obj <- readRDS(file = "../label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds")
DefaultAssay(seurat_obj) <- "RNA"


S14_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S14-CD34neg_S13_R1_001_Tm3_BCG"
)

S17_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S17-CD34neg_S16_R1_001_Tm3_BCG"
)

S18_Tm3_obj <- subset(
  x = seurat_obj,
  subset = UNIQUE_ID %in% "SNG-LB-SS-1S-RS-S18-CD34neg_S17_R1_001_Tm3_CTL"
)


cells_S14_Tm3 <- Cells(S14_Tm3_obj)
cells_S14_Tm3 <- gsub("capture14_", "capture14_output:", cells_S14_Tm3)
cells_S14_Tm3 <- gsub("-1", "x", cells_S14_Tm3)
write.csv(cells_S14_Tm3, file = paste0(OUT_DIR, "cellID_obs_S14_Tm3.csv"), row.names = FALSE)

cells_S17_Tm3 <- Cells(S17_Tm3_obj)
cells_S17_Tm3 <- gsub("capture14_", "capture14_output:", cells_S17_Tm3)
cells_S17_Tm3 <- gsub("-1", "x", cells_S17_Tm3)
write.csv(cells_S17_Tm3, file = paste0(OUT_DIR, "cellID_obs_S17_Tm3.csv"), row.names = FALSE)

cells_S18_Tm3 <- Cells(S18_Tm3_obj)
cells_S18_Tm3 <- gsub("capture14_", "capture14_output:", cells_S18_Tm3)
cells_S18_Tm3 <- gsub("-1", "x", cells_S18_Tm3)
write.csv(cells_S18_Tm3, file = paste0(OUT_DIR, "cellID_obs_S18_Tm3.csv"), row.names = FALSE)

embed_S14_Tm3 <- Embeddings(S14_Tm3_obj, reduction = "umap")
embed_names_S14_Tm3 <- row.names(embed_S14_Tm3)
embed_names_S14_Tm3 <- gsub("capture14_", "capture14_output:", embed_names_S14_Tm3)
embed_names_S14_Tm3 <- gsub("-1", "x", embed_names_S14_Tm3)
row.names(embed_S14_Tm3) <- embed_names_S14_Tm3
write.csv(embed_S14_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S14_Tm3.csv"))

embed_S17_Tm3 <- Embeddings(S17_Tm3_obj, reduction = "umap")
embed_names_S17_Tm3 <- row.names(embed_S17_Tm3)
embed_names_S17_Tm3 <- gsub("capture14_", "capture14_output:", embed_names_S17_Tm3)
embed_names_S17_Tm3 <- gsub("-1", "x", embed_names_S17_Tm3)
row.names(embed_S17_Tm3) <- embed_names_S17_Tm3
write.csv(embed_S17_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S17_Tm3.csv"))

embed_S18_Tm3 <- Embeddings(S18_Tm3_obj, reduction = "umap")
embed_names_S18_Tm3 <- row.names(embed_S18_Tm3)
embed_names_S18_Tm3 <- gsub("capture14_", "capture14_output:", embed_names_S18_Tm3)
embed_names_S18_Tm3 <- gsub("-1", "x", embed_names_S18_Tm3)
row.names(embed_S18_Tm3) <- embed_names_S18_Tm3
write.csv(embed_S18_Tm3, file = paste0(OUT_DIR, "cell_embeddings_S18_Tm3.csv"))


clusters_S14_Tm3 <- data.frame(S14_Tm3_obj@meta.data$clust_names)
row.names(clusters_S14_Tm3) <- embed_names_S14_Tm3
write.csv(clusters_S14_Tm3, file = paste0(OUT_DIR, "clusters_S14_Tm3.csv"))

clusters_S17_Tm3 <- data.frame(S17_Tm3_obj@meta.data$clust_names)
row.names(clusters_S17_Tm3) <- embed_names_S17_Tm3
write.csv(clusters_S17_Tm3, file = paste0(OUT_DIR, "clusters_S17_Tm3.csv"))

clusters_S18_Tm3 <- data.frame(S18_Tm3_obj@meta.data$clust_names)
row.names(clusters_S18_Tm3) <- embed_names_S18_Tm3
write.csv(clusters_S18_Tm3, file = paste0(OUT_DIR, "clusters_S18_Tm3.csv"))
