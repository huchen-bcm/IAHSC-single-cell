#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

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
library(stringr)
library(ggpubr)

# In this script we make a UMAP of the HSC subclusters and plot various variables for QC
# We also calculate an MPP score for each HSC subcluster based on the average expression of "THY1", "PTPRC", and "ITGA6"

# setup -------------------------------------------------------------------
options(future.globals.maxSize = 1000 * 1024^2)
setwd("/scRNA_analyses/2_main_analyses/HSC_subcluster_analyses")
IN_DIR <- "0_SCT_cluster/"
OUT_DIR <- "1_UMAP_label_clusters/"
dir.create(OUT_DIR)


# Read in merged_scaled_integrated_clustered HSC object -----------------------
allCells_integrated <- readRDS(file = paste0(IN_DIR, "HSC_obj_SCT_integrated_withClusters.rds"))


# print UMAP --------------------------------------------------------------
pdf(paste0(OUT_DIR,"UMAP.pdf"), width=7, height=5)
print(DimPlot(allCells_integrated, reduction = "umap", label = TRUE) + NoLegend())
dev.off()


# Plot categorical covariates ---------------------------------------------
covariates <- c("vaccination", "timepoint", "batchID")
for(i in 1:length(covariates))
{
  covar <- covariates[i]
  pdf(paste0(OUT_DIR,"UMAP_",covar,".pdf"), width=7, height=5)
  print(DimPlot(allCells_integrated, reduction = "umap", group.by = covar))
  dev.off()
}

# Plot numeric covariates -------------------------------------------------
covariates <- c("nCount_RNA", "nFeature_RNA", "percent.mt")
for(i in 1:length(covariates)){
  covar <- covariates[i]
  pdf(paste0(OUT_DIR,"UMAP_",covar,".pdf"), width=7, height=5)
  print(FeaturePlot(allCells_integrated, reduction = "umap", features = covar))
  dev.off()
}

pdf(paste0(OUT_DIR,"UMAP_split_CTL_BCG.pdf"), width = 10, height = 4)
DimPlot(allCells_integrated, reduction = "umap", split.by = "vaccination")
dev.off()

pdf(paste0(OUT_DIR,"UMAP_split_timepoint.pdf"), width = 10, height = 4)
DimPlot(allCells_integrated, reduction = "umap", split.by = "timepoint")
dev.off()

png(paste0(OUT_DIR,"single_cell_stats.png"), width = 2160, height = 720)
VlnPlot(allCells_integrated, features=c("nCount_RNA","nFeature_RNA","percent.mt"))
dev.off()


# Label cells with UNIQUE sample label combining donor+timepoint+vaccination ------------------------------------
UNIQUE_ID <- vector(length=length(row.names(allCells_integrated@meta.data)))

for (i in 1:length(row.names(allCells_integrated@meta.data)))
{
  UNIQUE_ID[i] <- paste0(allCells_integrated$BEST[i],"_",allCells_integrated$timepoint[i],"_",allCells_integrated$vaccination[i])
}

allCells_integrated$UNIQUE_ID <- UNIQUE_ID
saveRDS(allCells_integrated, file=paste0(OUT_DIR,"HSC_final_UNIQUE_IDS.rds"))


# Find markers for every cluster ------------------------------------------
HSC.markers <- FindAllMarkers(allCells_integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

HSC_markers_data <- data.frame(HSC.markers %>%
  group_by(cluster) %>%
  slice_max(n = 20, order_by = avg_log2FC))

write.table(HSC_markers_data, file=paste0(OUT_DIR, "HSC_markers20.txt"))

DefaultAssay(allCells_integrated) <- "SCT"
FeaturePlot(allCells_integrated, features = c("THY1"), reduction = "umap", cols=c('lightgrey', '#BA68C8')) #CD90
FeaturePlot(allCells_integrated, features = c("PTPRC"), reduction = "umap", cols=c('lightgrey', '#BA68C8')) #CD45RA
FeaturePlot(allCells_integrated, features = c("ITGA6"), reduction = "umap", cols=c('lightgrey', '#BA68C8')) #CD49f


######################
###   MPP score    ###
######################

# look for MPP expression signatures --------------------------------------
allCells_integrated <- readRDS(file=paste0(OUT_DIR,"HSC_final_UNIQUE_IDS.rds"))

MPP_average_expression <- data.frame(AverageExpression(allCells_integrated, assays='SCT',features = c("THY1", "PTPRC", "ITGA6"), group.by="seurat_clusters")[[1]])

MPP_scaled_data <- data.frame(scale(t(MPP_average_expression), scale=TRUE, center=TRUE))

MPP_scaled_data$average <- rowMeans(MPP_scaled_data)
row.names(MPP_scaled_data) <- str_replace(row.names(MPP_scaled_data),"X","")


# plot linegraph of average MPP score per subcluster ----------------------------------------- 
MPP_scaled_data_ordered <- MPP_scaled_data[order(MPP_scaled_data$average, decreasing = TRUE),]

df <- data.frame(
  cluster <- paste0("c",row.names(MPP_scaled_data_ordered)),
  MPP_score <- MPP_scaled_data_ordered$average
)
colnames(df) <- c("cluster", "MPP_score")
df$cluster <- factor(df$cluster, c(paste0("c",row.names(MPP_scaled_data_ordered))))

myColorScale <- c('c0'='#EF5350','c1'='#CC0C00','c2'='#FF7F0E','c3'='#FFD147','c4'='#CC9900',
                  'c5'='#33CC00', 'c6' ='#99CC00', 'c7'='#8C5782', 'c8'='#990080', 'c9'='#BA68C8')

tiff(paste0(OUT_DIR, "MPP_score_linegraph.tiff"), units="in", width=3.5, height=3.1, res=300)
ggplot(df, aes(x=cluster, y=MPP_score)) + geom_point(shape=21,size=4.5, color="grey", aes(fill=cluster)) + theme_test() + scale_fill_manual(values=myColorScale)+
  theme(legend.position = "none") + labs(x="", y="CD90/CD49f/CD45RA Score")
dev.off()


# GGplot version of the UMAP ----------------------------------------------
embed_query_BM_new <- data.frame(Embeddings(allCells_integrated, reduction = "umap"))
#check that cells are in the same order
which(row.names(embed_query_BM_new) != row.names(allCells_integrated@meta.data) )
df <- data.frame(allCells_integrated@meta.data)

embed_query_BM_new$cluster <- as.character(df$seurat_clusters)

myColorScale <- c('0'='#EF5350','1'='#CC0C00','2'='#FF7F0E','3'='#FFD147','4'='#CC9900',
                  '5'='#33CC00', '6' ='#99CC00', '7'='#8C5782', '8'='#990080', '9c'='#BA68C8')

tiff(paste0(OUT_DIR, "subclusters_newUMAP_ggplot.tiff"), units="in", width=3.5, height=3.4, res=300)
ggplot(embed_query_BM_new, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + geom_point(aes(color=cluster), alpha=0.6) + theme_classic() + theme(legend.position="none")+
  scale_color_manual(name= "", values=myColorScale) + labs(title="HSCs sub-clustered")
dev.off()



# MPP score on new UMAP ---------------------------------------------------

embed_query_BM <- data.frame(Embeddings(allCells_integrated, reduction = "umap"))
#check that cells are in the same order
which(row.names(embed_query_BM) != row.names(allCells_integrated@meta.data) )
df <- data.frame(allCells_integrated@meta.data)

MPP_average_score <- vector(length=length(row.names(df)))
for(i in 1:length(row.names(MPP_scaled_data)))
{
  MPP_average_score[which(df$seurat_clusters == row.names(MPP_scaled_data)[i])] <- MPP_scaled_data$average[i]
}
df$MPP_avg_score <- MPP_average_score

embed_query_BM$avg_score <- as.numeric(df$MPP_avg_score)

tiff(paste0(OUT_DIR, "MPP_score_newUMAP_ggplot.tiff"), units="in", width=3.7, height=3, res=300)
ggplot(embed_query_BM, aes(x=UMAP_1, y=UMAP_2, color=avg_score)) + geom_point(aes(color=avg_score), alpha=0.6) + theme_classic() + 
  scale_color_gradient(low="#FFF3E0", high="#E65100") + labs(title="CD90/CD49f/CD45RA Score") + theme(legend.position = "none")
dev.off()



# View the clusters in the original UMAP ----------------------------------

#S9 Td0 needs to be filtered out from this object!
all <- readRDS(file=paste0("../Analysis12_label_transfer_emmreml_edited/label_transfer/allCells_integrated_label_transfer_UNIQUE_IDS.rds"))
DefaultAssay(all) <- "RNA"

#Filter out S9 Td0
unique_ids <- unique(all$UNIQUE_ID)
unique_ids = unique_ids[!(unique_ids %in% "SNG-LB-SS-1S-RS-S9-CD34neg_S8_R1_001_Td0_BCG" )]
all<- subset(
  x = all,
  subset = UNIQUE_ID %in% unique_ids
)


# subset HSCs only ----------------------------------------------

all_meta <- all@meta.data

HSC_barcodes <- row.names(all_meta)[which(all_meta$clust_names %in% c('HSC_a', 'HSC_b'))]

all$temp <- row.names(all@meta.data)

HSC_full_obj <- subset(
  x = all,
  subset = temp %in% HSC_barcodes
)
HSC_full_obj$temp <- NULL


# Add subcluster labels to meta data -------------------------------------
subcluster <- vector(length=length(row.names(HSC_full_obj@meta.data)))

unique(allCells_integrated$seurat_clusters) #0-9
c0_barcodes <- row.names(allCells_integrated@meta.data)[which(allCells_integrated$seurat_clusters == 0)]
subcluster[which(row.names(HSC_full_obj@meta.data) %in% c0_barcodes)] <- 0
c1_barcodes <- row.names(allCells_integrated@meta.data)[which(allCells_integrated$seurat_clusters == 1)]
subcluster[which(row.names(HSC_full_obj@meta.data) %in% c1_barcodes)] <- 1
c2_barcodes <- row.names(allCells_integrated@meta.data)[which(allCells_integrated$seurat_clusters == 2)]
subcluster[which(row.names(HSC_full_obj@meta.data) %in% c2_barcodes)] <- 2
c3_barcodes <- row.names(allCells_integrated@meta.data)[which(allCells_integrated$seurat_clusters == 3)]
subcluster[which(row.names(HSC_full_obj@meta.data) %in% c3_barcodes)] <- 3
c4_barcodes <- row.names(allCells_integrated@meta.data)[which(allCells_integrated$seurat_clusters == 4)]
subcluster[which(row.names(HSC_full_obj@meta.data) %in% c4_barcodes)] <- 4
c5_barcodes <- row.names(allCells_integrated@meta.data)[which(allCells_integrated$seurat_clusters == 5)]
subcluster[which(row.names(HSC_full_obj@meta.data) %in% c5_barcodes)] <- 5
c6_barcodes <- row.names(allCells_integrated@meta.data)[which(allCells_integrated$seurat_clusters == 6)]
subcluster[which(row.names(HSC_full_obj@meta.data) %in% c6_barcodes)] <- 6
c7_barcodes <- row.names(allCells_integrated@meta.data)[which(allCells_integrated$seurat_clusters == 7)]
subcluster[which(row.names(HSC_full_obj@meta.data) %in% c7_barcodes)] <- 7
c8_barcodes <- row.names(allCells_integrated@meta.data)[which(allCells_integrated$seurat_clusters == 8)]
subcluster[which(row.names(HSC_full_obj@meta.data) %in% c8_barcodes)] <- 8
c9_barcodes <- row.names(allCells_integrated@meta.data)[which(allCells_integrated$seurat_clusters == 9)]
subcluster[which(row.names(HSC_full_obj@meta.data) %in% c9_barcodes)] <- 9

HSC_full_obj$subcluster <- subcluster

embed_query_BM <- data.frame(Embeddings(HSC_full_obj, reduction = "umap"))
#check that cells are in the same order
which(row.names(embed_query_BM) != row.names(HSC_full_obj@meta.data) )
df <- data.frame(HSC_full_obj@meta.data)

embed_query_BM$cluster <- as.character(df$subcluster)

myColorScale <- c('0'='#EF5350','1'='#CC0C00','2'='#FF7F0E','3'='#FFD147','4'='#CC9900',
                  '5'='#33CC00', '6' ='#99CC00', '7'='#8C5782', '8'='#990080', '9c'='#BA68C8')

tiff(paste0(OUT_DIR, "subclusters_origUMAP_ggplot.tiff"), units="in", width=4, height=3, res=300)
ggplot(embed_query_BM, aes(x=UMAP_1, y=UMAP_2, color=cluster)) + geom_point(aes(color=cluster), alpha=0.6) + theme_classic() + theme(legend.position="none")+
  scale_color_manual(name= "", values=myColorScale)+xlim(-5,5)
dev.off()



# View MPP score on original UMAP -----------------------------------------

embed_query_BM <- data.frame(Embeddings(HSC_full_obj, reduction = "umap"))
#check that cells are in the same order
which(row.names(embed_query_BM) != row.names(HSC_full_obj@meta.data) )
df <- data.frame(HSC_full_obj@meta.data)

MPP_average_score <- vector(length=length(row.names(df)))
for(i in 1:length(row.names(MPP_scaled_data)))
{
  MPP_average_score[which(df$subcluster == row.names(MPP_scaled_data)[i])] <- MPP_scaled_data$average[i]
}
df$MPP_avg_score <- MPP_average_score

embed_query_BM$avg_score <- as.numeric(df$MPP_avg_score)


tiff(paste0(OUT_DIR, "MPP_score_origUMAP_ggplot.tiff"), units="in", width=3.7, height=3, res=300)
ggplot(embed_query_BM, aes(x=UMAP_1, y=UMAP_2, color=avg_score)) + geom_point(aes(color=avg_score), alpha=0.6) + theme_classic() + 
  scale_color_gradient(low="#FFF3E0", high="#E65100") + labs(title="CD90/CD49f/CD45RA Score") + theme(legend.position = "none") + xlim(-5,5)
dev.off()




