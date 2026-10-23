library(ArchR)
addArchRThreads(threads = 16) 
addArchRGenome("hg38")
library(readr)
library(harmony)
library(Seurat)
library(stringr)
library(textTinyR)
library(pbapply)

# setup ----------------------------------------------------------------------------------------
setwd("/scATAC_analyses/2_main_analyses")
OUT_DIR <- "pseudobulk_clusters_harmony_all_combined/"
dir.create(OUT_DIR)

projBM <- loadArchRProject(path = "Unedited_with_Clusters_harmony_peaks_combined", force = FALSE, showLogo = TRUE)
getAvailableMatrices(projBM)


# Get peak matrix ---------------------------------------------------------

peakmat <- getMatrixFromProject(
  ArchRProj = projBM,
  useMatrix = "PeakMatrix",
  useSeqnames = NULL,
  verbose = TRUE,
  binarize = FALSE,
  threads = getArchRThreads(),
  logFile = createLogFile("getMatrixFromProject")
)

sparse_mat <- peakmat@assays@data@listData$PeakMatrix
sparse_dims <- dim(sparse_mat)
rownames(sparse_mat) <- paste0(1:sparse_dims[1])

df <- getCellColData(projBM)


# Split by clusters_harmony ----------------------------------------------
celltypes <- unique(df$Clusters_harmony)
for(i in 1:length(celltypes)){
  name <- celltypes[i]
  idxSample <- BiocGenerics::which(projBM$Clusters_harmony %in% paste0(name))
  cellsSample <- projBM$cellNames[idxSample]
  proj <- projBM[cellsSample, ]
  samples <- unique(proj$cellNames)
  cluster_cols <- which(colnames(sparse_mat) %in% samples)
  sparse_subset <- sparse_mat[,cluster_cols]
  
  for(k in 1:length(colnames(sparse_subset)))
  {
    name_original <- colnames(sparse_subset)[k]
    colnames(sparse_subset)[k] <- subset(df, row.names(df)==name_original)$Sample
  }
  IDs <- colnames(sparse_subset)
  unique_ID_list <- as.list(unique(IDs))
  
  pseudobulk <- as.data.frame(pblapply(unique_ID_list, FUN = function(x){sparse_Sums(sparse_subset[,IDs == x, drop = FALSE], rowSums = TRUE)}))
  cellcount <- pblapply(unique_ID_list, FUN = function(x){ncol(sparse_subset[,IDs == x, drop = FALSE])})
  
  colnames(pseudobulk) <- names(cellcount) <- unique_ID_list
  rownames(pseudobulk) <- rownames(sparse_mat)
  saveRDS(pseudobulk, file = paste0(OUT_DIR,name,"_pseudobulk.rds"))
  
  frags_in_peaks <- colSums(pseudobulk)
  
  stats_df <- data.frame(
    peak_frags <- frags_in_peaks,
    cellcount <- unlist(cellcount)
  )
  colnames(stats_df) <- c("peak_frags", "cellcount")
  assign(paste0(name, "_stats_df"), stats_df)
  saveRDS(stats_df, file=paste0(OUT_DIR, name, "_stats_df.rds"))
  print(i)
  print(name)
  print(length(cellsSample))
  
}
