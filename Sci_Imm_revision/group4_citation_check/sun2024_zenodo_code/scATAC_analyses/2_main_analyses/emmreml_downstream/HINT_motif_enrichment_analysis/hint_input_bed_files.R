library(ArchR)
addArchRThreads(threads = 16) 
addArchRGenome("hg38")
library(readr)
library(harmony)
library(Seurat)
library(stringr)
library(textTinyR)
library(pbapply)

setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE")
OUT_DIR <- "emmreml_downstream/HINT_motif_enrichment_analysis/hint_input_bed_files/"
dir.create(OUT_DIR)

projBM <- loadArchRProject(path = "../analysis1_Clusters_harmony/Unedited_with_Clusters_harmony_peaks_combined", force = FALSE, showLogo = TRUE)
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


assayNames(peakmat)
peakmat@colData
dim(peakmat)


# Add rownames to peakmat -------------------------------------------------

rowRanges(peakmat)$rownames <- as.character(1:dim(peakmat)[1])

df_peakmat <- as.data.frame(rowRanges(peakmat))

homer_background_peaks <- df_peakmat[,c("seqnames", "start", "end", "rownames")]
colnames(homer_background_peaks) <- c("Chr", "Start", "End", "Geneid")
fwrite(homer_background_peaks, file=paste0(OUT_DIR,"homer_background_peaks.bed"), row.names=FALSE, col.names=FALSE, sep = "\t") 

clusters <- c("C3", "C4", "C5", "C6", "C7","C9", "C10", "C12", "C14", "C15", "C17", "C18",  "C21", "C22","C23", "C24")
clusters_renamed <- c("CMP1", "CMP2", "HSC1", "HSC2", "MEP1","Pre-BNK", "GMP1","CLP", "MEP2", "MEP3", "CMP3", "MEP4", "CMP4", "CMP5", "GMP2", "GMP3")



# Get Significant peaks for CMPs --------------------------------------------------

CMP_clusters <- clusters[which(clusters_renamed %in% c("CMP1", "CMP2", "CMP3", "CMP4", "CMP5"))]
options <- c(1,1,2.5,2.25,1)
CMP_peaks_all <- vector()

for(i in 1:length(CMP_clusters)){
  name <- CMP_clusters[i]
  
  res_full <- read.csv(file=paste0("0X_test_filter_thresholds/res_full/", name, "_res_full_option",options[i] ,".csv"), header=TRUE)
  peaks_sig <- subset(res_full, res_full[,"BH_corrected"] < 0.1)$X
  
  CMP_peaks_all <- append(CMP_peaks_all, peaks_sig)
}

CMP_peaks_unique <- unique(CMP_peaks_all)

peakmat_subset <- subset(peakmat, rowRanges(peakmat)$rownames %in% as.character(CMP_peaks_unique))

df_peakmat_subset <- as.data.frame(rowRanges(peakmat_subset))

homer_peaks <- df_peakmat_subset[,c("seqnames", "start", "end", "rownames")]
colnames(homer_peaks) <- c("Chr", "Start", "End", "Geneid")
fwrite(homer_peaks, file=paste0(OUT_DIR,"CMP_homer_DE_peaks.bed"), row.names=FALSE, col.names=FALSE, sep = "\t") 

# Get Significant peaks for HSCs --------------------------------------------------

HSC_clusters <- clusters[which(clusters_renamed %in% c("HSC1", "HSC2"))]
options <- c(1,1)
HSC_peaks_all <- vector()

for(i in 1:length(HSC_clusters)){
  name <- HSC_clusters[i]
  
  res_full <- read.csv(file=paste0("0X_test_filter_thresholds/res_full/", name, "_res_full_option",options[i] ,".csv"), header=TRUE)
  peaks_sig <- subset(res_full, res_full[,"BH_corrected"] < 0.1)$X
  
  HSC_peaks_all <- append(HSC_peaks_all, peaks_sig)
}

HSC_peaks_unique <- unique(HSC_peaks_all)

peakmat_subset <- subset(peakmat, rowRanges(peakmat)$rownames %in% as.character(HSC_peaks_unique))

df_peakmat_subset <- as.data.frame(rowRanges(peakmat_subset))

homer_peaks <- df_peakmat_subset[,c("seqnames", "start", "end", "rownames")]
colnames(homer_peaks) <- c("Chr", "Start", "End", "Geneid")
fwrite(homer_peaks, file=paste0(OUT_DIR,"HSC_homer_DE_peaks.bed"), row.names=FALSE,col.names=FALSE, sep = "\t") 


# Get Significant peaks for GMPs --------------------------------------------------

GMP_clusters <- clusters[which(clusters_renamed %in% c("GMP1","GMP2", "GMP3"))]
options <- c(1,2,1)
GMP_peaks_all <- vector()

for(i in 1:length(GMP_clusters)){
  name <- GMP_clusters[i]
  
  res_full <- read.csv(file=paste0("0X_test_filter_thresholds/res_full/", name, "_res_full_option",options[i] ,".csv"), header=TRUE)
  peaks_sig <- subset(res_full, res_full[,"BH_corrected"] < 0.1)$X
  
  GMP_peaks_all <- append(GMP_peaks_all, peaks_sig)
}

GMP_peaks_unique <- unique(GMP_peaks_all)

peakmat_subset <- subset(peakmat, rowRanges(peakmat)$rownames %in% as.character(GMP_peaks_unique))

df_peakmat_subset <- as.data.frame(rowRanges(peakmat_subset))

homer_peaks <- df_peakmat_subset[,c("seqnames", "start", "end", "rownames")]
colnames(homer_peaks) <- c("Chr", "Start", "End", "Geneid")
fwrite(homer_peaks, file=paste0(OUT_DIR,"GMP_homer_DE_peaks.bed"), row.names=FALSE,col.names=FALSE, sep = "\t")

# Get Significant peaks for MEPs --------------------------------------------------

MEP_clusters <- clusters[which(clusters_renamed %in% c("MEP1", "MEP2", "MEP3", "MEP4"))]
options <- c(1,1,1.25,1.75)
MEP_peaks_all <- vector()

for(i in 1:length(MEP_clusters)){
  name <- MEP_clusters[i]
  
  res_full <- read.csv(file=paste0("0X_test_filter_thresholds/res_full/", name, "_res_full_option",options[i] ,".csv"), header=TRUE)
  peaks_sig <- subset(res_full, res_full[,"BH_corrected"] < 0.1)$X
  
  MEP_peaks_all <- append(MEP_peaks_all, peaks_sig)
}

MEP_peaks_unique <- unique(MEP_peaks_all)

peakmat_subset <- subset(peakmat, rowRanges(peakmat)$rownames %in% as.character(MEP_peaks_unique))

df_peakmat_subset <- as.data.frame(rowRanges(peakmat_subset))

homer_peaks <- df_peakmat_subset[,c("seqnames", "start", "end", "rownames")]
colnames(homer_peaks) <- c("Chr", "Start", "End", "Geneid")
fwrite(homer_peaks, file=paste0(OUT_DIR,"MEP_homer_DE_peaks.bed"), row.names=FALSE, col.names=FALSE, sep = "\t")


# Get Significant peaks for CLP --------------------------------------------------

CLP_clusters <- clusters[which(clusters_renamed %in% c("CLP"))]
options <- c(1.25)
CLP_peaks_all <- vector()

for(i in 1:length(CLP_clusters)){
  name <- CLP_clusters[i]
  
  res_full <- read.csv(file=paste0("0X_test_filter_thresholds/res_full/", name, "_res_full_option",options[i] ,".csv"), header=TRUE)
  peaks_sig <- subset(res_full, res_full[,"BH_corrected"] < 0.1)$X
  
  CLP_peaks_all <- append(CLP_peaks_all, peaks_sig)
}

CLP_peaks_unique <- unique(CLP_peaks_all)

peakmat_subset <- subset(peakmat, rowRanges(peakmat)$rownames %in% as.character(CLP_peaks_unique))

df_peakmat_subset <- as.data.frame(rowRanges(peakmat_subset))

homer_peaks <- df_peakmat_subset[,c("seqnames", "start", "end", "rownames")]
colnames(homer_peaks) <- c("Chr", "Start", "End", "Geneid")
fwrite(homer_peaks, file=paste0(OUT_DIR,"CLP_homer_DE_peaks.bed"), row.names=FALSE,col.names=FALSE, sep = "\t")

# Get Significant peaks for PreBNK --------------------------------------------------

PreBNK_clusters <- clusters[which(clusters_renamed %in% c("Pre-BNK"))]
options <- c(1.75)
PreBNK_peaks_all <- vector()

for(i in 1:length(PreBNK_clusters)){
  
  name <- PreBNK_clusters[i]
  
  res_full <- read.csv(file=paste0("0X_test_filter_thresholds/res_full/", name, "_res_full_option",options[i] ,".csv"), header=TRUE)
  peaks_sig <- subset(res_full, res_full[,"BH_corrected"] < 0.1)$X
  
  PreBNK_peaks_all <- append(PreBNK_peaks_all, peaks_sig)
}

PreBNK_peaks_unique <- unique(PreBNK_peaks_all)

peakmat_subset <- subset(peakmat, rowRanges(peakmat)$rownames %in% as.character(PreBNK_peaks_unique))

df_peakmat_subset <- as.data.frame(rowRanges(peakmat_subset))

homer_peaks <- df_peakmat_subset[,c("seqnames", "start", "end", "rownames")]
colnames(homer_peaks) <- c("Chr", "Start", "End", "Geneid")
fwrite(homer_peaks, file=paste0(OUT_DIR,"PreBNK_homer_DE_peaks.bed"), row.names=FALSE, col.names=FALSE, sep = "\t")
