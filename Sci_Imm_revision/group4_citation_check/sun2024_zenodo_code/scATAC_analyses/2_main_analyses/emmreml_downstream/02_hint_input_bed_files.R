library(ArchR)
addArchRThreads(threads = 16) 
addArchRGenome("hg38")
library(readr)
library(harmony)
library(Seurat)
library(stringr)
library(textTinyR)
library(pbapply)

# setup -----------------------------------------------------------------------------------
setwd("/scATAC_analyses/2_main_analyses/emmreml_downstream")
OUT_DIR <- "hint_input_bed_files/"
dir.create(OUT_DIR)

projBM <- loadArchRProject(path = "../Unedited_with_Clusters_harmony_peaks_combined", force = FALSE, showLogo = TRUE)
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




# Get sig peaks for each cluster ------------------------------------------


clusters <- c("C3", "C4", "C5", "C6", "C7","C9", "C10", "C12", "C14", "C15", "C17", "C18",  "C21", "C22","C23", "C24")
clusters_renamed <- c("CMP1", "CMP2", "HSC1", "HSC2", "MEP1","Pre-BNK", "GMP1","CLP", "MEP2", "MEP3", "CMP3", "MEP4", "CMP4", "CMP5", "GMP2", "GMP3")
options <- c(1,1,1,1,1,1.75, 1, 1.25, 1, 1.25, 2.5, 1.75, 2.25, 1, 2, 1)

peaks_all <- vector()


for(i in 1:length(clusters)){
  name <- clusters[i]
  
  res_full <- read.csv(file=paste0("../0X_test_filter_thresholds/res_full/",name,"_res_full_option",options[i] ,".csv"), header=TRUE)
  peaks_sig <- subset(res_full, res_full[,"BH_corrected"] < 0.1)$X
  
  peakmat_subset <- subset(df_peakmat, df_peakmat$rownames %in% peaks_sig)
  
  df_peakmat_subset <- as.data.frame(peakmat_subset)
  
  homer_peaks <- df_peakmat_subset[,c("seqnames", "start", "end", "rownames")]
  colnames(homer_peaks) <- c("Chr", "Start", "End", "Geneid")
  fwrite(homer_peaks, file=paste0(OUT_DIR, clusters_renamed[i],"_homer_DE_peaks.bed"), row.names=FALSE, col.names=FALSE, sep = "\t") 
  
}

