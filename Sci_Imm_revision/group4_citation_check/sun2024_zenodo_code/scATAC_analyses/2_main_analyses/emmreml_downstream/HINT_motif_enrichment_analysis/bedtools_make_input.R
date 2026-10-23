library(stringr)
library(textTinyR)
library(pbapply)
library(readr)
library(data.table)

FILE_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/match/"
OUT_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/bedtools_input/"
dir.create(OUT_DIR)

clusters <- c("HSC", "CMP", "GMP", "MEP", "CLP", "PreBNK")

for(i in 1:length(clusters)){
  
  file <- read.table(file=paste0(FILE_DIR, clusters[i],"_homer_DE_peaks_mpbs.bed"))
  
  peaks_file <- file[,c("V1", "V2", "V3", "V4")]
  
  
  colnames(peaks_file) <- c("Chr", "Start", "End", "Geneid")
  
  fwrite(peaks_file, file=paste0(OUT_DIR,clusters[i],"_HINT_mpbs.bed"), row.names=FALSE,sep = "\t") 
  
}

##Background peaks
file <- read.table(file=paste0(FILE_DIR,"homer_background_peaks_mpbs.bed"))

peaks_file <- file[,c("V1", "V2", "V3", "V4")]


colnames(peaks_file) <- c("Chr", "Start", "End", "Geneid")

fwrite(peaks_file, file=paste0(OUT_DIR,"background_HINT_mpbs.bed"), row.names=FALSE,sep = "\t") 
