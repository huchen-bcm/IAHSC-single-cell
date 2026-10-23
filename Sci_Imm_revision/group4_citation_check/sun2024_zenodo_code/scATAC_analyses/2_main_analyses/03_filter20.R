
library(ArchR)
addArchRThreads(threads = 16) 
addArchRGenome("hg38")
library(readr)
library(harmony)
library(Seurat)
library(stringr)
library(textTinyR)
library(pbapply)

# setup --------------------------------------------------------------
setwd("/scATAC_analyses/2_main_analyses")
OUT_DIR <- "filter20/"
dir.create(OUT_DIR)



# Read in pseudobulk stats ------------------------------------------------

celltypes <- c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8", "C9", "C10", "C11", "C12", "C13", "C14", "C15", "C16", "C17", "C18", "C19", "C20", "C21", "C22", "C23", "C24")
for(i in 1:length(celltypes))
{
  df <- readRDS(file=paste0("pseudobulk_clusters_harmony_all_combined/", celltypes[i], "_stats_df.rds"))
  df_subset <- subset(df, as.numeric(df$cellcount) <20)
  names <- rownames(df_subset)
  assign(paste0(celltypes[i],"_filter20"), names)
}

filter_20_list <- list(C1_filter20, C2_filter20, C3_filter20, C4_filter20, C5_filter20, C6_filter20, C7_filter20, C8_filter20, C9_filter20, C10_filter20, C11_filter20, C12_filter20, C13_filter20, C14_filter20, C15_filter20, C16_filter20, C17_filter20, C18_filter20, C19_filter20, C20_filter20, C21_filter20, C22_filter20, C23_filter20, C24_filter20)
saveRDS(filter_20_list, file=paste0(OUT_DIR,"filter_20_list.rds"))
