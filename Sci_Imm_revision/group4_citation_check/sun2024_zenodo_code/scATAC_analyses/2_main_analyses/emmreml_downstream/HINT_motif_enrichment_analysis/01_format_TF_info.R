
#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

.libPaths("/project/lbarreiro/USERS/sarah/Rlibs_new")

library(stringr)
library(textTinyR)
library(pbapply)
library(readr)
library(data.table)
library(ggplot2)
library(ggrepel)
library(ggpubr)

setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/")
OUT_DIR <- "emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/"
dir.create(OUT_DIR)

# Cluster and TF info -----------------------------------------------
clusters <- c("HSC", "CMP", "GMP", "MEP", "CLP", "PreBNK")
TF_info <- data.frame(read.csv(file=paste0("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis1_Clusters_harmony/scHINT_subtypes_combined/interaction_model/TF_cluster_information.csv"), header=TRUE))
TF_clusters <- unique(TF_info$clusters)
#Add a name to specific TFs of interest
TF_info$TF_clusters_name <- rep("", length(TF_info$clusters))
TF_info$TF_clusters_name[which(TF_info$clusters == 103)] <- "ZNF/MAZ"
TF_info$TF_clusters_name[which(TF_info$clusters == 16)] <- "KLF/SP_1"
TF_info$TF_clusters_name[which(TF_info$clusters == 39)] <- "KLF/SP_2"
TF_info$TF_clusters_name[which(TF_info$clusters == 40)] <- "ETS_1"
TF_info$TF_clusters_name[which(TF_info$clusters == 27)] <- "ETS_2"
TF_info$TF_clusters_name[which(TF_info$clusters == 57)] <- "KLF17/ZNF354C"
TF_info$TF_clusters_name[which(TF_info$clusters == 49)] <- "TCF/SNA/MYB/ZEB1/FIGLA"
TF_info$TF_clusters_name[which(TF_info$clusters == 188)] <- "ZNF135/460"
TF_info$TF_clusters_name[which(TF_info$clusters == 66)] <- "ZNF24/EWSR1-FL1"
TF_info$TF_clusters_name[which(TF_info$clusters == 20)] <- "IRF1/2"
TF_info$TF_clusters_name[which(TF_info$clusters == 173)] <- "KLF2/3/6"
TF_info$TF_clusters_name[which(TF_info$clusters == 37)] <- "ETS_3"
TF_info$TF_clusters_name[which(TF_info$clusters == 62)] <- "CTCF"
TF_info$TF_clusters_name[which(TF_info$clusters == 74)] <- "EGR1/2/3/4"
TF_info$TF_clusters_name[which(TF_info$clusters == 122)] <- "IRF3-9"
TF_info$TF_clusters_name[which(TF_info$clusters == 170)] <- "FOS/JUN_1"
TF_info$TF_clusters_name[which(TF_info$clusters == 7)] <- "CEBP_1"
TF_info$TF_clusters_name[which(TF_info$clusters == 13)] <- "GATA"
TF_info$TF_clusters_name[which(TF_info$clusters == 51)] <- "NFkB"
TF_info$TF_clusters_name[which(TF_info$clusters == 98)] <- "RUNX"
TF_info$TF_clusters_name[which(TF_info$clusters == 109)] <- "HOX"
TF_info$TF_clusters_name[which(TF_info$clusters == 48)] <- "FOS/JUN_2"
TF_info$TF_clusters_name[which(TF_info$clusters == 107)] <- "ATF3/FOS/JUN/BACH"
TF_info$TF_clusters_name[which(TF_info$clusters == 43)] <- "NFE2"
TF_info$TF_clusters_name[which(TF_info$clusters == 81)] <- "CEBP_2"
TF_info$TF_clusters_name[which(TF_info$clusters == 5)] <- "CREB_JUN"
TF_info$TF_clusters_name[which(TF_info$clusters == 9)] <- "ERG/ETV/ETS"
TF_info$TF_clusters_name[which(TF_info$clusters == 31)] <- "PAX"
TF_info$TF_clusters_name[which(TF_info$clusters == 61)] <- "REST"
TF_info$TF_clusters_name[which(TF_info$clusters == 84)] <- "FOS/JUN_3"
TF_info$TF_clusters_name[which(TF_info$clusters == 104)] <- "HOXA9"
TF_info$TF_clusters_name[which(TF_info$clusters == 123)] <- "NFE2/BACH2"
TF_info$TF_clusters_name[which(TF_info$clusters == 152)] <- "TGIF"
TF_info$TF_clusters_name[which(TF_info$clusters == 46)] <- "USF1/2"
TF_info$TF_clusters_name[which(TF_info$clusters == 114)] <- "TCFL5"
TF_info$TF_clusters_name[which(TF_info$clusters == 83)] <- "E2F6"
TF_info$TF_clusters_name[which(TF_info$clusters == 93)] <- "NRF1"
TF_info$TF_clusters_name[which(TF_info$clusters == 101)] <- "TFAP"
TF_info$TF_clusters_name[which(TF_info$clusters == 147)] <- "ZNF740"
TF_info$TF_clusters_name[which(TF_info$clusters == 119)] <- "ELF"
TF_info$TF_clusters_name[which(TF_info$clusters == 144)] <- "GLIS"
TF_info$TF_clusters_name[which(TF_info$clusters == 1)] <- "TFAP2A/B/C/E"
TF_info$TF_clusters_name[which(TF_info$clusters == 47)] <- "YY1/2, ZFP42"
TF_info$TF_clusters_name[which(TF_info$clusters == 58)] <- "HINFP"
TF_info$TF_clusters_name[which(TF_info$clusters == 75)] <- "PLAG1/L2"
TF_info$TF_clusters_name[which(TF_info$clusters == 102)] <- "ZBTB33/ZBED1"
TF_info$TF_clusters_name[which(TF_info$clusters == 110)] <- "HEY1/2, HES1/2"
TF_info$TF_clusters_name[which(TF_info$clusters == 194)] <- "ZNF682"
TF_info$TF_clusters_name[which(TF_info$clusters == 155)] <- "PROX1"
TF_info$TF_clusters_name[which(TF_info$clusters == 55)] <- "NR2/RXR"
TF_info$TF_clusters_name[which(TF_info$clusters == 69)] <- "INSM1"
TF_info$TF_clusters_name[which(TF_info$clusters == 137)] <- "ZIC1/3/4/5"
TF_info$TF_clusters_name[which(TF_info$clusters == 157)] <- "TFAP2"
TF_info$TF_clusters_name[which(TF_info$clusters == 196)] <- "ZKSCAN5"
TF_info$TF_clusters_name[which(TF_info$clusters == 23)] <- "MZF1"
TF_info$TF_clusters_name[which(TF_info$clusters == 197)] <- "ZNF16"
write.csv(TF_info, file=paste0(OUT_DIR, "TF_cluster_information_with_names.csv"))


