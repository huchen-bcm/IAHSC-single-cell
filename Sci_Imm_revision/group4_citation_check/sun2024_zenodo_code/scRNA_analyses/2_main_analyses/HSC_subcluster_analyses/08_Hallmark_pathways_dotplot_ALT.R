
#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

library(ggplot2)
library(DESeq2)
library(statmod)
library(RColorBrewer)
library(edgeR)
library(devtools)
library(cluster)
library(ggrepel)
library(EMMREML)
library(ggpointdensity)
library(ashr)
library(mashr)
library(dplyr)
library(plyr)
library(grid)
library(gridExtra)
library(Seurat)

# Here we make a scatterplot that correlates HSC subcluster MPP score with gsea pathway enrichments 
# gsea pathways are classified on a spectrum from HSC specific to universal depending on whether they were enriched only in HSCs or in the majority of 
HSPC clusters 
# We show that universally shared pathways are primarily enriched within the MPP-like subtypes of HSCs 


# setup ----------------------------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/HSC_subcluster_analyses/")
OUT_DIR <- c("Hallmark_pathways_dotplot/")
dir.create(OUT_DIR)

# read in gsea results for each cluster ------------------------------------------------
clusters <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", 'c7', "c8")

target_pathways <- c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_APOPTOSIS",
                     "HALLMARK_ALLOGRAFT_REJECTION","HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_DNA_REPAIR",
                     "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_UV_RESPONSE_UP", "HALLMARK_HEME_METABOLISM", "HALLMARK_KRAS_SIGNALING_DN")

MPP_data <- read.csv(file="MPP_cluster_data.csv")


# subset out pathways that are target pathways --------------------------
# "Target pathways" are those that were significant within at least one cluster of the HSPC gsea analysis (see Fig 2C of the paper)
for(i in 1:length(clusters)){
  dat <- data.frame(read.csv(file=paste0("9_gsea/NOT_SORTED_gsea_result_",clusters[i],"_no_correlations")))
  dat_subset <- dat[which(dat$pathways %in% target_pathways),]
  names <- dat_subset$pathways
  assign(paste0(clusters[i]), as.numeric(dat_subset$padj))
}
final_data <- cbind(c0, c1, c2, c3, c4, c5, c6, c7, c8)
row.names(final_data) <- names


gsea_file <- "gsea_padj_plots/Figure_2B_gsea.rds"
gsea <- readRDS(gsea_file)
gsea_original <- gsea$data


score <- vector(length=length(row.names(final_data)))
count <- vector(length=length(row.names(final_data)))

for(i in 1:length(row.names(final_data))){
  pathway_name <- row.names(final_data)[i]
  pathway_name_final <- gsub("HALLMARK_", "", pathway_name)
  values <- head(sort(final_data[i,]),3)
  clusters <- colnames(final_data)[which(final_data[i,] %in% values)]
  clusters <- gsub("c", "g", clusters)
  averages <- MPP_data$average[which(MPP_data$X %in% clusters)]
  score[i] <- mean(averages)
  
  temp <- subset(gsea_original, gsea_original$pathway_name %in% pathway_name_final)
  count[i] <- length(which(temp$padj < 0.1))
}

plot_data <- data.frame(
  p <- gsub("HALLMARK_", "", row.names(final_data)),
  c <- count,
  s <- score
)
colnames(plot_data) <- c("p", "c", "s")



final_data_ordered <- final_data[c("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MYC_TARGETS_V2", "HALLMARK_DNA_REPAIR",
                                "HALLMARK_HYPOXIA", "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_HEME_METABOLISM", "HALLMARK_APOPTOSIS", "HALLMARK_IL2_STAT5_SIGNALING", 
                                "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_UV_RESPONSE_UP",  "HALLMARK_KRAS_SIGNALING_DN"),c("c0", "c4", "c5", "c7", "c3", "c1", "c6", "c2", "c8")]
final_data_copy <- final_data_ordered
for (i in 1:length(row.names(final_data_ordered))){
  index <- as.numeric(which(final_data_copy[i,] == min(final_data_copy[i,])))
  final_data_copy[i,] <- 0
  final_data_copy[i, index] <- 1
}


plot_data <- data.frame(
  cluster_name <- rep(colnames(final_data_copy),each = length(row.names(final_data_copy))),
  pathway_name <- rep(row.names(final_data_copy), length(colnames(final_data_copy))),
  vals <- c(final_data_copy[,"c0"], final_data_copy[,"c4"],final_data_copy[,"c5"],final_data_copy[,"c7"],final_data_copy[,"c3"],
            final_data_copy[,"c1"],final_data_copy[,"c6"],final_data_copy[,"c2"],final_data_copy[,"c8"])
)
colnames(plot_data) <- c("cluster_name", "pathway_name", "vals")

plot_data$pathway_name <- gsub("HALLMARK_", "", plot_data$pathway_name)
plot_data$pathway_name <- factor(plot_data$pathway_name, c("OXIDATIVE_PHOSPHORYLATION", "TNFA_SIGNALING_VIA_NFKB", "REACTIVE_OXYGEN_SPECIES_PATHWAY", "MYC_TARGETS_V1", "MYC_TARGETS_V2", "DNA_REPAIR",
                                   "HYPOXIA", "ALLOGRAFT_REJECTION", "HEME_METABOLISM", "APOPTOSIS", "IL2_STAT5_SIGNALING", 
                                   "EPITHELIAL_MESENCHYMAL_TRANSITION", "UV_RESPONSE_UP",  "KRAS_SIGNALING_DN"))

plot_data$cluster_name <- factor(plot_data$cluster_name, c("c0", "c4", "c5", "c7", "c3", "c1", "c6", "c2", "c8"))



p <- ggplot(plot_data, aes(x=pathway_name, y = cluster_name)) + geom_point(shape=21, color="white", aes(size=vals, fill=vals)) + theme_test() +
  scale_fill_gradient(low="white", high="#A50021") + labs(x="", y = "") + theme(legend.position = "none") + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=5))


ggsave(filename = file.path(OUT_DIR, "HSC_subcluster_dotplot_simple.pdf"), 
       plot = p, 
       height = 4 , width = 4)
