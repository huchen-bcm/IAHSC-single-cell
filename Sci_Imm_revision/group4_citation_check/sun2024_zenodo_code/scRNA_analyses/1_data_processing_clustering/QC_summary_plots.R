
#R version 4.1.0
library("ggplot2")
library ("DESeq2")
library(statmod)
library("RColorBrewer")
library(edgeR)
library(devtools)
library(cluster)
library(ggrepel)
library(ggpubr)

setwd("scRNA_analyses/1_data_processing_clustering/")
OUT_DIR <- "QC_summary_plots/"
dir.create(OUT_DIR)


# Read in QC stats tables -------------------------------------------------
percent_mt_stats <- read.table(file="percent_mt_stats_unfiltered.txt")
nRNA_nFeature_stats <- read.table(file="nRNA_nFeature_stats_unfiltered.txt")
cell_nums_stats <- read.table(file="cell_numbers_after_QC_summary.txt")


# Plot percent mt stats ---------------------------------------------------
percent_mt_df <- data.frame(
  capture <- c("capture1", "capture2", "capture3", "capture4", "capture5", "capture6", "capture7", "capture8", "capture9", "capture10", "capture11", "capture12", "capture13", "capture14"),
  timepoint <- c("mixed", "Tm3", "Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3"),
  percent_mt <- percent_mt_stats$average_percent_mt
)
colnames(percent_mt_df) <- c("capture","timepoint", "percent_mt")

percent_mt_df$timepoint <- factor(percent_mt_df$timepoint, c("Td0","Tm3","mixed"))
tiff(file=paste0(OUT_DIR,"percent_mt_stats_before_filtering_barplot.tiff"), units="in", width=6, height=6, res=350)
ggbarplot(percent_mt_df, x = "timepoint", y="percent_mt", fill="timepoint", color="black",
            add = c("mean_sd","dotplot"))+ theme_test()+labs(main="Percent mt stats before filtering", x="",y="Mean percent mt reads")
dev.off()


# Plot percent nRNA stats -------------------------------------------------
nRNA_df <- data.frame(
  capture <- c("capture1", "capture2", "capture3", "capture4", "capture5", "capture6", "capture7", "capture8", "capture9", "capture10", "capture11", "capture12", "capture13", "capture14"),
  timepoint <- c("mixed", "Tm3", "Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3"),
  nRNA <- nRNA_nFeature_stats$average_nRNA
)
colnames(nRNA_df) <- c("capture","timepoint", "nRNA")

nRNA_df$timepoint <- factor(nRNA_df$timepoint, c("Td0","Tm3","mixed"))
tiff(file=paste0(OUT_DIR,"nRNA_stats_before_filtering_barplot.tiff"), units="in", width=6, height=6, res=350)
ggbarplot(nRNA_df, x = "timepoint", y="nRNA", fill="timepoint", color="black",
          add = c("mean_sd","dotplot"))+ theme_test()+labs(main="nRNA stats before filtering", x="",y="Mean nRNA Count")
dev.off()


# Plot nFeature stats -----------------------------------------------------
nFeature_df <- data.frame(
  capture <- c("capture1", "capture2", "capture3", "capture4", "capture5", "capture6", "capture7", "capture8", "capture9", "capture10", "capture11", "capture12", "capture13", "capture14"),
  timepoint <- c("mixed", "Tm3", "Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3"),
  nFeature <- nRNA_nFeature_stats$average_nFeature
)
colnames(nFeature_df) <- c("capture","timepoint", "nFeature")

nFeature_df$timepoint <- factor(nFeature_df$timepoint, c("Td0","Tm3","mixed"))
tiff(file=paste0(OUT_DIR,"nFeature_stats_before_filtering_barplot.tiff"), units="in", width=6, height=6, res=350)
ggbarplot(nFeature_df, x = "timepoint", y="nFeature", fill="timepoint", color="black",
          add = c("mean_sd","dotplot"))+ theme_test()+labs(main="nFeature stats before filtering", x="",y="Mean nFeature Count")

dev.off()



# Plot cell nums before and after -----------------------------------------
cell_nums_df <- data.frame(
  capture <- rep(c("capture1", "capture2", "capture3", "capture4", "capture5", "capture6", "capture7", "capture8", "capture9", "capture10", "capture11", "capture12", "capture13", "capture14"),3),
  value <- c(cell_nums_stats$cell_nums_unfiltered,cell_nums_stats$cell_nums_singlets,cell_nums_stats$cell_nums_singlets_QCfiltered),
  category <- c(rep("before",14), rep("before_single_cell",14), rep("after_QC",14)),
  timepoint <- c(rep("",28),"mixed", "Tm3", "Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3","Td0", "Tm3")
)
colnames(cell_nums_df) <- c("capture","value", "category","timepoint")

cell_nums_df$category <- factor(cell_nums_df$category, c("before", "before_single_cell", "after_QC"))
tiff(file=paste0(OUT_DIR,"cell_nums_stats_before_after.tiff"), units="in", width=6, height=6, res=350)
ggplot(cell_nums_df, aes(x=category, y=value, group=capture)) +
  geom_line(aes(color=capture))+
  geom_point(aes(color=capture))+geom_label_repel(aes(label = timepoint),
                                                  nudge_x = 1,
                                                  na.rm = TRUE)
dev.off()
