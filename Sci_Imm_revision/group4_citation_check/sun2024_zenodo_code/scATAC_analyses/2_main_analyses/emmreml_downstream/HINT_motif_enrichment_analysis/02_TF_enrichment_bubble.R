
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
TF_info <- data.frame(read.csv(file=paste0("emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/TF_cluster_information_with_names.csv"), header=TRUE))
TF_clusters <- unique(TF_info$clusters)


TF <- vector()
cell <- vector()
percentage <- vector()
fdr <- vector()
count <- 1
for(i in 1:length(clusters)){
  #Read in statistical enrichment data
  enrichment_data <- read.table(file=paste0("emmreml_downstream/HINT_motif_enrichment_analysis/hint_Enrichment/",clusters[i],"_homer_DE_peaks/fulltest_statistics.txt"), header=TRUE, fill=TRUE)
  
  #Read in bedtools overlap data
  bedtools_overlap <- read.table(file=paste0("emmreml_downstream/HINT_motif_enrichment_analysis/bedtools_output/peaks_with_motifs/", clusters[i],"_HINT_out.bed"))
  
  #Get the percent DA peaks bound by each enriched TF
  tfs <- enrichment_data$FACTOR
  for(j in 1:length(TF_clusters)){
    #Get the most enriched member of each TF cluster from TF info
    motifs <- TF_info$TF.motif[which(TF_info$clusters == TF_clusters[j])]
    subset <- subset(enrichment_data, enrichment_data$FACTOR %in% motifs)
    tf_specific <- subset$FACTOR[which.min(subset$CORR.P.VALUE)]
    fdr[count] <- subset$CORR.P.VALUE[which.min(subset$CORR.P.VALUE)] #min pval for TFs in the group
    
    num_DA_peaks <- unique(bedtools_overlap$V4[which(bedtools_overlap$V8 %in% tf_specific)])
    total_peaks <- unique(bedtools_overlap$V4)
    percentage[count] <- length(num_DA_peaks)/length(total_peaks)
    
    TF[count] <- TF_clusters[j]
    if(clusters[i]=="CLP"){
      cell[count] <- "MLP"
    }
    else{
      cell[count] <- clusters[i]
    }
    
    count <- count + 1
  }
  
  hint_enrich <- data.frame(cbind(TF, fdr))
  colnames(hint_enrich) <- c("TF","fdr")
  write.csv(hint_enrich, file=paste0(OUT_DIR,clusters[i],"_hint_enrich_fdr_vals.csv"))
  
  print(i)
}

TF_renamed <- vector(length=length(TF))
for(i in 1:length(TF)){
  TF_renamed[i] <- unique(TF_info$TF_clusters_name[which(TF_info$clusters == TF[i])])
}

df_table <- data.frame(
  X <- as.character(TF_renamed),
  Y <- cell,
  p <- -log10(fdr),
  perc <- percentage,
  tf_num <- TF
)
colnames(df_table) <- c("X", "Y", "p", "perc", "tf_num")
write.csv(df_table, file=paste0(OUT_DIR, "FDR_percentage_FULL_table.csv"))


###########################################################
### get the top enriched TF families within CMP and GMP ###
###########################################################

data_ordered <- df_table[order(df_table$p, decreasing=TRUE),]





subset1 <- subset(df_table,df_table$p > (-log10(0.0005)) )
tfs_to_keep <- subset(subset1,subset1$perc >=0.15 )$tf_num

final_df <- subset(df_table, df_table$tf_num %in% unique(tfs_to_keep))
#colnames(df_table) <- c("X", "Y", "p", "perc")

p_size <- final_df$p
p_size[which(p_size < (-log10(0.0005)))] <- 0

p_color <- rep("black", length(final_df$p))
p_color[which(p_size < (-log10(0.0005)))] <- "white"

final_df$p_size <- p_size
final_df$p_color <- p_color

final_df$Y <- factor(final_df$Y, c("HSC", "CMP", "GMP", "MEP","PreBNK", "MLP"))

# tiff(paste0(OUT_DIR,"TF_enrichemnts_all_clusters_bubble.tiff"), units="in", width=2.5, height=6, res=300)
# ggplot(final_df, aes(x=X, y=Y, size=p_size, fill=p_size, color=p_color)) + geom_point(shape=21) + theme_classic() + labs(x="", y="") + scale_fill_gradient(high="red", low="white") + 
#   theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90))+ coord_flip() + scale_color_manual(values=c("black","white"))
# dev.off()
# 
# #same plot but with legend
# tiff(paste0(OUT_DIR,"TF_enrichemnts_all_clusters_bubble_wlegend.tiff"), units="in", width=3.8, height=6, res=350)
# ggplot(final_df, aes(x=X, y=Y,  fill=p_size,size=p_size, color=p_color)) + geom_point(shape=21) + theme_classic() + labs(x="", y="", fill="-log10(FDR)", size="-log10(FDR)") + scale_fill_gradient(high="red", low="white") + 
#   theme(axis.text.x = element_text(angle = 90))+ coord_flip() + scale_color_manual(values=c("black","white")) + guides(color = FALSE)   
# dev.off()

#Horizontal version
final_df$X <- as.character(final_df$X)
final_df$X[which(final_df$tf_num == 49)] <- "TCF/SNA/MYB"
tiff(paste0(OUT_DIR,"TF_enrichemnts_all_clusters_bubble_wlegend_HORIZONTAL.tiff"), units="in", width=6.7, height=3.6, res=400)
ggplot(final_df, aes(x=X, y=Y,  fill=p_size,size=p_size, color=p_color)) + geom_point(shape=21) + theme_test() + labs(x="", y="", fill="-log10(FDR)", size="-log10(FDR)") + scale_fill_gradient(high="red", low="white") + 
  theme(axis.text.x = element_text(size=10, angle = 45, hjust=1, vjust=1))+ scale_color_manual(values=c("black","white")) + guides(color = FALSE)  + scale_y_discrete(limits=rev) + theme(legend.position="bottom")
dev.off()

#save same plot as RDS
bubble_data <- final_df ##copy of final_df for this version of the plot
colnames(bubble_data) <- c("X", "Y", "-log10FDR", "percent_sites", "TF_ID", "point_size", "point_color")
p <- ggplot(bubble_data, aes(x=X, y=Y,  fill=point_size,size=point_size, color=point_color)) + geom_point(shape=21) + theme_test() + labs(x="", y="", fill="-log10(FDR)", size="-log10(FDR)") + scale_fill_gradient(high="red", low="white") + 
  theme(axis.text.x = element_text(size=10, angle = 45, hjust=1, vjust=1))+ scale_color_manual(values=c("black","white")) + guides(color = FALSE)  + scale_y_discrete(limits=rev) + theme(legend.position="bottom")

saveRDS(p, file=paste0(OUT_DIR,"Figure_5A_bubble.rds"))



# pca based on HSC regulon enrichments ------------------------------------

#Read in HSC regulon data
REGULON_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scRNA/Rprojects/projects_version2_Rv4.1/Analysis12_label_transfer_emmreml_edited/pyscenic/get_diff_regulons_outs/"
regulon_HSCa <- read.csv(file=paste0(REGULON_DIR,"HSC_a_wilcox_data_no_outliers.csv"))
regulon_HSCb <- read.csv(file=paste0(REGULON_DIR,"HSC_b_wilcox_data_no_outliers.csv"))
sig_regulons <- c(as.character(regulon_HSCa$X[which(regulon_HSCa$pval_wilcox < 0.15)]), as.character(regulon_HSCb$X[which(regulon_HSCb$pval_wilcox < 0.15)]))
sig_regulons_final <- unique(sig_regulons)

TF_info_copy <- TF_info
str_sub(TF_info_copy$TF.motif, start = 1, end = 9, omit_na = FALSE) <- ""
TF_matched_group <- vector()
count <- 1
for(j in 1:length(sig_regulons_final)){
  if(sig_regulons_final[j] %in% TF_info_copy$TF.motif)
  {
    TF_matched_group[count] <- unique(TF_info_copy$clusters[which(TF_info_copy$TF.motif == sig_regulons_final[j] )])
    count <- count + 1
  }
}

df_pca <- data.frame(
  TF <- TF_clusters,
  c <- cell,
  p <- -log10(fdr)
)
colnames(df_pca) <- c("TF", "c", "p")

df_pca_subset <- subset(df_pca, df_pca$TF %in% TF_matched_group)

pca_matrix <- data.frame(df_pca_subset$p[which(df_pca_subset$c == "CMP")])
pca_matrix <- cbind(pca_matrix, df_pca_subset$p[which(df_pca_subset$c == "GMP")])
pca_matrix <- cbind(pca_matrix, df_pca_subset$p[which(df_pca_subset$c == "MEP")])
pca_matrix <- cbind(pca_matrix, df_pca_subset$p[which(df_pca_subset$c == "MLP")])
pca_matrix <- cbind(pca_matrix, df_pca_subset$p[which(df_pca_subset$c == "PreBNK")])
pca_matrix <- cbind(pca_matrix, df_pca_subset$p[which(df_pca_subset$c == "HSC")])

pca_matrix <- cbind(pca_matrix, rep(0, length(row.names(data.frame(df_pca_subset$p[which(df_pca_subset$c == "CMP")])))))
colnames(pca_matrix) <- c( "CMP", "GMP", "MEP", "MLP", "PreBNK","HSC", "No enrichment")
pca_mat_final <- t(pca_matrix)


pca = prcomp(pca_mat_final)
summary(pca)

pca_values <- pca$x
df <- data.frame(
  pc1 <- as.numeric(pca_values[,1]),
  pc2 <- as.numeric(pca_values[,2]),
  type <- row.names(pca_mat_final)
)
row.names(df) <- row.names(pca_values)
df$type <- factor(df$type, c("HSC", "CMP", "GMP", "MEP", "MLP", "PreBNK", "No enrichment"))

tiff(paste0(OUT_DIR,"PCA_hsc_reg_TF_enrichments.tiff"), units="in", width=3.5, height=3.5, res=300)
ggplot(df,aes(x=pc1,y=pc2)) + theme_classic()+geom_point(aes( fill=type),size = 5,shape=21, stroke = 0.8, alpha=0.8, color="black")+ labs(x="PC1 (93.83%)", y="PC2 (5.56%)")+
  scale_fill_manual(values=c('#EF5350','#FF7F0E', '#33CC00', '#8C5782', '#46B8DA', '#1A0099', "lightgrey")) +  geom_text_repel(aes(label=type), color=c('#FF7F0E', '#33CC00', '#8C5782', '#46B8DA', '#1A0099','#EF5350', "grey")) +
  theme(legend.position = "none")
dev.off()





#PCA based on all enrichments

pca_matrix <- data.frame(df_pca$p[which(df_pca$c == "CMP")])
pca_matrix <- cbind(pca_matrix, df_pca$p[which(df_pca$c == "GMP")])
pca_matrix <- cbind(pca_matrix, df_pca$p[which(df_pca$c == "MEP")])
pca_matrix <- cbind(pca_matrix, df_pca$p[which(df_pca$c == "MLP")])
pca_matrix <- cbind(pca_matrix, df_pca$p[which(df_pca$c == "PreBNK")])


pca_matrix <- cbind(pca_matrix, rep(0, length(row.names(data.frame(df_pca$p[which(df_pca$c == "CMP")])))))
colnames(pca_matrix) <- c( "CMP", "GMP", "MEP", "MLP", "PreBNK","No enrichment")
pca_mat_final <- t(pca_matrix)


pca = prcomp(pca_mat_final)
summary(pca)

pca_values <- pca$x
df <- data.frame(
  pc1 <- as.numeric(pca_values[,1]),
  pc2 <- as.numeric(pca_values[,2]),
  type <- row.names(pca_mat_final)
)
row.names(df) <- row.names(pca_values)
df$type <- factor(df$type, c( "CMP", "GMP", "MEP", "MLP", "PreBNK", "No enrichment"))

tiff(paste0(OUT_DIR,"PCA_TF_enrichments_all.tiff"), units="in", width=2.5, height=2.5, res=300)
ggplot(df,aes(x=pc1,y=pc2)) + theme_classic()+geom_point(aes( fill=type),size = 4,shape=21, stroke = 0.8, alpha=0.8, color="black")+ labs(x="PC1 (95.12%)", y="PC2 (3.03%)")+
  scale_fill_manual(values=c('#FF7F0E', '#33CC00', '#8C5782', '#46B8DA', '#1A0099', "lightgrey")) +  geom_text_repel(aes(label=type), color=c('#FF7F0E', '#33CC00', '#8C5782', '#46B8DA', '#1A0099', "grey")) +
  theme(legend.position = "none")
dev.off()

