#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

library("ggplot2")
library ("DESeq2")
library(statmod)
library("RColorBrewer")
library(edgeR)
library(devtools)
library(cluster)
library(ggrepel)
library(UpSetR)
library(dplyr)
library("gg.gap")
library("stringr")

# In this script we make a preliminary heatmap to summarize the gsea findings 

# setup -----------------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/4_mash_downstream")
OUT_DIR <- "gsea_padj_plots/"
dir.create(OUT_DIR)

# Read in gsea results for all clusters  -----------------------
clusters <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_b1", "GMP_b2", "MEP_a1", "MEP_a2", "MEP_a3", "mixed_a", "mixed_b", "PreBNK")
for(i in 1:length(clusters)){
  name <- clusters[i]
  gsea_data <- read.csv(file=paste0("gsea_padj/NOT_SORTED_gsea_result_",name,"_NO_correlations_with_custom"), header=TRUE)
  assign(paste0(name,"_data"), gsea_data)
}

gsea_result_list <- list(HSC_a_data, HSC_b_data, CMP_a_data, CMP_b_data, CMP_c_data, GMP_b1_data, GMP_b2_data, MEP_a1_data, MEP_a2_data, MEP_a3_data, mixed_a_data, mixed_b_data, PreBNK_data)
gsea_result_names <- list("HSC_a_data", "HSC_b_data", "CMP_a_data", "CMP_b_data", "CMP_c_data", "GMP_a_data", "GMP_b_data", "MEP_a_data", "MEP_b_data", "MEP_c_data", "MLP_a_data", "MLP_b_data", "PreBNK_data")


# Subset out pathways with pval<0.05 --------------------------------------
for(i in 1:length(gsea_result_list)){
  name <- gsea_result_names[i]
  dat <- gsea_result_list[[i]]
  dat_subset <- subset(dat, dat$pval<=0.05)
  assign(paste0(name,"_sig"),dat_subset)
}
gsea_sig_subset_list <- list(HSC_a_data_sig, HSC_b_data_sig, CMP_a_data_sig, CMP_b_data_sig, CMP_c_data_sig, GMP_a_data_sig, GMP_b_data_sig, MEP_a_data_sig, MEP_b_data_sig, MEP_c_data_sig, MLP_a_data_sig, MLP_b_data_sig, PreBNK_data_sig)


# Make a data frame combining results across all conditions ---------------
pathways_all <- vector()
for(i in 1:length(gsea_sig_subset_list)){
  dat <- gsea_sig_subset_list[[i]]
  pathways_all <- append(pathways_all, dat$pathways)
}
pathways <- unique(pathways_all)


df_pval <- data.frame(
  pathway <- HSC_a_data$pathways,
  HSCa_data <- HSC_a_data$NES,
  HSCb_data <- HSC_b_data$NES,
  CMPa_data <- CMP_a_data$NES,
  CMPb_data <- CMP_b_data$NES,
  CMPc_data <- CMP_c_data$NES,
  GMPb1_data <- GMP_b1_data$NES,
  GMPb2_data <- GMP_b2_data$NES,
  MEPa1_data <- MEP_a1_data$NES,
  MEPa2_data <- MEP_a2_data$NES,
  MEPa3_data <- MEP_a3_data$NES,
  mixeda_data <- mixed_a_data$NES,
  mixedb_data <- mixed_b_data$NES,
  Prebnk_data <- PreBNK_data$NES,
  HSCa_pval <- HSC_a_data$pval,
  HSCb_pval <- HSC_b_data$pval,
  CMPa_pval <- CMP_a_data$pval,
  CMPb_pval <- CMP_b_data$pval,
  CMPc_pval <- CMP_c_data$pval,
  GMPb1_pval <- GMP_b1_data$pval,
  GMPb2_pval <- GMP_b2_data$pval,
  MEPa1_pval <- MEP_a1_data$pval,
  MEPa2_pval <- MEP_a2_data$pval,
  MEPa3_pval <- MEP_a3_data$pval,
  mixeda_pval <- mixed_a_data$pval,
  mixedb_pval <- mixed_b_data$pval,
  Prebnk_pval <- PreBNK_data$pval,
  HSCa_padj <- HSC_a_data$padj,
  HSCb_padj <- HSC_b_data$padj,
  CMPa_padj <- CMP_a_data$padj,
  CMPb_padj <- CMP_b_data$padj,
  CMPc_padj <- CMP_c_data$padj,
  GMPb1_padj <- GMP_b1_data$padj,
  GMPb2_padj <- GMP_b2_data$padj,
  MEPa1_padj <- MEP_a1_data$padj,
  MEPa2_padj <- MEP_a2_data$padj,
  MEPa3_padj <- MEP_a3_data$padj,
  mixeda_padj <- mixed_a_data$padj,
  mixedb_padj <- mixed_b_data$padj,
  Prebnk_padj <- PreBNK_data$padj
)

df_pval_final <- subset(df_pval, df_pval$pathway %in% pathways)


colnames(df_pval_final) <- c("pathways", "HSCa_data", "HSCb_data", "CMPa_data", "CMPb_data", "CMPc_data", "GMPb1_data", "GMPb2_data", "MEPa1_data", "MEPa2_data", "MEPa3_data","mixeda_data", "mixedb_data", "Prebnk_data",
                             "HSCa_pval", "HSCb_pval", "CMPa_pval", "CMPb_pval", "CMPc_pval", "GMPb1_pval", "GMPb2_pval", "MEPa1_pval", "MEPa2_pval", "MEPa3_pval","mixeda_pval", "mixedb_pval", "Prebnk_pval",
                             "HSCa_padj", "HSCb_padj", "CMPa_padj", "CMPb_padj", "CMPc_padj", "GMPb1_padj", "GMPb2_padj", "MEPa1_padj", "MEPa2_padj", "MEPa3_padj","mixeda_padj", "mixedb_padj", "Prebnk_padj")


## Clean up pathway names 
df_pval_final$pathways <- str_replace( df_pval_final$pathways,"HALLMARK_", "")


## set pathway text colors and order pathways by group
label_cols <- rep("black", length(df_pval_final$pathways))
path_order <- rep(0, length(df_pval_final$pathways))
order_count <- 1
metabolism_pathways <- c("PI3K","CHOLESTEROL", "HEME", "XENOBIOTIC", "OXIDATIVE", "FATTY", "ADIPOGENESIS", "GLYCOLYSIS", "BILE", "MTORC1", "HYPOXIA", "OXYGEN")
for(i in 1:length(metabolism_pathways)){
  label_cols[which(grepl(metabolism_pathways[i], df_pval_final$pathways)==TRUE)] <- "#D43F3A"
  path_order[which(grepl(metabolism_pathways[i], df_pval_final$pathways)==TRUE)] <- order_count
  order_count <- order_count + 1
}
immune_pathways <- c("COMPLEMENT", "IL2", "TNFA", "INFLAMMATORY", "IL6", "INTERFERON", "TGF_BETA")
for(i in 1:length(immune_pathways)){
  label_cols[which(grepl(immune_pathways[i], df_pval_final$pathways)==TRUE)] <- "#357EBD"
  path_order[which(grepl(immune_pathways[i], df_pval_final$pathways)==TRUE)] <- order_count
  order_count <- order_count + 1
}
proliferation_pathways <- c("MYC", "KRAS", "P53",  "E2F", "APOPTOSIS")
for(i in 1:length(proliferation_pathways)){
  label_cols[which(grepl(proliferation_pathways[i], df_pval_final$pathways)==TRUE)] <- "#5CB85C"
  path_order[which(grepl(proliferation_pathways[i], df_pval_final$pathways)==TRUE)] <- order_count
  order_count <- order_count + 1
}

df_pval_final$labs <- label_cols
df_pval_final$order <- path_order
df_pval_final <- df_pval_final[order(df_pval_final$order),]  



# Visualize results -------------------------------------------------------

df_combined <- data.frame(
  X <- rep(1:length(gsea_sig_subset_list),each=length(pathways)),
  Y <- rep(1:length(pathways),length(gsea_sig_subset_list)),
  NES <- c(as.numeric(df_pval_final$HSCa_data), as.numeric(df_pval_final$HSCb_data), as.numeric(df_pval_final$CMPa_data), as.numeric(df_pval_final$CMPb_data),as.numeric(df_pval_final$CMPc_data), as.numeric(df_pval_final$GMPb1_data), as.numeric(df_pval_final$GMPb2_data), as.numeric(df_pval_final$MEPa1_data), as.numeric(df_pval_final$MEPa2_data), as.numeric(df_pval_final$MEPa3_data), as.numeric(df_pval_final$mixeda_data), as.numeric(df_pval_final$mixedb_data), as.numeric(df_pval_final$Prebnk_data)),
  pval <- c(as.numeric(df_pval_final$HSCa_pval), as.numeric(df_pval_final$HSCb_pval), as.numeric(df_pval_final$CMPa_pval), as.numeric(df_pval_final$CMPb_pval),as.numeric(df_pval_final$CMPc_pval), as.numeric(df_pval_final$GMPb1_pval), as.numeric(df_pval_final$GMPb2_pval), as.numeric(df_pval_final$MEPa1_pval), as.numeric(df_pval_final$MEPa2_pval), as.numeric(df_pval_final$MEPa3_pval), as.numeric(df_pval_final$mixeda_pval), as.numeric(df_pval_final$mixedb_pval), as.numeric(df_pval_final$Prebnk_pval)),
  padj <- c(as.numeric(df_pval_final$HSCa_padj), as.numeric(df_pval_final$HSCb_padj), as.numeric(df_pval_final$CMPa_padj), as.numeric(df_pval_final$CMPb_padj),as.numeric(df_pval_final$CMPc_padj), as.numeric(df_pval_final$GMPb1_padj), as.numeric(df_pval_final$GMPb2_padj), as.numeric(df_pval_final$MEPa1_padj), as.numeric(df_pval_final$MEPa2_padj), as.numeric(df_pval_final$MEPa3_padj), as.numeric(df_pval_final$mixeda_padj), as.numeric(df_pval_final$mixedb_padj), as.numeric(df_pval_final$Prebnk_padj)),
  celltype <- c(rep("HSC_a", length(pathways)), rep("HSC_b", length(pathways)), rep("CMP_a", length(pathways)), rep("CMP_b", length(pathways)), rep("CMP_c", length(pathways)), rep("GMP_a", length(pathways)), rep("GMP_b", length(pathways)), rep("MEP_a", length(pathways)), rep("MEP_b", length(pathways)), rep("MEP_c", length(pathways)), rep("MLP_a", length(pathways)), rep("MLP_b", length(pathways)), rep("PreBNK", length(pathways))),
  pathway_names <- rep(df_pval_final$pathways, length(gsea_sig_subset_list))
)
colnames(df_combined) <- c("X","Y","NES","pval", "padj", "celltype", "pathway_name")

NES_upd_combined <- vector(length=length(df_combined$NES))

for (i in 1:length(df_combined$NES)){
  if(is.na(df_combined$pval[i])==TRUE )
  {
    NES_upd_combined[i]<- 0
  }
  else
  {
    if(as.numeric(df_combined$pval[i]) > 0.05)
    {
      NES_upd_combined[i]<- 0
    }
    else
    {
      NES_upd_combined[i]<- df_combined$NES[i]
    }
  } 
  
  print(i)
}

padj_sig <- rep(0, length(df_combined$padj))
padj_sig[which(df_combined$padj <= 0.1)] <- 1

df_combined_upd1 <- data.frame(cbind(df_combined, NES_upd_combined))
df_combined_upd <- data.frame(cbind(df_combined_upd1, padj_sig))
df_combined_upd$celltype_fill <- df_combined_upd$celltype
df_combined_upd$celltype_fill[which(NES_upd_combined==0)] <- "none"

celltype_color <- df_combined_upd$celltype_fill
celltype_color[which(padj_sig == 0)] <- "none"
df_combined_upd <- data.frame(cbind(df_combined_upd, celltype_color))
colnames(df_combined_upd) <- c("X","Y","NES","pval","padj","celltype","pathway_name", "NES_plot", "padj_plot","celltype_fill", "celltype_color")

myColorScale <- c('HSC_a'='#EF5350','HSC_b'='#CC0C00','CMP_a'='#FF7F0E','CMP_b'='#FFD147','CMP_c'='#CC9900',
                  'GMP_a'='#33CC00', 'GMP_b' ='#99CC00', 'MEP_a'='#8C5782', 'MEP_b'='#990080', 'MEP_c'='#BA68C8',
                  'MLP_a'='#46B8DA', 'MLP_b'='#0288D1', 'PreBNK'='#1A0099', 'none'='white')


p <- ggplot(df_combined_upd, aes(x=X, y=Y, fill=celltype_fill)) + geom_tile(col="black", fill="white", size=0.1)  + 
  geom_point(aes(size = abs(NES_plot), fill=celltype_fill, color=celltype_color), shape=22, alpha=0.6, stroke=2) + theme_test() + scale_fill_manual(values=myColorScale) + scale_color_manual(values=myColorScale)+
  scale_x_discrete(expand=c(0,0), limits = c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_a", "GMP_b", "MEP_a", "MEP_b", "MEP_c", "MLP_a", "MLP_b", "PreBNK"), labels=c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_a", "GMP_b", "MEP_a", "MEP_b", "MEP_c", "MLP_a", "MLP_b", "PreBNK")) +
  scale_y_discrete(expand=c(0,0), limits = df_combined_upd$pathway_name[1:length(unique(df_combined_upd$pathway_name))], labels=df_combined_upd$pathway_name[1:length(unique(df_combined_upd$pathway_name))]) +  theme(axis.text.x = element_text( color="black", angle=45, hjust=1))+
  labs(title = "Gene set enrichment analysis",x="",y = "") + scale_size(range=c(0,5), guide=NULL) + theme(axis.text.y = element_text(colour=df_pval_final$labs)) + theme(legend.position="none")

saveRDS(p, file=paste0(OUT_DIR,"Figure_2B_gsea.rds"))



########
# PCA  #
########

# Here we make a PCA that is based on gsea results for each cluster -----------------------------------

set.seed(831)
pca_data <- df_pval_final[,which(grepl("data", colnames(df_pval_final))==TRUE)]
pval_data <- df_pval_final[,which(grepl("pval", colnames(df_pval_final))==TRUE)]
for(i in 1:length(colnames(pval_data))){
  set_nonsig <- which(pval_data[,i] > 0.05)
  pca_data[set_nonsig,i] <- 0
}
pca_data[,(length(colnames(pca_data))+1)] <- rep(0, length(row.names(pca_data)))

meta_data <- data.frame(
  sample <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_a", "GMP_b", "MEP_a", "MEP_b", "MEP_c", "MLP_a", "MLP_b", "PreBNK", "ref"),
  type <- c("HSC", "HSC", "CMP", "CMP", "CMP", "GMP", "GMP", "MEP", "MEP", "MEP", "mixed", "mixed", "PreBNK", "ref"),
  category <- c(rep("CD34+", 13), "reference"),
  sample_new_name <- c("HSC c1", "HSC c2", "CMP c1", "CMP c2", "CMP c3", "GMP c1", "GMP c2", "MEP c1", "MEP c2", "MEP c3",
                       "MLP c1", "MLP c2", "PreBNK", "ref")
)
colnames(meta_data) <- c("sample", "type", "category", "new_name")

pca = prcomp(t(pca_data))
summary(pca)

pca_values <- pca$x
df <- data.frame(
  pc1 <- as.numeric(pca_values[,1]),
  pc2 <- as.numeric(pca_values[,2]),
  type <- meta_data$type
)
row.names(df) <- row.names(pca_values)

df$type <- factor(df$type, c("HSC", "CMP", "GMP", "MEP", "mixed", "PreBNK", "ref"))

myColorScale <- c('HSC_a'='#EF5350','HSC_b'='#CC0C00','CMP_a'='#FF7F0E','CMP_b'='#FFD147','CMP_c'='#CC9900',
                  'GMP_a'='#33CC00', 'GMP_b' ='#99CC00', 'MEP_a'='#8C5782', 'MEP_b'='#990080', 'MEP_c'='#BA68C8',
                  'MLP_a'='#46B8DA', 'MLP_b'='#0288D1', 'PreBNK'='#1A0099', 'ref'='grey')

tiff(paste0(OUT_DIR,"gsea_PCA_nonsig0_p0.05_NO_correlations_with_custom.tiff"), units="in", width=4, height=3.2, res=300)
ggplot(df,aes(x=pc1,y=pc2, group=type, label=sample_new_name, shape=category)) + theme_classic()+geom_point(aes(fill = sample, shape=category), color="black", size = 3.5, stroke = 1)+ labs(x="PC1 (35.03%)", y="PC2 (22.13%)")+
  geom_text_repel(size=2.5) + scale_fill_manual(values=myColorScale) + scale_color_manual(values=c("#CC0C00", "#FF7F0E", "#33CC00", "#990080", "#46B8DA", "#1A0099", "grey"))+ scale_shape_manual(values=c(21,23)) +
  theme(legend.position="none")
dev.off()
