
#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

library(ggplot2)
library (DESeq2)
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


# The purpose of this code is to make a bargraph that shows the number of significant genes within each HSPC cluster

setwd("/scRNA_analyses/2_main_analyses/4_mash_downstream/")
OUT_DIR <- c("DE_genes_bargraphs/")
dir.create(OUT_DIR)


# Read in MASH results without correlations with custom matrices ----------------------------------------------------

PM <- read.table(file=paste0("../MASH/mash_results/posteriorMeans_wcustom.txt"), header=TRUE)
lfsr <- read.table(file=paste0("../MASH/mash_results/lfsr_wcustom_output.txt"), header=TRUE)

#check that genes are in the same order
which(row.names(PM) != row.names(lfsr))



# Make bar graph of num DE genes per cluster lfsr < 0.01 ------------------------------

PM <- read.table(file=paste0("../MASH/mash_results/posteriorMeans_wcustom.txt"), header=TRUE)
lfsr <- read.table(file=paste0("../MASH/mash_results/lfsr_wcustom_output.txt"), header=TRUE)

#check that genes are in the same order
which(row.names(PM) != row.names(lfsr))


clusters <- colnames(PM)
new_names <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_a", "GMP_b", "MEP_a", "MEP_b", "MEP_c", "MLP_a", "MLP_b", "PreBNK")
for(i in 1:length(clusters)){
  direction_vec <- vector(length=length(row.names(clusters[i])))
  sig_vec <- vector(length=length(row.names(clusters[i])))
  combined_vec <- vector(length=length(row.names(clusters[i])))
  
  vec_PM <- PM[,paste0(clusters[i])]
  direction_vec[which(vec_PM > 0)] <- "up"
  direction_vec[which(vec_PM <= 0)] <- "down"
  
  vec_lfsr <- lfsr[,paste0(clusters[i])]
  sig_vec[which(vec_lfsr <= 0.01)] <- "sig"
  sig_vec[which(vec_lfsr > 0.01)] <- "not_sig"
  
  combined_vec <- paste0(direction_vec, "_", sig_vec)
  
  count <- c(length(which(combined_vec == "up_sig")), length(which(combined_vec == "down_sig")))
  assign(paste0(new_names[i],"_counts"), count)
}

count_list <- list(HSC_a_counts, HSC_b_counts, CMP_a_counts, CMP_b_counts, CMP_c_counts, GMP_a_counts, GMP_b_counts, MEP_a_counts, MEP_b_counts, MEP_c_counts, MLP_a_counts, MLP_b_counts, PreBNK_counts)
val <- vector()
clust <- vector()
for(i in 1:length(count_list)){
  val <- append(val,as.numeric(count_list[[i]]))
  clust <- append(clust, rep(new_names[i],2))
}

df_0.01 <- data.frame(
  values <- val,
  group <- c(rep(c("up", "down"),length(count_list))),
  cluster <- clust
)
colnames(df_0.01) <- c("values", "group", "cluster")
df_0.01$values2 <- ifelse(df_0.01$group == "down", -1 * df_0.01$values, df_0.01$values)

df_0.01$group <- factor(df_0.01$group, c("up", "down"))
df_0.01$cluster <- factor(df_0.01$cluster, c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_a", "GMP_b", "MEP_a", "MEP_b", "MEP_c", "MLP_a", "MLP_b", "PreBNK"))


myColorScale <- c('HSC_a'='#EF5350','HSC_b'='#CC0C00','CMP_a'='#FF7F0E','CMP_b'='#FFD147','CMP_c'='#CC9900',
                  'GMP_a'='#33CC00', 'GMP_b' ='#99CC00', 'MEP_a'='#8C5782', 'MEP_b'='#990080', 'MEP_c'='#BA68C8',
                  'MLP_a'='#46B8DA', 'MLP_b'='#0288D1', 'PreBNK'='#1A0099')


# save same plot as RDS -------------------------------------------------
p <- ggplot(df_0.01, aes(x = cluster, y = values2, fill = cluster)) + 
  geom_bar(stat = "identity", aes(fill=cluster), alpha=0.7, color="black") +
  scale_fill_manual(values = myColorScale)+theme_classic()+ 
  labs(main="DR genes in all clusters post MASH with custom", x="",y="Number of DR genes (lfsr<0.01)", fill="")+
  theme(axis.text.x = element_text(angle = 45, hjust=1, size=10))+ theme(axis.text.y = element_text(angle = 45, size=10)) +
  theme(axis.title.y = element_text( size=11)) + guides(color=FALSE) + scale_y_continuous(labels=abs)+
  theme(legend.position="none")

saveRDS(p, file=paste0(OUT_DIR,"Figure_2A_bargraph.rds"))



# number of sig genes overlapping between HSC c1 and HSC c2 -----------------------------------------

PM <- read.table(file=paste0("../MASH/mash_results/posteriorMeans_wcustom.txt"), header=TRUE)
lfsr <- read.table(file=paste0("../MASH/mash_results/lfsr_wcustom_output.txt"), header=TRUE)
#check that genes are in the same order
which(row.names(PM) != row.names(lfsr))

row.names(lfsr)[which(lfsr[,'HSC_a'] <0.01)]
row.names(lfsr)[which(lfsr[,'HSC_b'] <0.01)]

length(intersect(row.names(lfsr)[which(lfsr[,'HSC_a'] <0.01)], row.names(lfsr)[which(lfsr[,'HSC_b'] <0.01)]))
length(unique(c(row.names(lfsr)[which(lfsr[,'HSC_a'] <0.01)], row.names(lfsr)[which(lfsr[,'HSC_b'] <0.01)])))
