#R version 4.1.0
install.packages("VennDiagram")
library(ggplot2)
library(cowplot)
library(reticulate)
library(tidyr)
library(reshape2)
library(plyr)
library(dplyr)
library(grid)
library(gridExtra)
library(Seurat)
library(sctransform)
library(data.table)
library(VennDiagram)

# This code generates a table for each cell type with log2FC (D90-D0) expression for each gene in each sample 
# We correlate FC cytokine secretion with these log2FC values to determine genes who expression in HSCs correlates with cytokine secretion by PBMCs 


# setup ---------------------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/4_mash_downstream/")
OUT_DIR <- "DE_genes_cytokines_correlation/"
dir.create(OUT_DIR)


# Read in cytokine data ---------------------------------------------------

cytokine_data <- data.frame(fread("../PBMC_cytokine_data_CHM.txt"))


# Read in the corrected expression matrices -------------------------------
# create a list of expression tables of genes with lfsr < 0.1 ----------------------

celltypes <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_b1", "GMP_b2", "MEP_a1", "MEP_a2", "MEP_a3", "mixed_a", "mixed_b", "PreBNK")
data_genes <- vector()

for(i in 1:length(celltypes)){
  c <- celltypes[i]
  data <- read.table(file=paste0("../emmreml_betas/raw_expression/", c, "_raw_expression.txt"))
  data$X <- row.names(data)
  
  # read in mash results 
  PM <- read.table(file=paste0("../MASH/mash_results/posteriorMeans_wcustom.txt"), header=TRUE)
  lfsr <- read.table(file=paste0("../MASH/mash_results/lfsr_wcustom_output.txt"), header=TRUE)
  
  DE_df <- data.frame(cbind(PM[,paste0(celltypes[i])], lfsr[,paste0(celltypes[i])]))
  row.names(DE_df) <- row.names(PM)
  colnames(DE_df) <- c("beta", "lfsr")

  # get genes with lfsr < 0.1
  DE_subset_genes <- row.names(subset(DE_df, DE_df$lfsr < 0.1))
  
  data_subset <- subset(data, data$X %in% DE_subset_genes)
  
  assign(paste0(c), data_subset)
}

expression_data <- list(HSC_a, HSC_b, CMP_a, CMP_b, CMP_c, GMP_b1, GMP_b2, MEP_a1, MEP_a2, MEP_a3, mixed_a, mixed_b, PreBNK)
expression_data_names <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_b1", "GMP_b2", "MEP_a1", "MEP_a2", "MEP_a3", "mixed_a", "mixed_b","PreBNK")

# read in meta data ---------------------------------------------------
meta_data <- read.csv(file="../all_samples_meta_data.csv")
meta_data$Sample <- gsub('-','\\.', meta_data$Sample)


# Loop through each cell type to create a matrix of log2FC (D90-D0) values for each gene and donor --------------------------------

Tm3_all <- c("SNG.LB.SS.1S.RS.S1.CD34neg_S1_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S2.CD34neg_S2_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S5.CD34neg_S4_R1_001_Tm3_BCG",
             "SNG.LB.SS.1S.RS.S6.CD34neg_S5_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S8.CD34neg_S7_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S9.CD34neg_S8_R1_001_Tm3_BCG",
             "SNG.LB.SS.1S.RS.S10.CD34neg_S9_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S12.CD34neg_S11_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S14.CD34neg_S13_R1_001_Tm3_BCG",
             "SNG.LB.SS.1S.RS.S15.CD34neg_S14_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S16.CD34neg_S15_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S17.CD34neg_S16_R1_001_Tm3_BCG",
             "SNG.LB.SS.1S.RS.S19.CD34neg_S18_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S20.CD34neg_S19_R1_001_Tm3_BCG", "SNG.LB.SS.1S.RS.S21.CD34neg_S20_R1_001_Tm3_BCG", 
             "SNG.LB.SS.1S.RS.S3.CD34neg_S3_R1_001_Tm3_CTL", "SNG.LB.SS.1S.RS.S7.CD34neg_S6_R1_001_Tm3_CTL", "SNG.LB.SS.1S.RS.S11.CD34neg_S10_R1_001_Tm3_CTL",
             "SNG.LB.SS.1S.RS.S13.CD34neg_S12_R1_001_Tm3_CTL", "SNG.LB.SS.1S.RS.S18.CD34neg_S17_R1_001_Tm3_CTL")


Td0_all <- c("SNG.LB.SS.1S.RS.S1.CD34neg_S1_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S2.CD34neg_S2_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S5.CD34neg_S4_R1_001_Td0_BCG",
             "SNG.LB.SS.1S.RS.S6.CD34neg_S5_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S8.CD34neg_S7_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S9.CD34neg_S8_R1_001_Td0_BCG",
             "SNG.LB.SS.1S.RS.S10.CD34neg_S9_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S12.CD34neg_S11_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S14.CD34neg_S13_R1_001_Td0_BCG",
             "SNG.LB.SS.1S.RS.S15.CD34neg_S14_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S16.CD34neg_S15_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S17.CD34neg_S16_R1_001_Td0_BCG",
             "SNG.LB.SS.1S.RS.S19.CD34neg_S18_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S20.CD34neg_S19_R1_001_Td0_BCG", "SNG.LB.SS.1S.RS.S21.CD34neg_S20_R1_001_Td0_BCG", 
             "SNG.LB.SS.1S.RS.S3.CD34neg_S3_R1_001_Td0_CTL", "SNG.LB.SS.1S.RS.S7.CD34neg_S6_R1_001_Td0_CTL", "SNG.LB.SS.1S.RS.S11.CD34neg_S10_R1_001_Td0_CTL",
             "SNG.LB.SS.1S.RS.S13.CD34neg_S12_R1_001_Td0_CTL", "SNG.LB.SS.1S.RS.S18.CD34neg_S17_R1_001_Td0_CTL")

logfc_mat_final_list <- list()
for(k in 1:length(expression_data)){
  #split each corrected expression matrix into 2 parts Td0 and Tm3
  celltype <- celltypes[k]
  mat <- expression_data[[k]]
  Td0_cols <- which(grepl("Td0", colnames(mat))==TRUE)
  Td0_subset <- mat[Td0_cols]
  row.names(Td0_subset) <- mat$X
  
  Tm3_cols <- which(grepl("Tm3", colnames(mat))==TRUE)
  Tm3_subset <- mat[Tm3_cols]
  row.names(Tm3_subset) <- mat$X
  
  #for each Tm3 sample get the corresponsing Td0 sample and subtract columns 
  
  logfc_mat <- data.frame()
  count <- 0
  name_short <- vector()
  name_count <- 1
  for(i in 1:length(Tm3_all))
  {
    name <- Tm3_all[i]
    if(name %in% colnames(Tm3_subset))
    {
      d <- name
      Td0_pair_name <- Td0_all[which(Tm3_all == d)]
      Td0_col <- Td0_subset[Td0_pair_name]
      Tm3_col <- Tm3_subset[d]
      log2fc <- Tm3_col - Td0_col
      name_short[name_count] <- subset(meta_data, meta_data$Sample %in% as.character(d))$donor
      name_count <- name_count+1
      if(count==0)
      {
        logfc_mat <- log2fc[1]
      }
      else
      {
        logfc_mat <- cbind(logfc_mat, log2fc[1])
      }
      
      count <- count + 1
    }
    
  }
  colnames(logfc_mat) <- name_short
  logfc_mat_final <- t(logfc_mat)
  logfc_mat_final_list[[k]] <- logfc_mat_final
  names(logfc_mat_final_list)[k] <- expression_data_names[k]
  saveRDS(logfc_mat_final_list, file=paste0(OUT_DIR, "logfc_mat_final_list_DR_only.rds"))

  

  # record fold change cytokine secretion for each cytokine in matched samples -------------------------------------------------
  # for each cytokine determine the number of genes whose log2FC expression correlates with FC cytokine secretion ---------------

  il10 <- vector(length=length(name_short))
  il1b <- vector(length=length(name_short))
  il1ra <- vector(length=length(name_short))
  il6 <- vector(length=length(name_short))
  ifna <- vector(length=length(name_short))
  ifng <- vector(length=length(name_short))
  tnfa <- vector(length=length(name_short))
   
  for(i in 1:length(name_short))
  {
    cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == name_short[i])
    cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
    il10[i] <- as.numeric(cyt_subset2$IL10_FC)
    il1b[i] <- as.numeric(cyt_subset2$IL1B_FC)
    il1ra[i] <- as.numeric(cyt_subset2$IL1RA_FC)
    il6[i] <- as.numeric(cyt_subset2$IL6_FC)
    ifna[i] <- as.numeric(cyt_subset2$IFNA_FC)
    ifng[i] <- as.numeric(cyt_subset2$IFNG_FC)
    tnfa[i] <- as.numeric(cyt_subset2$TNF_FC)
    print(i)
  }
  
  cyt_mat <- cbind(il10, il1b, il1ra, il6, ifna, ifng, tnfa)
  row.names(cyt_mat) <- name_short

  
  count <- 1
  gene_name <- vector()
  cyt_name <- vector()
  pval <- vector()
  rho <- vector()
  for(i in 1:length(colnames(logfc_mat_final)))
  {
    gene_vals <- logfc_mat_final[,i]
    for(j in 1:length(colnames(cyt_mat)))
    {
      gene_name[count] <- colnames(logfc_mat_final)[i]
      cyt_vals <- cyt_mat[,j]
      cyt_name[count] <- colnames(cyt_mat)[j]
      val <- cor.test(cyt_vals, gene_vals, method="spearman")
      pval[count] <- val$p.value
      rho[count] <- val$estimate
      
      count <- count + 1
    }
  }
  
  df <- data.frame(cbind(gene_name, cyt_name, pval, rho))
  df_list <- split(df, f = df$cyt_name)
  for(j in 1:length(df_list))
  {
    cyt_name <- names(df_list)[j]
    dat_cyt <- df_list[[j]]
    padj <- p.adjust(dat_cyt[,3], method="BH")
    dat_cyt <- cbind(dat_cyt, padj)
    assign(paste0(celltype, "_",cyt_name, "_result"), dat_cyt)
  }
  print(k)
}

ifna_result_list <- list(HSC_a_ifna_result, HSC_b_ifna_result, CMP_a_ifna_result, CMP_b_ifna_result, CMP_c_ifna_result, GMP_b1_ifna_result, GMP_b2_ifna_result, MEP_a1_ifna_result, MEP_a2_ifna_result, MEP_a3_ifna_result, mixed_a_ifna_result, mixed_b_ifna_result, PreBNK_ifna_result)
saveRDS(ifna_result_list, file=paste0(OUT_DIR, "ifna_result_list_DE_genes_only.rds"))
ifng_result_list <- list(HSC_a_ifng_result, HSC_b_ifng_result, CMP_a_ifng_result, CMP_b_ifng_result, CMP_c_ifng_result, GMP_b1_ifng_result, GMP_b2_ifng_result, MEP_a1_ifng_result, MEP_a2_ifng_result, MEP_a3_ifng_result, mixed_a_ifng_result, mixed_b_ifng_result, PreBNK_ifng_result)
saveRDS(ifng_result_list, file=paste0(OUT_DIR, "ifng_result_list_DE_genes_only.rds"))
il10_result_list <- list(HSC_a_il10_result, HSC_b_il10_result, CMP_a_il10_result, CMP_b_il10_result, CMP_c_il10_result, GMP_b1_il10_result, GMP_b2_il10_result, MEP_a1_il10_result, MEP_a2_il10_result, MEP_a3_il10_result, mixed_a_il10_result, mixed_b_il10_result, PreBNK_il10_result)
saveRDS(il10_result_list, file=paste0(OUT_DIR, "il10_result_list_DE_genes_only.rds"))
il1b_result_list <- list(HSC_a_il1b_result, HSC_b_il1b_result, CMP_a_il1b_result, CMP_b_il1b_result, CMP_c_il1b_result, GMP_b1_il1b_result, GMP_b2_il1b_result, MEP_a1_il1b_result, MEP_a2_il1b_result, MEP_a3_il1b_result, mixed_a_il1b_result, mixed_b_il1b_result, PreBNK_il1b_result)
saveRDS(il1b_result_list, file=paste0(OUT_DIR, "il1b_result_list_DE_genes_only.rds"))
il1ra_result_list <- list(HSC_a_il1ra_result, HSC_b_il1ra_result, CMP_a_il1ra_result, CMP_b_il1ra_result, CMP_c_il1ra_result, GMP_b1_il1ra_result, GMP_b2_il1ra_result, MEP_a1_il1ra_result, MEP_a2_il1ra_result, MEP_a3_il1ra_result, mixed_a_il1ra_result, mixed_b_il1ra_result, PreBNK_il1ra_result)
saveRDS(il1ra_result_list, file=paste0(OUT_DIR, "il1ra_result_list_DE_genes_only.rds"))
il6_result_list <- list(HSC_a_il6_result, HSC_b_il6_result, CMP_a_il6_result, CMP_b_il6_result, CMP_c_il6_result, GMP_b1_il6_result, GMP_b2_il6_result, MEP_a1_il6_result, MEP_a2_il6_result, MEP_a3_il6_result, mixed_a_il6_result, mixed_b_il6_result, PreBNK_il6_result)
saveRDS(il6_result_list, file=paste0(OUT_DIR, "il6_result_list_DE_genes_only.rds"))
tnfa_result_list <- list(HSC_a_tnfa_result, HSC_b_tnfa_result, CMP_a_tnfa_result, CMP_b_tnfa_result, CMP_c_tnfa_result, GMP_b1_tnfa_result, GMP_b2_tnfa_result, MEP_a1_tnfa_result, MEP_a2_tnfa_result, MEP_a3_tnfa_result, mixed_a_tnfa_result, mixed_b_tnfa_result, PreBNK_tnfa_result)
saveRDS(tnfa_result_list, file=paste0(OUT_DIR, "tnfa_result_list_DE_genes_only.rds"))



# Bargraph summarizing all results ----------------------------------------

celltypes <- c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_b1", "GMP_b2", "MEP_a1", "MEP_a2", "MEP_a3", "mixed_a", "mixed_b", "PreBNK")
cytokines <- c("ifna", "ifng", "il10", "il1b", "il1ra", "il6", "tnfa")


##ifna
num_genes_ifna <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(ifna_result_list[[j]])
  num_genes_ifna[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##ifng
num_genes_ifng <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(ifng_result_list[[j]])
  num_genes_ifng[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##il10
num_genes_il10 <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(il10_result_list[[j]])
  num_genes_il10[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##il1b
num_genes_il1b <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(il1b_result_list[[j]])
  num_genes_il1b[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##il1ra
num_genes_il1ra <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(il1ra_result_list[[j]])
  num_genes_il1ra[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##il6
num_genes_il6 <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(il6_result_list[[j]])
  num_genes_il6[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

##tnfa
num_genes_tnfa <- vector(length=length(celltypes))
for(j in 1:length(celltypes))
{
  dat <- data.frame(tnfa_result_list[[j]])
  num_genes_tnfa[j] <- length(row.names(subset(dat, dat$padj < 0.1)))
}

df_summary <- data.frame(
  values <- c(num_genes_ifna, num_genes_ifng, num_genes_il10, num_genes_il1b, num_genes_il1ra, num_genes_il6, num_genes_tnfa),
  cytokine_name <- c(rep("ifna", length(celltypes)), rep("ifng", length(celltypes)), rep("il10", length(celltypes)), rep("il1b", length(celltypes)), rep("il1ra", length(celltypes)), rep("il6", length(celltypes)), rep("tnfa", length(celltypes))),
  celltype_name <- c(rep(c("HSC c1", "HSC c2", "CMP c1", "CMP c2", "CMP c3", "GMP c1", "GMP c2", "MEP c1", "MEP c2", "MEP c3", "MLP c1", "MLP c2", "PreBNK"),7))
)
colnames(df_summary) <- c("values", "cytokine_name", "celltype_name")
df_summary$celltype_name <- factor(df_summary$celltype_name, c("HSC c1", "HSC c2", "CMP c1", "CMP c2", "CMP c3", "GMP c1", "GMP c2", "MEP c1", "MEP c2", "MEP c3", "MLP c1", "MLP c2", "PreBNK"))

#tiff(paste0(OUT_DIR, "summary_df_DE_genes_only.tiff"), units="in", width=4, height=3, res=350)
p <- ggplot(data=df_summary, aes(x=celltype_name, y=values, fill=cytokine_name)) +
  geom_bar(stat="identity", size=0.5)+theme_linedraw()+ylab("Correlated genes (padj < 0.1)")+xlab("")+ labs(fill="")+
  theme(axis.text.x = element_text( angle=45, hjust=1, size=10)) +  theme(axis.title.y = element_text(size=12)) + scale_fill_manual(values=c("#00468B", "#ED0000", "#42B540","#0099B4", "#925E9F", "#FDAF91", "#AD002A")) +
  theme (legend.key.size = unit (0.5, 'cm'))
#dev.off()

final_file <- file.path(OUT_DIR,paste0("summary_df_DE_genes_only.pdf"))
ggsave(filename = final_file, plot = p, height = 3, width = 4)






# General code to make scatterplots for individual gene-cytokine pairs -----------------------------------
# enter celltype, gene, and cytokine values where indicated -----------------------------------------------

celltype <- "HSC_a" ###CHOOSE celltype
mat <- expression_data[[1]] ##CHOOSE data
Td0_cols <- which(grepl("Td0", colnames(mat))==TRUE)
Td0_subset <- mat[Td0_cols]
row.names(Td0_subset) <- mat$X

Tm3_cols <- which(grepl("Tm3", colnames(mat))==TRUE)
Tm3_subset <- mat[Tm3_cols]
row.names(Tm3_subset) <- mat$X

#for each Tm3 sample get the corresponsing Td0 sample and subtract columns 

logfc_mat <- data.frame()
count <- 0
name_short <- vector()
name_count <- 1
for(i in 1:length(Tm3_all)){
  name <- Tm3_all[i]
  if(name %in% colnames(Tm3_subset))
  {
    d <- name
    Td0_pair_name <- Td0_all[which(Tm3_all == d)]
    Td0_col <- Td0_subset[Td0_pair_name]
    Tm3_col <- Tm3_subset[d]
    log2fc <- Tm3_col - Td0_col
    name_short[name_count] <- subset(meta_data, meta_data$Sample %in% as.character(d))$donor
    name_count <- name_count+1
    if(count==0)
    {
      logfc_mat <- log2fc[1]
    }
    else
    {
      logfc_mat <- cbind(logfc_mat, log2fc[1])
    }
    
    count <- count + 1
  }
  
}
colnames(logfc_mat) <- name_short
logfc_mat_final <- t(logfc_mat)


gene_vec <- logfc_mat_final[,which(colnames(logfc_mat_final) == "FOSB")] ##CHOOSE GENE NAME

cytokine <- vector(length=length(name_short))
treatment <- vector(length=length(name_short))
meta_data <- read.csv(file="../all_samples_meta_data.csv")
meta_data$Sample <- gsub('-','\\.', meta_data$Sample)
for(i in 1:length(name_short)){
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == name_short[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  cytokine[i] <- as.numeric(cyt_subset2$IL1B_FC)  ##CHOOSE CYTOKINE
  
  #GET BCG/CTL info for each donor in name_short
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == name_short[i])])
  print(i)
}

val <- cor.test(cytokine, gene_vec, method="spearman")

df <- data.frame(
  g <- gene_vec,
  c <- cytokine,
  t <- treatment
)
colnames(df) <- c("g", "c", "Group")

#tiff(paste0(OUT_DIR, "HSC_a_KLF6_IL1B.tiff"), units="in", width=4, height=3.2, res=350)
p <- ggplot(df, aes(x=g, y=c)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_linedraw() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 0, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("Log2FC (Tm3 vs Td0)") +
  ylab("Fold Change IL1B") + ggtitle("FOSB in HSC c1")+
  scale_fill_manual(values=c("#D91B60", "grey")) + ##EDIT GRAPH LABS
  theme(axis.title = element_text(size=14)) + theme(plot.title=element_text(size=15))
#dev.off()

final_file <- file.path(OUT_DIR,paste0("HSC_a_FOSB_IL1B.pdf")) ##change the file name accordingly 
ggsave(filename = final_file, plot = p, height = 3.2, width = 4)






# correlation of IL1B and IL6 production by PBMCs -------------------------------------------------
cytokine_data <- data.frame(fread("PBMC_cytokine_data_CHM.txt"))
meta_data <- read.csv(file="all_samples_meta_data.csv")
donors <- unique(cytokine_data$Donor)
il1b <- vector(length=length(donors))
il6 <- vector(length=length(donors))
treatment <- vector(length=length(donors))

for(i in 1:length(donors))
{
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) == donors[i])
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  il1b[i] <- as.numeric(cyt_subset2$IL1B_FC)
  il6[i] <- as.numeric(cyt_subset2$IL6_FC)
  treatment[i] <- unique(meta_data$vaccination[which(meta_data$donor == donors[i])])
  print(i)
}

val <- cor.test(il1b, il6, method="spearman")
pval <- val$p.value

df <- data.frame(
  IL1B <- il1b,
  IL6 <- il6,
  t <- treatment,
  lab <- donors
)
row.names(df) <- donors
colnames(df) <- c("IL1B", "IL6", "Group", "Label")

tiff(paste0(OUT_DIR, "IL6_IL1B_corr.tiff"), units="in", width=5, height=4, res=250)
ggplot(df, aes(x=IL1B, y=IL6)) + geom_point(shape=21, size=4,aes(fill=Group)) + theme_test() + geom_smooth(method=lm, color="black", alpha=0.2, fullrange=TRUE, size=1)+
  geom_vline(xintercept = 1, color="grey", linetype="dashed") + geom_hline(yintercept = 1, color="grey", linetype="dashed") + xlab("Fold Change IL1B production of PBMCs") +
  ylab("Fold Change IL6 production of PBMCs") + ggtitle("IL1B vs. IL6")+
  annotate("text", x=0.9, y=2.9, label= paste0("Spearman corr = ",round(val$estimate,2),"\npvalue = ",round(pval,2))) + scale_fill_manual(values=c("#D91B60", "#1E88E5")) +
  geom_text_repel(aes(label=Label))
dev.off()









