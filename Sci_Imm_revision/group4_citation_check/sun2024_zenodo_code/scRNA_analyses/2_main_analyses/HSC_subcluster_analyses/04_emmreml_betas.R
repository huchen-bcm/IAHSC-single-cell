#R version 4.1.0
library(ggplot2)
library (DESeq2)
library(statmod)
library(RColorBrewer)
library(edgeR)
library(devtools)
library(cluster)
library(ggrepel)
library(EMMREML)
library(qvalue)


# We run emmreml to determine DE genes within each HSC subcluster 


# setup -------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/HSC_subcluster_analyses/") 
OUT_DIR <- "4_emmreml_betas/"
dir.create(OUT_DIR)


# Read in the filter 15 samples list of samples with fewer than 15 cells ---------------------------------------------
# These samples will not be included in downstream analyses
filter_15_list <- readRDS(file=paste0("3_filter15/filter_15_list.rds"))
names(filter_15_list) <- c("c0", "c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8", "c9")
clusters <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9)
DE_genes_count_pval <- vector(length=length(clusters))
DE_genes_count_qBH <- vector(length=length(clusters))
DE_genes_count_qST <- vector(length=length(clusters))
DE_0.1_genes_count_pval <- vector(length=length(clusters))
DE_0.1_genes_count_qBH <- vector(length=length(clusters))
DE_0.1_genes_count_qST <- vector(length=length(clusters))
num_genes <- vector(length=length(clusters))
threshold <- 0

expression_threshold=1
for(i in 1:length(clusters)){
  
  name <- clusters[i]
  {
    print("############")
    print(paste(i,name))
    print("############")
  }
  
  #### 1. Read and format pseudobulk data & meta_data
  pseudobulk <- readRDS(file = paste0("2_pseudobulk/",name,"_pseudobulk.rds"))
  meta_data <- read.csv(file="all_samples_meta_data_emmreml.csv", header=TRUE)
  meta_data$detailed_group=paste0(meta_data$group,"_",meta_data$timepoint)
  meta_data$group=factor(meta_data$group,levels=c("CTL","BCG"))
  
  #### 2. Make the orders of both objects match:
  rownames(meta_data)=meta_data$Sample
  meta_data=meta_data[order(
  factor(meta_data$group,levels=c("CTL","BCG")),
  meta_data$timepoint,
  factor(meta_data$donor,levels=c("S1","S2","S3","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20","S21"))),]
  pseudobulk=pseudobulk[,order(match(colnames(pseudobulk),meta_data$Sample))]
  
  
  ## Check orders of elements match
  print(paste("Number of incongruent elements before filteirng steps",length(which(colnames(pseudobulk)!=rownames(meta_data)))))
  
  #### 3. Filter out undesired samples:
  ## 3.1 Filter out individuals with fewer than 15 cells
  filter15_samples <- filter_15_list[[i]]
  samples_unfilt <- colnames(pseudobulk)[colnames(pseudobulk) %in% filter15_samples==FALSE]
  
  pseudobulk <- pseudobulk[colnames(pseudobulk) %in% samples_unfilt]
  meta_data=meta_data[which(meta_data$Sample %in% samples_unfilt),]
  
  ## 3.2 Filter out samples from unpaired individuals:
  donors_present_in_Td0=unique(meta_data$donor[which(meta_data$timepoint==0)])
  donors_present_in_Tm3=unique(meta_data$donor[which(meta_data$timepoint==1)])
  donors_present_in_both=donors_present_in_Td0[which(donors_present_in_Td0 %in% donors_present_in_Tm3)]
  samples_filt=meta_data$Sample[which(meta_data$donor %in% donors_present_in_both)]
  
  pseudobulk <- pseudobulk[colnames(pseudobulk) %in% samples_filt]
  meta_data=meta_data[which(meta_data$Sample %in% samples_filt),]
  
  ## Redefine donor factor after dropping levels:
  donors=c("S1","S2","S3","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20","S21")
  donors=donors[which(donors %in% unique(meta_data$donor))]
  meta_data$donor=factor(meta_data$donor,levels=donors)
  
  ## Check orders of elements match
  print(paste("Number of incongruent elements after filteirng steps",length(which(colnames(pseudobulk)!=rownames(meta_data)))))
  
  ## analyze only if there are more than 2 unique controls (= at least 4 total CTL samples)
  if(length(meta_data$group[which(meta_data$group == "CTL")])>=4){
    if(length(meta_data$group[which(meta_data$group == "BCG")])>=4){
      
      #### 4. Filter out undesired (lowly expressed) genes: I change the strategy of filtering, keeping genes above threshold in @ least 1 detailed group (CTL_0, CTL_1, BCG_0 or BCG_1)
      y <- DGEList(pseudobulk)
      y <- calcNormFactors(y)
      v <- voom(y, plot=TRUE,span=0.1)
      cpm <- v$E
      medians=cpm[,1:length(unique(meta_data$detailed_group))]
      colnames(medians)=unique(meta_data$detailed_group)
      for(j in 1:ncol(medians))
        medians[,j]=apply(cpm[,which(meta_data$detailed_group==colnames(medians)[j])],1,median)
      
      max_median_across_experimental_units=apply(medians,1,max)
      medians=cbind(medians,max_median_across_experimental_units)
      
      genes_to_keep=names(max_median_across_experimental_units)[which(max_median_across_experimental_units>expression_threshold)]
      pseudobulk=pseudobulk[genes_to_keep,]
      
      
      ## Check orders of elements match after filtering
      length(which(colnames(pseudobulk)==rownames(meta_data)))
      length(which(colnames(pseudobulk)!=rownames(meta_data)))
      
      ## Check numbers of people in each group
      print(summary(factor(meta_data$detailed_group)))
      
      ##  Now, model DE. First, using EMMREML (weights do not matter here, so the design neither does).
      y <- DGEList(pseudobulk)
      y <- calcNormFactors(y)
      v <- voom(y, plot=TRUE,span=0.1)
      exp <- v$E
      
      K=matrix(rep(0,length(levels(meta_data$donor))*length(levels(meta_data$donor))),nrow=length(levels(meta_data$donor)),ncol=length(levels(meta_data$donor)))
      row.names(K) <- levels(meta_data$donor)
      colnames(K) <- levels(meta_data$donor)
      
      for(j in 1:length(colnames(K)))
      {
        K[j,j] <- 1
      }
      
      ##Design with fixed effects
      design = model.matrix(~group+ timepoint + group:timepoint, data=meta_data)
      v <- voom(y,design,plot=TRUE)
      exp=v$E
      
      ##Set up res_full matrix with rownames matching "exp" rownames
      res_full=exp[,1:(4*ncol(design))]
      colnames(res_full)[1:ncol(design)]=paste0("beta_",colnames(design))
      colnames(res_full)[(ncol(design)+1):(2*ncol(design))]=paste0("sdev_",colnames(design))
      colnames(res_full)[((2*ncol(design))+1):(3*ncol(design))]=paste0("p_value_",colnames(design))
      colnames(res_full)[((3*ncol(design))+1):(4*ncol(design))]=paste0("q_BH_",colnames(design))
      
      random_effects=exp[,1:length(levels(meta_data$donor))]
      colnames(random_effects)=levels(meta_data$donor)
      
      ##Create Z matrix assigning donor id to each sample
      Z=matrix(rep(0,nrow(meta_data)*ncol(K)),nrow=nrow(meta_data),ncol=ncol(K))
      rownames(Z)=rownames(meta_data)
      colnames(Z)=colnames(K)
      for(j in 1:ncol(Z))
      {
        set=which(meta_data$donor == colnames(Z)[j])
        Z[set,j]=1
      }
      
      for(j in 1:nrow(exp))
      {
        emma=emmreml(y=exp[j,],X=design,Z=as.matrix(Z),K=as.matrix(K),varbetahat=T,varuhat=T,PEVuhat=T,test=T)
        res_full[j,]=t(c(emma$betahat,emma$varbetahat,emma$pvalbeta[,"none"],emma$pvalbeta[,"BH"]))
        random_effects[j,]=t(emma$uhat)
        if(j%%3000==0)print(j)
      }
      
      
      ##corrected q values
      res_full <- data.frame(res_full)
      p <- res_full[,"p_value_groupBCG.timepoint"]
      BH_corrected <- p.adjust(p, method = "BH", n = length(p))
      qvalues=qvalue(p)$qvalues
      res_full$BH_corrected <- BH_corrected
      res_full$ST_corrected <- qvalues
      
      write.csv(res_full, file=paste0(OUT_DIR, name, "_res_full.csv"))
      write.csv(random_effects, file=paste0(OUT_DIR,  name, "_random_effects.csv"))
      
      DE_genes_count_pval[i] <- length(row.names(subset(res_full, res_full[,"p_value_groupBCG.timepoint"] < 0.05)))
      DE_genes_count_qBH[i] <- length(row.names(subset(res_full, res_full[,"BH_corrected"] < 0.05)))
      DE_0.1_genes_count_pval[i] <- length(row.names(subset(res_full, res_full[,"p_value_groupBCG.timepoint"] < 0.1)))
      DE_0.1_genes_count_qBH[i] <- length(row.names(subset(res_full, res_full[,"BH_corrected"] < 0.1)))
      DE_genes_count_qST[i]=length(which(qvalues<0.05))
      DE_0.1_genes_count_qST[i]=length(which(qvalues<0.1))
      
      num_genes[i] <- length(row.names(res_full))
      
    }
    else{
      DE_genes_count_pval[i] <- "N/A"
      DE_genes_count_qBH[i] <- "N/A"
      DE_0.1_genes_count_pval[i] <- "N/A"
      DE_0.1_genes_count_qBH[i] <- "N/A"
      DE_genes_count_qST[i]="N/A"
      DE_0.1_genes_count_qST[i]="N/A"
      
      num_genes[i] <- "N/A"
    }
  }
  
  else{
    DE_genes_count_pval[i] <- "N/A"
    DE_genes_count_qBH[i] <- "N/A"
    DE_0.1_genes_count_pval[i] <- "N/A"
    DE_0.1_genes_count_qBH[i] <- "N/A"
    DE_genes_count_qST[i]="N/A"
    DE_0.1_genes_count_qST[i]="N/A"
    
    num_genes[i] <- "N/A"
  }
  
  print(i)
  
}

DE_genes_count <- cbind(num_genes, DE_genes_count_pval, DE_genes_count_qBH, DE_genes_count_qST, DE_0.1_genes_count_pval, DE_0.1_genes_count_qBH,DE_0.1_genes_count_qST)
row.names(DE_genes_count) <- clusters
colnames(DE_genes_count ) <- c("total_genes", "total_pval_0.05", "total_qBH_0.05", "total_qST_0.05","total_pval_0.1", ",total_qBH_0.1", "total_qST_0.1")
write.csv(DE_genes_count, file=paste(paste0(OUT_DIR, "DE_genes_count_table.csv")))


