#R version 4.1.0
#install.packages("EMMREML")
library(ggplot2)
library(statmod)
library(RColorBrewer)
library(edgeR)
library(devtools)
library(cluster)
library(ggrepel)
library(EMMREML)
library(qvalue)
library(dplyr)


# setup -------------------------------------------------------------------
setwd("/scATAC_analyses/2_main_analyses")
OUT_DIR <- paste0("0X_test_filter_thresholds/")
dir.create(OUT_DIR)

# Import pseudobulk objects -----------------------------------------------
filter_20_list <- readRDS(file=paste0("filter20/filter_20_list.rds"))
filter_20_index <- c(1,3,4,5,6,7,9,10,12,14,15,17,18,21,22,23,24)
##main cluster only
clusters <- c("C1", "C3", "C4", "C5", "C6", "C7","C9", "C10", "C12", "C14", "C15", "C17", "C18",  "C21", "C22","C23", "C24")
clusters_renamed <- c("unknown1", "CMP1", "CMP2", "HSC1", "HSC2", "MEP1","Pre-BNK", "GMP1","CLP", "MEP2", "MEP3", "CMP3", "MEP4", "CMP4", "CMP5", "GMP2", "GMP3")

threshold <- 0


expression_threshold <- 1

for(i in 1:length(clusters)){
  
  name <- clusters[i]
  new_name <- clusters_renamed[i]
  {
    print("############")
    print(paste(i,new_name))
    print("############")
  }
  
  plot_list <- list()
  options <- c(1, 1.25, 1.5, 1.75, 2, 2.25, 2.5)
  for(k in 1:length(options)){
    expression_threshold <- options[k]
    
    #### 1. Read and format pseudobulk data & meta_data
    pseudobulk <- readRDS(file = paste0("pseudobulk_clusters_harmony_all_combined/",name,"_pseudobulk.rds"))
    meta_data <- read.csv(file="meta_data_limma_betas_interaction.csv", header=TRUE)
    meta_data$detailed_group=paste0(meta_data$group,"_",meta_data$timepoint)
    meta_data$group=factor(meta_data$group,levels=c("CTL","BCG"))
    
    #### 2. Make the orders of both objects match:
    rownames(meta_data)=meta_data$sample
    meta_data=meta_data[order(
      factor(meta_data$group,levels=c("CTL","BCG")),
      meta_data$timepoint,
      factor(meta_data$donor,levels=c("S1","S2","S3","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20","S21"))),]
    pseudobulk=pseudobulk[,order(match(colnames(pseudobulk),meta_data$sample))]
    
    ## Check orders of elements match
    print(paste("Number of incongruent elements before filteirng steps",length(which(colnames(pseudobulk)!=rownames(meta_data)))))
    
    
    #### 3. Filter out undesired samples:
    ## 3.1 Filter out individuals with fewer than 20 cells
    filter20_samples <- filter_20_list[[filter_20_index[i]]]
    samples_unfilt <- colnames(pseudobulk)[colnames(pseudobulk) %in% filter20_samples==FALSE]
    
    pseudobulk <- pseudobulk[colnames(pseudobulk) %in% samples_unfilt]
    meta_data=meta_data[which(meta_data$sample %in% samples_unfilt),]
    
    ## 3.2 Filter out samples from unpaired individuals:
    donors_present_in_Td0=unique(meta_data$donor[which(meta_data$timepoint==0)])
    donors_present_in_Tm3=unique(meta_data$donor[which(meta_data$timepoint==1)])
    donors_present_in_both=donors_present_in_Td0[which(donors_present_in_Td0 %in% donors_present_in_Tm3)]
    samples_filt=meta_data$sample[which(meta_data$donor %in% donors_present_in_both)]
    
    pseudobulk <- pseudobulk[colnames(pseudobulk) %in% samples_filt]
    meta_data=meta_data[which(meta_data$sample %in% samples_filt),]
    
    ## Redefine donor factor after dropping levels:
    donors=c("S1","S2","S3","S5","S6","S7","S8","S9","S10","S11","S12","S13","S14","S15","S16","S17","S18","S19","S20","S21")
    donors=donors[which(donors %in% unique(meta_data$donor))]
    meta_data$donor=factor(meta_data$donor,levels=donors)
    
    ## Check orders of elements match
    print(paste("Number of incongruent elements before filteirng steps",length(which(colnames(pseudobulk)!=rownames(meta_data)))))
    
    write.csv(meta_data, file=paste0(OUT_DIR,"meta_data/", name, "_meta_data_options",options[k], ".csv"))
    
    ## analyze only if there are more than 2 unique controls (= at least 4 total CTL samples)
    if(length(meta_data$group[which(meta_data$group == "CTL")])>=2){
      if(length(meta_data$group[which(meta_data$group == "BCG")])>=2){
        
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
        res_full=data.frame(matrix(,nrow=length(row.names(exp)),ncol=length(1:(4*ncol(design)))))
        row.names(res_full) <- row.names(exp)
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
        
        ### plot betas ###
        expression_table <- exp
        betas <- data.frame(res_full)
        betas$genes <- row.names(betas)
        name <- clusters[i]
        
        expression_avg <- data.frame(gene = row.names(expression_table), expression = rowMeans(expression_table))
        expression_avg_binned <- expression_avg %>% mutate(quantile = ntile(expression, 20))
        colnames(expression_avg_binned)[3] <- "binned_expression"
        expression_avg_binned$binned_expression <- as.factor(expression_avg_binned$binned_expression)
        
        betas_subset <- subset(betas, select = c("genes", "beta_groupBCG.timepoint","BH_corrected")); colnames(betas_subset)[2] <- "betas"
        
        order <- match(as.character(expression_avg_binned$gene),as.character(betas_subset$genes))
        betas_subset <- betas_subset[order, ]
        which(row.names(betas_subset) != expression_avg_binned$gene)
        expression_avg_binned$betas <- betas_subset$betas
        
        print("plot for binned betas")
        ##Plot
        p <- ggplot(expression_avg_binned, aes(binned_expression, betas)) +
          geom_boxplot() +
          geom_point(size = 0.5) + 
          theme_bw() +
          geom_hline(yintercept = 0) +
          #ylim(c(-0.5,0.5)) +
          theme(panel.border = element_rect(colour = "black", fill = NA, size = 1)) +
          ggtitle(paste0(name,"_binned_betas_plot_option_",options[k]))
        plot_list[[k]] <- p
        
        
        
        write.csv(res_full, file=paste0(OUT_DIR,"res_full/", name, "_res_full_option", options[k], ".csv"))
        write.csv(random_effects, file=paste0(OUT_DIR,"random_effects/",  name, "_random_effects_option",options[k], ".csv"))
        

      }
 
    }
 
  }
  
  if(length(plot_list)>0){
    for(s in 1:length(plot_list)) {
      print("step1")
      file_name = paste0(OUT_DIR,"betas_plots/", clusters_renamed[i], "_binned_betas_plot_option", options[s],".tiff")
      print("step2")
      tiff(file_name)
      print(plot_list[[s]])
      dev.off()
      
    }
  }
  
  
}


