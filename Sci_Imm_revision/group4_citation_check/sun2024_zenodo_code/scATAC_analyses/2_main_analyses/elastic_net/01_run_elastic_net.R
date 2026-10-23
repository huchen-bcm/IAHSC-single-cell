##load libraries
library(glmnet)
library(dplyr)
library(stringi)
library(stringr)
library(readxl)
library(data.table)

set.seed(831)  # Set seed for reproducibility

## Datatype (set to RNA or ATAC)
data <- "ATAC"


## Set directory structure
setwd("/scATAC_analyses/2_main_analyses/elastic_net/")

cluster_names <- c("CMP1", "CMP2", "PreBNK", "GMP1", "CMP3",  "CMP4", "CMP5",  "GMP3", "GMP2")

for(j in 1:length(cluster_names)){
  
  print(cluster_names[j])
  
  out_dir <- paste0(cluster_names[j],"_EN_outputs/")
  dir.create(out_dir)
  
  #create output directories
  dir.create(file.path(paste0(out_dir,"predicted")))
  dir.create(file.path(paste0(out_dir,"weights")))
  pred_out = file.path(paste0(out_dir,"predicted"))
  weights_out = file.path(paste0(out_dir,"weights"))
  
  #Load in ATAC meta data
  meta_data <- read.csv(file="../meta_data_limma_betas_interaction.csv")
  
  #Load in a table of RNA or ATAC logFC values for each gene/peak (columns are donors, rows are genes/peaks)
  logfc_mat_final_list <- readRDS(file=paste0("../emmreml_downstream/DApeaks_cytokines_correlation/0.05_logfc_mat_final_list_DR_only.rds"))
  #choose the celltype
  cts <- logfc_mat_final_list[[cluster_names[j]]]
  cts <- t(cts)
  colnames(cts) <- paste0(colnames(cts))
  
  #Load in cytokine data and format
  cytokine_data <- data.frame(fread("../PBMC_cytokine_data_CHM.txt"))
  cyt_subset1 <- subset(cytokine_data, as.character(cytokine_data$Donor) %in% colnames(cts))
  cyt_subset2 <- subset(cyt_subset1, as.character(cyt_subset1$Timepoint) == "D90")
  #final cytokine scores
  cyt_FC_final <- cyt_subset2[,c("Donor", "Condition", "IL10_FC", "IL1B_FC", "IL1RA_FC", "IL6_FC", "IFNA_FC", "IFNG_FC", "TNF_FC")]
  write.table(cyt_FC_final, file=paste0(cluster_names[j],"_cyt_FC_final.txt"))
  
  cts <- cts[,cyt_FC_final$Donor]
  
  ## set alpha parameters to try
  alpha=c(0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0)
  NF <- dim(cts)[2]-1 # (number of samples) - 1
  print(paste0("number of samples-1 = ", NF))
  models <- c("IL10_FC_score", "IL1B_FC_score", "IL1RA_FC_score", "IL6_FC_score", "IFNA_FC_score", "IFNG_FC_score", "TNF_FC_score")
  
  
  ################################
  ## Elastic net model building ##
  ################################
  
  for (m in 1:7){
    print(models[m])
    
    for(i in 1:length(alpha)){
      
      print(alpha[i])
      
      #create table to store results
      pred_samp <- data.frame(samp=character(0),ind=character(0),predicted=numeric(0), weight=numeric(0), stringsAsFactors = FALSE)
      
      for(SAMP in 1:dim(cts)[2]){
        
        #################
        ## qqnorm data ##
        #################
        
        ## make sure that cytokine donors and cts match (and remove any samples there is not count information for)
        colnames(cts)==cyt_FC_final$Donor ## make sure == TRUE
        
        #Using alpha
        #Quantile normalize gene exp by column (sample)
        norm_counts<-as.matrix(apply(cts,2,function(x){return(qqnorm(x,plot=F)$x)}))
        
        #Remove test subject(s)
        #SAMP indexes from 1 to n samples
        norm_train<-norm_counts[,-SAMP]
        norm_test<-norm_counts[,SAMP]
        
        #Quantile normalize training samples by row (features)
        trainreads_norm<-as.matrix(apply(norm_train,1,function(x){return(qqnorm(x,plot=F)$x)}))
        
        if (m==1){
          train_model <- cyt_FC_final$IL10_FC[-SAMP]
        }else if (m==2){
          train_model <- cyt_FC_final$IL1B_FC[-SAMP]
        }else if (m==3){
          train_model <- cyt_FC_final$IL1RA_FC[-SAMP]
        }else if (m==4){
          train_model <- cyt_FC_final$IL6_FC[-SAMP]
        }else if (m==5){
          train_model <- cyt_FC_final$IFNA_FC[-SAMP]
        }else if (m==6){
          train_model <- cyt_FC_final$IFNG_FC[-SAMP]
        }else if (m==7){
          train_model <- cyt_FC_final$TNF_FC[-SAMP]
        }
        
        #QQ normalize each row in the test sample
        #Note that this could be much more efficient using parallelR (e.g. parSapply)
        
        #Store a new vector for normalizing the test sample
        testreads_normalized<-norm_test
        
        #For each feature
        for (d in 1:length(testreads_normalized)){
          #Define the ECDF (depending on training sample size, this can be replaced with simply the training data -- norm_train)
          a<-ecdf(norm_counts[d,])
          #From this ECDF, return the probability of values being less than the training sample
          probs<-a(norm_test[d])
          #To avoid extreme values outside the range of the training samples, give 0's and 1's a manageable quantile (e.g. 0.99 and .01)
          #Note depending on the sample size, consider changing this number (e.g. to 1/N and 1-1/N respectively)
          probs[probs==1]<-.99
          probs[probs==0]<-.01
          #Given this probability, return the associated value from a standard normal that falls into the same quantile
          testreads_normalized[d]<-qnorm(probs)
        }
        
        ################################
        ## Elastic-net model building ##
        ################################
        
        #Using N-fold internal CV, train the elastic net model using the training data
        #Note with larger sample sizes, N-fold internal CV becomes intractable
        model<-cv.glmnet(trainreads_norm,train_model,nfolds=NF,alpha=alpha[i],standardize=F)
        
        #Predict Ab using the test sample from parameters that minimized MSE during internal CV
        predicted_SAMP<-predict(model,newx=t(testreads_normalized),s="lambda.min")
        
        #Extract weights for this model
        weights_SAMP<-unlist(coef(model,lambda="lambda.min"))[,1]
        
        #Write out results for later concatenation
        write.table(predicted_SAMP,paste0(out_dir,"predicted/",models[m], "_predicted_QQ_n",NF,"_a",alpha[i],"_s",SAMP,".txt"),quote=F,row.names=F,col.names=F)
        write.table(weights_SAMP,paste0(out_dir,"weights/", models[m], "_weights_QQ_n",NF,"_a",alpha[i],"_s",SAMP,".txt"),quote=F,row.names=F,col.names=F)
        
        pred_samp[nrow(pred_samp)+1,] <- c(SAMP,(as.character(colnames(cts)[SAMP])),predicted_SAMP, weights_SAMP)
      }
      
      write.table(pred_samp, paste0(out_dir, models[m], "_predicted_ATAC_",cluster_names[j],"_n",NF,"_a",alpha[i],".txt"), quote = FALSE, sep = ",", row.names = F)
    }
  }
}
