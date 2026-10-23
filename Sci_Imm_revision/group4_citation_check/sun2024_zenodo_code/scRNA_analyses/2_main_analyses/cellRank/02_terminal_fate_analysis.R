#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded


library(devtools)
library(velocyto.R)
library(pagoda2)
library(Seurat)
library(hdf5r)
library(SeuratDisk)
library(SingleCellExperiment)
library(ggplot2)
library(ggrepel)
library(stringr)

# setup -------------------------------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/cellRank/") 
IN_DIR <- "wVelocity_input/"
IN_DIR2 <- "terminal_state_probs/"
OUT_DIR <- "terminal_state_analysis_output/"
dir.create(OUT_DIR)


# vector of donors --------------------------------------------------------

donors <- c("S1", "S2", "S3", "S5", "S6", "S7", "S8", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S20", "S21")


# create vectors for D0 and D90 lineage biases  ---------------------------

cmpa_d90 <- vector(length=length(donors))
cmpb_d90 <- vector(length=length(donors))
cmpc_d90 <- vector(length=length(donors))
gmpa_d90 <- vector(length=length(donors))
gmpb_d90 <- vector(length=length(donors))
mlpa_d90 <- vector(length=length(donors))
mlpb_d90 <- vector(length=length(donors))
mepa_d90 <- vector(length=length(donors))
mepb_d90 <- vector(length=length(donors))
mepc_d90 <- vector(length=length(donors))
cmp_d90 <- vector(length=length(donors))
preBNK_d90 <- vector(length=length(donors))
lymphoid_d90 <- vector(length=length(donors))
myeloid_d90 <- vector(length=length(donors))
erythroid_d90 <- vector(length=length(donors))


cmpa_d0 <- vector(length=length(donors))
cmpb_d0 <- vector(length=length(donors))
cmpc_d0 <- vector(length=length(donors))
gmpa_d0 <- vector(length=length(donors))
gmpb_d0 <- vector(length=length(donors))
mlpa_d0 <- vector(length=length(donors))
mlpb_d0 <- vector(length=length(donors))
mepa_d0 <- vector(length=length(donors))
mepb_d0 <- vector(length=length(donors))
mepc_d0 <- vector(length=length(donors))
preBNK_d0 <- vector(length=length(donors))
erythroid_d0 <- vector(length=length(donors))
lymphoid_d0 <- vector(length=length(donors))
myeloid_d0 <- vector(length=length(donors))


cmpa_change <- vector(length=length(donors))
cmpb_change <- vector(length=length(donors))
cmpc_change <- vector(length=length(donors))
gmpa_change <- vector(length=length(donors))
gmpb_change <- vector(length=length(donors))
mlpa_change <- vector(length=length(donors))
mlpb_change <- vector(length=length(donors))
mepa_change <- vector(length=length(donors))
mepb_change <- vector(length=length(donors))
mepc_change <- vector(length=length(donors))
preBNK_change <- vector(length=length(donors))
lymphoid_change <- vector(length=length(donors))
myeloid_change <- vector(length=length(donors))
cmp_change <- vector(length=length(donors))
erythroid_change <- vector(length=length(donors))

plot_list <- list()
k <- 1
plot_list2 <- list()


################################################################################
####### Big for loop to record terminal fate at D0 and D90 for each donor ######
################################################################################

for(i in 1:length(donors)){
  
  # two timepoints
  time <- c("Td0", "Tm3")
  time_label <- c("D0", "D90")
  
  # initialize vectors to record % HSCs biased to each lineage at each timepoint 
  cmpa_perc <- vector(length=length(time))
  cmpb_perc <- vector(length=length(time))
  cmpc_perc <- vector(length=length(time))
  gmpa_perc <- vector(length=length(time))
  gmpb_perc <- vector(length=length(time))
  mlpa_perc <- vector(length=length(time))
  mlpb_perc <- vector(length=length(time))
  mepa_perc <- vector(length=length(time))
  mepb_perc <- vector(length=length(time))
  mepc_perc <- vector(length=length(time))
  preBNK_perc <- vector(length=length(time))
  lymphoid_perc <- vector(length=length(time))
  myeloid_perc <- vector(length=length(time))
  cmp_perc <- vector(length=length(time))
  erythroid_perc <- vector(length=length(time))
  
  # for each donor loop through each timepoint to determine lineage biases at D0 and D90
  for(j in 1:length(time)){
    # coordinate contains the UMAP coordinates for each cell
    coordinates <- read.csv(file=paste0(IN_DIR, "cell_embeddings_", donors[i], "_", time[j],".csv"), header=TRUE)
    
    # states the cluster membership of each cell 
    id <- read.csv(file=paste0(IN_DIR, "clusters_", donors[i], "_", time[j], ".csv"), header=TRUE)
    
    # predicted terminal state for each cell
    terminal_states <- read.csv(file=paste0(IN_DIR2, donors[i],"_",time[j],"_to_terminal_states.csv"), header=TRUE)
    
    # make sure all cells match between ID and terminal states files 
    id <- id[which(id$X %in% terminal_states$Cell.ID),]
    id <- id[order(match(as.character(id$X), as.character(terminal_states$Cell.ID))),]
    which(id$X != terminal_states$Cell.ID)
    colnames(id) <- c("X", "clust_names")
    
    # make sure all cells match bewteen coordinates and terminal states 
    row.names(coordinates) <- coordinates$X
    coordinates <- coordinates[which(coordinates$X %in% terminal_states$Cell.ID),]
    coordinates <- coordinates[order(match(as.character(coordinates$X), as.character(terminal_states$Cell.ID))),]
    which(coordinates$X != terminal_states$Cell.ID)
    
    
    # for each cell recored the most likely terminal state 
    terminal_states <- subset(terminal_states, select = -X)
    terminal_states <- subset(terminal_states, select = -Cell.ID)
    most_likely <- colnames(terminal_states)[max.col(terminal_states, ties.method = "first")]
    
    # organize all this data in a data frame 
    df <- data.frame(
      name <- id$X,
      clust <- id$clust_names,
      x <- coordinates$UMAP_1,
      y <- coordinates$UMAP_2,
      col <- most_likely
    )
    colnames(df) <- c("name", "clust", "x", "y", "col")
    
    # include only fate predictions for HSCs 
    df_final <- subset(df, df$clust %in% c("HSC_a", "HSC_b"))
    df_final$col <- factor(df_final$col, colnames(terminal_states))
    

    
    myColorScale <- c("#006DDB", "#FF7F0E", "#00BB00", "#FF0000", "#BB00BB", "#924900", "#FF6DB6", "#749B58", "#00A9CC", "#AEC7E8", "#FFCC65")
    names(myColorScale) <- colnames(terminal_states)
    
    # plot UMAP
    p <- ggplot(df_final, aes(x=x, y=y, color=col)) + geom_point(aes(fill=col),shape=21, size=2, alpha=1, stroke=0.4, color="black") + theme_void() +
      labs(title=paste0("Donor ", donors[i]," ",time_label[j]," \nPredicted terminal state"), x= "UMAP1", y="UMAP2") + scale_fill_manual(values=myColorScale) + theme(legend.position = "none")   
    
    plot_list[[k]] <- p
    k <- k + 1
    
    # record percentages of HSCs that are biased towards each terminal fate 
    cmpa_perc[j] <- length(df_final$col[which(df_final$col == "CMP_a")])/length(df_final$col)
    cmpb_perc[j] <- length(df_final$col[which(df_final$col == "CMP_b")])/length(df_final$col)
    cmpc_perc[j] <- length(df_final$col[which(df_final$col == "CMP_c")])/length(df_final$col)
    gmpa_perc[j] <- length(df_final$col[which(df_final$col == "GMP_a")])/length(df_final$col)
    gmpb_perc[j] <- length(df_final$col[which(df_final$col == "GMP_b")])/length(df_final$col)
    mlpa_perc[j] <- length(df_final$col[which(df_final$col == "MLP_a")])/length(df_final$col)
    mlpb_perc[j] <- length(df_final$col[which(df_final$col == "MLP_b")])/length(df_final$col)
    preBNK_perc[j] <- length(df_final$col[which(df_final$col == "PreBNK")])/length(df_final$col)
    mepa_perc[j] <- length(df_final$col[which(df_final$col == "MEP_a")])/length(df_final$col)
    mepb_perc[j] <- length(df_final$col[which(df_final$col == "MEP_b")])/length(df_final$col)
    mepc_perc[j] <- length(df_final$col[which(df_final$col == "MEP_c")])/length(df_final$col)
    myeloid_perc[j] <- cmpa_perc[j] + cmpb_perc[j] + cmpc_perc[j] + gmpa_perc[j] + gmpb_perc[j]
    cmp_perc[j] <- cmpa_perc[j] + cmpb_perc[j] + cmpc_perc[j]
    lymphoid_perc[j] <- mlpa_perc[j] + mlpb_perc[j] + preBNK_perc[j]
    erythroid_perc[j] <- mepa_perc[j] + mepb_perc[j] + mepc_perc[j]
  }
  
  cmpa_change[i] <- cmpa_perc[2] - cmpa_perc[1]
  cmpb_change[i] <- cmpb_perc[2] - cmpb_perc[1]
  cmpc_change[i] <- cmpc_perc[2] - cmpc_perc[1]
  gmpa_change[i] <- gmpa_perc[2] - gmpa_perc[1]
  gmpb_change[i] <- gmpb_perc[2] - gmpb_perc[1]
  mlpa_change[i] <- mlpa_perc[2] - mlpa_perc[1]
  mlpb_change[i] <- mlpb_perc[2] - mlpb_perc[1]
  preBNK_change[i] <- preBNK_perc[2] - preBNK_perc[1]
  mepa_change[i] <- mepa_perc[2] - mepa_perc[1]
  mepb_change[i] <- mepb_perc[2] - mepb_perc[1]
  mepc_change[i] <- mepc_perc[2] - mepc_perc[1]
  myeloid_change[i] <- myeloid_perc[2] - myeloid_perc[1]
  cmp_change[i] <- cmp_perc[2] - cmp_perc[1]
  lymphoid_change[i] <- lymphoid_perc[2] - lymphoid_perc[1]
  erythroid_change[i] <- erythroid_perc[2] - erythroid_perc[1]

  
  #record day90 percentages
  cmpa_d90[i] <- cmpa_perc[2] 
  cmpb_d90[i] <- cmpb_perc[2] 
  cmpc_d90[i] <- cmpc_perc[2] 
  gmpa_d90[i] <- gmpa_perc[2] 
  gmpb_d90[i] <- gmpb_perc[2] 
  mlpa_d90[i] <- mlpa_perc[2] 
  mlpb_d90[i] <- mlpb_perc[2] 
  mepa_d90[i] <- mepa_perc[2] 
  mepb_d90[i] <- mepb_perc[2] 
  mepc_d90[i] <- mepc_perc[2] 
  preBNK_d90[i] <- preBNK_perc[2] 
  myeloid_d90[i] <- myeloid_perc[2] 
  lymphoid_d90[i] <- lymphoid_perc[2] 
  erythroid_d90[i] <- erythroid_perc[2] 
  cmp_d90[i] <- cmp_perc[2] 
  
  #record day0 percentages
  cmpa_d0[i] <- cmpa_perc[1] 
  cmpb_d0[i] <- cmpb_perc[1] 
  cmpc_d0[i] <- cmpc_perc[1] 
  gmpa_d0[i] <- gmpa_perc[1] 
  gmpb_d0[i] <- gmpb_perc[1] 
  mlpa_d0[i] <- mlpa_perc[1] 
  mlpb_d0[i] <- mlpb_perc[1] 
  preBNK_d0[i] <- preBNK_perc[1] 
  mepa_d0[i] <- mepa_perc[1] 
  mepb_d0[i] <- mepb_perc[1] 
  mepc_d0[i] <- mepc_perc[1] 
  myeloid_d90[i] <- myeloid_perc[1] 
  lymphoid_d90[i] <- lymphoid_perc[1] 
  erythroid_d90[i] <- erythroid_perc[1] 
  
  # data frame for each donor that shows lineage biases (values column: percent of HSCs with bias) at D0 and D90
  df2 <- data.frame(
    timepoint <- c(rep("D0", length(colnames(terminal_states))), rep("D90", length(colnames(terminal_states))) ),
    Celltype <- rep(paste0(colnames(terminal_states)), 2),
    values <- c(cmpa_perc[1], cmpb_perc[1], cmpc_perc[1], gmpa_perc[1], gmpb_perc[1], mepa_perc[1], mepb_perc[1], mepc_perc[1], mlpa_perc[1], mlpb_perc[1], preBNK_perc[1], cmpa_perc[2], cmpb_perc[2], cmpc_perc[2], gmpa_perc[2], gmpb_perc[2], mepa_perc[2], mepb_perc[2], mepc_perc[2], mlpa_perc[2], mlpb_perc[2], preBNK_perc[2])
  )
  colnames(df2) <- c("timepoint", "Celltype", "values")
  
}



################################################################################
####################            Plot the results            #################### 
################################################################################


donors <- c("S1", "S2", "S3", "S5", "S6", "S7", "S8", "S10", "S11", "S12", "S13", "S14", "S15", "S16", "S17", "S18", "S20", "S21")
wilcox_pvals <- vector(length=11)
names(wilcox_pvals) <- c('CMP_a', 'CMP_b', 'CMP_c', 'GMP_a', "GMP_b", 'MEP_a', "MEP_b", "MEP_c", "MLP_a", "MLP_b", "PreBNK")


# filter donors with fewer than 20 cells for for the DE gene analysis 
filter_20 <- readRDS(file=paste0("label_transfer/filter_20_list.rds"))
donors_to_filter <- unique(c(filter_20[[1]], filter_20[[2]]))




# CMP_b boxplots ----------------------------------------------------------

df_cmpb <- data.frame(
  val <- cmpb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

# subset to exclude donors 2, 3, and 12 (these had fewer than 20 cells)
df_cmpb <- subset(df_cmpb, df_cmpb$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpb) <- c("val", "d", "group")

tiff(paste0(OUT_DIR, "cmp_b_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpb, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group ), color = "black") + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_b Terminal state prediction") + scale_fill_manual(values=c("#FFD320", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()


# compare BCG and CTL D90
df_cmpb_d90 <- data.frame(
  vals_d90 <- cmpb_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

# subset to exclude donors 2, 3, and 12 (these had fewer than 20 cells)
df_cmpb_d90 <- subset(df_cmpb_d90, df_cmpb_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))

colnames(df_cmpb_d90) <- c("vals_d90", "d", "group")

wilcox_pvals['CMP_b'] <- wilcox_test <- wilcox.test(df_cmpb_d90$vals_d90[which(df_cmpb_d90$group == "BCG")], df_cmpb_d90$vals_d90[which(df_cmpb_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "cmp_b_term_state_d90.tiff"), units="in", width=2.4, height=3, res=250)
ggplot(data=df_cmpb_d90, aes(x=group, y=vals_d90)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_classic()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent CMPb biased HSCs", fill="Group", title=paste0("CMP_b bias Day 90; \npval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#FFD320", "#004586"))+
  theme(legend.position = "none")
dev.off()



# CMP_a boxplots ----------------------------------------------------------


df_cmpa <- data.frame(
  val <- cmpa_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)


df_cmpa <- subset(df_cmpa, df_cmpa$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpa) <- c("val", "d", "group")

tiff(paste0(OUT_DIR, "cmp_a_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group , color="black")) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_a Terminal state prediction") + scale_fill_manual(values=c("#FF7F0E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()



#d90 BCG vs. placebo comparison plots

df_cmpa_d90 <- data.frame(
  vals_d90 <- cmpa_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

df_cmpa_d90 <- subset(df_cmpa_d90, df_cmpa_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpa_d90) <- c("vals_d90", "d", "group")

wilcox_pvals['CMP_a'] <- wilcox.test(df_cmpa_d90$vals_d90[which(df_cmpa_d90$group == "BCG")], df_cmpa_d90$vals_d90[which(df_cmpa_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "cmp_a_term_state_d90.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpa_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent cmpa biased HSCs", fill="Group", title=paste0("CMP_a bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#FF7F0E", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()




# CMP_c boxplots ----------------------------------------------------------

df_cmpc <- data.frame(
  val <- cmpc_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)


df_cmpc <- subset(df_cmpc, df_cmpc$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpc) <- c("val", "d", "group")

tiff(paste0(OUT_DIR, "cmp_c_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpc, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="CMP_c Terminal state prediction") + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()



#d90 BCG vs. placebo comparison plots

df_cmpc_d90 <- data.frame(
  vals_d90 <- cmpc_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

df_cmpc_d90 <- subset(df_cmpc_d90, df_cmpc_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_cmpc_d90) <- c("vals_d90", "d", "group")

wilcox_pvals['CMP_c'] <- wilcox.test(df_cmpc_d90$vals_d90[which(df_cmpc_d90$group == "BCG")], df_cmpc_d90$vals_d90[which(df_cmpc_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "cmp_c_term_state_d90.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_cmpc_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent cmpc biased HSCs", fill="Group", title=paste0("CMP_c bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#CC9900", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()


# mepa boxplots ----------------------------------------------------------

df_mepa <- data.frame(
  val <- mepa_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

df_mepa <- subset(df_mepa, df_mepa$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_mepa) <- c("val", "d", "group")

tiff(paste0(OUT_DIR, "mepa_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mepa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MEP_a Terminal state prediction") + scale_fill_manual(values=c("#CE6DBD", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()



#d90 BCG vs. placebo comparison plots

df_mepa_d90 <- data.frame(
  vals_d90 <- mepa_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)
wilcox_pvals['MEP_a'] <- wilcox.test(df_mepa_d90$vals_d90[which(df_mepa_d90$group == "BCG")], df_mepa_d90$vals_d90[which(df_mepa_d90$group == "CTL")])$p.value 

ggplot(data=df_mepa_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent mepa biased HSCs", fill="Group", title=paste0("mepa bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#CE6DBD", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)


# mepb boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_mepb <- data.frame(
  val <- mepb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
#d90 BCG vs. placebo comparison plots

df_mepb_d90 <- data.frame(
  vals_d90 <- mepb_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
wilcox_pvals['MEP_b'] <- wilcox.test(df_mepb_d90$vals_d90[which(df_mepb_d90$group == "BCG")], df_mepb_d90$vals_d90[which(df_mepb_d90$group == "CTL")])$p.value 


# mepc boxplots ----------------------------------------------------------

lab <- rep("no",length(donors))
lab[which(!(donors %in% c("S8", "S10", "S12", "S14", "S16")))] <- "yes"
df_mepc <- data.frame(
  val <- mepc_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
#d90 BCG vs. placebo comparison plots

df_mepc_d90 <- data.frame(
  vals_d90 <- mepc_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG"),
  l <- lab
)
wilcox_pvals['MEP_c'] <- wilcox.test(df_mepc_d90$vals_d90[which(df_mepc_d90$group == "BCG")], df_mepc_d90$vals_d90[which(df_mepc_d90$group == "CTL")])$p.value 




# GMP_a boxplots ----------------------------------------------------------

df_gmpa <- data.frame(
  val <- gmpa_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

df_gmpa <- subset(df_gmpa, df_gmpa$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_gmpa) <- c("val", "d", "group")

tiff(paste0(OUT_DIR, "gmp_a_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_gmpa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="GMP_a Terminal state prediction") + scale_fill_manual(values=c("#33CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()


#d90 BCG vs. placebo comparison plots

df_gmpa_d90 <- data.frame(
  vals_d90 <- gmpa_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

df_gmpa_d90 <- subset(df_gmpa_d90, df_gmpa_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_gmpa_d90) <- c("vals_d90", "d", "group")

wilcox_pvals['GMP_a'] <- wilcox.test(df_gmpa_d90$vals_d90[which(df_gmpa_d90$group == "BCG")], df_gmpa_d90$vals_d90[which(df_gmpa_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "gmpa_term_state_d90.tiff"), units="in", width=5, height=4, res=250)
ggplot(data=df_gmpa_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent gmpa biased HSCs", fill="Group", title=paste0("GMP_a bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#33CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()



# gmpb boxplots ----------------------------------------------------------

df_gmpb <- data.frame(
  val <- gmpb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)


df_gmpb <- subset(df_gmpb, df_gmpb$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_gmpb) <- c("val", "d", "group")

tiff(paste0(OUT_DIR, "gmpb_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_gmpb, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group )) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="GMP_b Terminal state prediction") + scale_fill_manual(values=c("#99CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()

#d90 BCG vs. placebo comparison plots

df_gmpb_d90 <- data.frame(
  vals_d90 <- gmpb_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

df_gmpb_d90 <- subset(df_gmpb_d90, df_gmpb_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_gmpb_d90) <- c("vals_d90", "d", "group")

wilcox_pvals['GMP_b'] <- wilcox.test(df_gmpb_d90$vals_d90[which(df_gmpb_d90$group == "BCG")], df_gmpb_d90$vals_d90[which(df_gmpb_d90$group == "CTL")])$p.value 

tiff(paste0(OUT_DIR, "gmpb_term_state_d90.tiff"), units="in", width=5, height=4, res=250)
ggplot(data=df_gmpb_d90, aes(x=group, y=vals_d90, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Percent gmpb biased HSCs", fill="Group", title=paste0("GMP_b bias Day 90; pval = ", round(wilcox_test,3))) + scale_fill_manual(values=c("#99CC00", "#004586"))+
  scale_color_manual(values=c( "#83CAFF", "black")) + geom_text_repel(size=2)
dev.off()


# mlpa boxplots ----------------------------------------------------------

df_mlpa <- data.frame(
  val <- mlpa_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

df_mlpa <- subset(df_mlpa, df_mlpa$d....donors %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_mlpa) <- c("val", "d", "group")

tiff(paste0(OUT_DIR, "mlpa_term_state_labels_2_3_12_removed.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_mlpa, aes(x=group, y=val, label=d)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="MLP_a Terminal state prediction") + scale_fill_manual(values=c("#83CAFF", "#004586"))+
  scale_color_manual(values=c( "#DD4477", "black")) + geom_text_repel(size=2)
dev.off()

#d90 BCG vs. placebo comparison plots

df_mlpa_d90 <- data.frame(
  vals_d90 <- mlpa_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

df_mlpa_d90 <- subset(df_mlpa_d90, df_mlpa_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_mlpa_d90) <- c("vals_d90", "d", "group")

wilcox_pvals['MLP_a'] <- wilcox.test(df_mlpa_d90$vals_d90[which(df_mlpa_d90$group == "BCG")], df_mlpa_d90$vals_d90[which(df_mlpa_d90$group == "CTL")])$p.value 



# mlpb boxplots ----------------------------------------------------------

df_mlpb <- data.frame(
  val <- mlpb_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)
#d90 BCG vs. placebo comparison plots

df_mlpb_d90 <- data.frame(
  vals_d90 <- mlpb_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

df_mlpb_d90 <- subset(df_mlpb_d90, df_mlpb_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_mlpb_d90) <- c("vals_d90", "d", "group")

wilcox_pvals['MLP_b'] <- wilcox.test(df_mlpb_d90$vals_d90[which(df_mlpb_d90$group == "BCG")], df_mlpb_d90$vals_d90[which(df_mlpb_d90$group == "CTL")])$p.value 



# pre-BNK boxplots ----------------------------------------------------------

df_preBNK <- data.frame(
  val <- preBNK_change,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

tiff(paste0(OUT_DIR, "preBNK_term_state.tiff"), units="in", width=3.5, height=3, res=250)
ggplot(data=df_preBNK, aes(x=group, y=val)) +
  geom_boxplot(aes(fill=group), alpha=0.5) + geom_jitter(shape=21, size=3, alpha=1,stroke=1, position=position_jitter(0.2), aes(fill=group)) + 
  theme_light()+geom_hline(yintercept=0, color="brown") + labs(x="", y="Change in percentage (Tm3-Td0)", fill="Group", title="PreBNK Terminal state prediction") + scale_fill_manual(values=c("#83CAFF", "#004586"))+
  scale_color_manual(values=c( "#DD4477", "black"))
dev.off()


#d90 BCG vs. placebo comparison plots

df_preBNK_d90 <- data.frame(
  vals_d90 <- preBNK_d90,
  d <- donors,
  group <- c("BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "BCG", "CTL", "BCG", "CTL", "BCG", "BCG", "BCG", "BCG", "CTL", "BCG", "BCG")
)

df_preBNK_d90 <- subset(df_preBNK_d90, df_preBNK_d90$d %in% c("S1", "S5" , "S6" , "S7" , "S8" , "S10" ,"S11" , "S13", "S14" ,"S15" ,"S16", "S17" ,"S18" ,"S20", "S21"))
colnames(df_preBNK_d90) <- c("vals_d90", "d", "group")

wilcox_pvals['PreBNK'] <- wilcox.test(df_preBNK_d90$vals_d90[which(df_preBNK_d90$group == "BCG")], df_preBNK_d90$vals_d90[which(df_preBNK_d90$group == "CTL")])$p.value 


