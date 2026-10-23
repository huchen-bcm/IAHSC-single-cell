library(ggplot2)
library(stringr)

setwd("/scATAC_analyses/2_main_analyses/emmreml_downstream/Homer_annotate")
OUT_DIR <- "plot_homer_annotated_enriched_pathways_subgroups_combined/"
dir.create(OUT_DIR)


# Reactome pathways  -------------------------------------------------------
clusters <- c("CMP", "GMP", "MEP", "CLP", "PreBNK")
pathways_final <- vector()
fdr_final <- vector()
p_final <- vector()
clust_final <- vector()
for(i in 1:length(clusters)){
  
  results <- read.delim(file=paste0("homer_annotate_subgroups_combined_output/",clusters[i],"_homer_annotate_GO/reactome.txt"), header=TRUE)
  results$pval <- 10^(results$logP)
  results$q <- p.adjust(results$pval, method = "BH")
  results_sub1 <- subset(results, results$q < 0.05)
  results_sub2 <- subset(results_sub1, results_sub1$Target.Genes.in.Term >= 15)
  sig_pathways <- results_sub2$Term
  
  select <- rep("no", length(sig_pathways))
  select[which(grepl("GM-CSF",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Interleukin",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("TLR",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("MyD88",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Cytokine",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("myeloid",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("immune",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Immune",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("STAT",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Toll-Like Receptors",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("neutrophil",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("granulocyte",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("inflammasome",sig_pathways)==TRUE)] <- "yes"
  
  if(i == 1)
  {
    pathways_final <- sig_pathways[which(select=="yes")]
    fdr_final <- results_sub2$q[which(select=="yes")]
    p_final <- results_sub2$pval[which(select=="yes")]
    clust_final <- rep(paste0(clusters[i]), length(sig_pathways[which(select=="yes")]))
  }
  else
  {
    pathways_final <- append(pathways_final, sig_pathways[which(select=="yes")])
    fdr_final <- append(fdr_final, results_sub2$q[which(select=="yes")])
    p_final <- append(p_final, results_sub2$pval[which(select=="yes")])
    clust_final <- append(clust_final, rep(paste0(clusters[i]), length(sig_pathways[which(select=="yes")])))
  }
}

df_table <- data.frame(
  X <- pathways_final,
  Y <- clust_final,
  p <- -log10(p_final)
)
colnames(df_table) <- c("X", "Y", "p")

#tiff(paste0(OUT_DIR,"reactome_all_clusters_bubble.tiff"), units="in", width=6, height=6, res=300)
p <- ggplot(df_table, aes(x=X, y=Y, size=p, fill=p)) + geom_point(shape=21, color="black") + theme_classic() + labs(x="", y="") + scale_fill_gradient(high="red", low="white") + 
  theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90)) + coord_flip()
#dev.off()

final_file <- file.path(OUT_DIR,paste0("reactome_all_clusters_bubble.pdf"))
ggsave(filename = final_file, plot = p, height = 6, width = 6)


## Reactome short version

clusters <- c("CMP", "GMP", "MEP", "CLP", "PreBNK")
pathways_final <- vector()
fdr_final <- vector()
p_final <- vector()
clust_final <- vector()
for(i in 1:length(clusters)){
  
  results <- read.delim(file=paste0("homer_annotate_subgroups_combined_output/",clusters[i],"_homer_annotate_GO/reactome.txt"), header=TRUE)
  results$pval <- 10^(results$logP)
  results$q <- p.adjust(results$pval, method = "BH")
  results_sub1 <- subset(results, results$q < 0.05)
  results_sub2 <- subset(results_sub1, results_sub1$Target.Genes.in.Term >= 15)
  sig_pathways <- results_sub2$Term
  
  select <- rep("no", length(sig_pathways))
  select[which(grepl("GM-CSF",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Interleukin",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Cytokine",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("STAT",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Toll-Like Receptors",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("neutrophil",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("granulocyte",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("inflammasome",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("TLR4",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("TLR2",sig_pathways)==TRUE)] <- "yes"
  
  if(i == 1)
  {
    pathways_final <- sig_pathways[which(select=="yes")]
    fdr_final <- results_sub2$q[which(select=="yes")]
    p_final <- results_sub2$pval[which(select=="yes")]
    clust_final <- rep(paste0(clusters[i]), length(sig_pathways[which(select=="yes")]))
  }
  else
  {
    pathways_final <- append(pathways_final, sig_pathways[which(select=="yes")])
    fdr_final <- append(fdr_final, results_sub2$q[which(select=="yes")])
    p_final <- append(p_final, results_sub2$pval[which(select=="yes")])
    clust_final <- append(clust_final, rep(paste0(clusters[i]), length(sig_pathways[which(select=="yes")])))
  }
}

df_table <- data.frame(
  X <- pathways_final,
  Y <- clust_final,
  p <- -log10(p_final)
)
colnames(df_table) <- c("X", "Y", "p")

#tiff(paste0(OUT_DIR,"Fig_6A_reactome_short_bubble.tiff"), units="in", width=4.3, height=4, res=300)
p <- ggplot(df_table, aes(x=X, y=Y, size=p, fill=p)) + geom_point(shape=21, color="black") + theme_linedraw() + labs(x="", y="") + scale_fill_gradient(high="#990F20", low="white") + 
   theme(axis.text.x = element_text(angle = 45, hjust=1, vjust=1)) + coord_flip()
#dev.off()

final_file <- file.path(OUT_DIR,paste0("Fig_6A_reactome_short_bubble.pdf"))
ggsave(filename = final_file, plot = p, height = 4, width = 4.5)






# Biological process pathways  -------------------------------------------------------
clusters <- c("CMP", "GMP", "MEP", "CLP", "PreBNK")
pathways_final <- vector()
fdr_final <- vector()
p_final <- vector()
clust_final <- vector()
for(i in 1:length(clusters)){
  
  results <- read.delim(file=paste0("homer_annotate_subgroups_combined_output/",clusters[i],"_homer_annotate_GO/biological_process.txt"), header=TRUE)
  results$pval <- 10^(results$logP)
  results$q <- p.adjust(results$pval, method = "BH")
  results_sub1 <- subset(results, results$q < 0.05)
  results_sub2 <- subset(results_sub1, results_sub1$Target.Genes.in.Term >= 15)
  sig_pathways <- results_sub2$Term
  
  select <- rep("no", length(sig_pathways))
  select[which(grepl("GM-CSF",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Interleukin",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("TLR",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("MyD88",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Cytokine",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("myeloid",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("immune",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Immune",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("STAT",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Toll-Like Receptors",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("neutrophil",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("MAPK",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("granulocyte",sig_pathways)==TRUE)] <- "yes"
  
  
  if(i == 1)
  {
    pathways_final <- sig_pathways[which(select=="yes")]
    fdr_final <- results_sub2$q[which(select=="yes")]
    p_final <- results_sub2$pval[which(select=="yes")]
    clust_final <- rep(paste0(clusters[i]), length(sig_pathways[which(select=="yes")]))
  }
  else
  {
    pathways_final <- append(pathways_final, sig_pathways[which(select=="yes")])
    fdr_final <- append(fdr_final, results_sub2$q[which(select=="yes")])
    p_final <- append(p_final, results_sub2$pval[which(select=="yes")])
    clust_final <- append(clust_final, rep(paste0(clusters[i]), length(sig_pathways[which(select=="yes")])))
  }
}

lab <- rep("none", length(pathways_final))
lab[which(grepl("neutrophil",pathways_final)==TRUE)] <- "yes"
lab[which(grepl("granulocyte",pathways_final)==TRUE)] <- "yes"
lab[which(grepl("myeloid",pathways_final)==TRUE)] <- "yes"
lab[which(grepl("MAPK",pathways_final)==TRUE)] <- "yes"

df_table <- data.frame(
  X <- pathways_final,
  Y <- clust_final,
  p <- -log10(p_final),
  l <- lab
)
colnames(df_table) <- c("X", "Y", "p", "l")

label_cols <- rep("black", length(unique(df_table$X)))
path_order <- rep(0, length(unique(df_table$X)))
order_count <- 1

df_table$X[which(df_table$X == "regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains")] <- "reg adaptive immunity based on somatic recomb."

immune_pathways <- c("immune", "immunity", "immune system")
for(i in 1:length(immune_pathways))
{
  label_cols[which(grepl(immune_pathways[i], unique(df_table$X))==TRUE)] <- "#357EBD"
  path_order[which(grepl(immune_pathways[i], unique(df_table$X))==TRUE)] <- order_count
  order_count <- order_count + 1
}
myeloid_pathways <- c("myeloid")
for(i in 1:length(myeloid_pathways))
{
  label_cols[which(grepl(myeloid_pathways[i], unique(df_table$X))==TRUE)] <- "#5CB85C"
  path_order[which(grepl(myeloid_pathways[i], unique(df_table$X))==TRUE)] <- order_count
  order_count <- order_count + 1
}
gran_pathways <- c("neutrophil", "granulocyte")
for(i in 1:length(gran_pathways))
{
  label_cols[which(grepl(gran_pathways[i], unique(df_table$X))==TRUE)] <- "#D43F3A"
  path_order[which(grepl(gran_pathways[i], unique(df_table$X))==TRUE)] <- order_count
  order_count <- order_count + 1
}

#tiff(paste0(OUT_DIR,"biological_process_all_clusters_bubble.tiff"), units="in", width=7, height=6.5, res=300)
p <- ggplot(df_table, aes(x=X, y=Y, size=p, fill=p, color=l)) + geom_point(shape=21, aes(color=l), stroke=1) + theme_classic() + labs(x="", y="") + scale_fill_gradient(high="blue", low="white") + 
  theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + scale_color_manual(values=c("black", "#EE4C97")) 
#dev.off()

final_file <- file.path(OUT_DIR,paste0("biological_process_all_clusters_bubble.pdf"))
ggsave(filename = final_file, plot = p, height = 6.5, width = 7)


##alternate

df_table$X <- factor(df_table$X, unique(df_table$X)[order(path_order)])  

#tiff(paste0(OUT_DIR,"biological_process_all_clusters_bubble_alternate.tiff"), units="in", width=7, height=6.5, res=300)
p <- ggplot(df_table, aes(x=X, y=Y, size=p, fill=p)) + geom_point(shape=21,color="black", stroke=1) + theme_classic() + labs(x="", y="") + scale_fill_gradient(high="blue", low="white") + 
  theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + scale_color_manual(values=c("black", "#EE4C97")) + theme(axis.text.y = element_text(colour=label_cols[order(path_order)])) 
#dev.off()

final_file <- file.path(OUT_DIR,paste0("biological_process_all_clusters_bubble_alternate.pdf"))
ggsave(filename = final_file, plot = p, height = 6, width = 7)


# Biological process pathways NEUTROPHIL only -------------------------------------------------------
clusters <- c("CMP", "GMP", "MEP", "CLP", "PreBNK")
pathways_final <- vector()
fdr_final <- vector()
p_final <- vector()
clust_final <- vector()

for(i in 1:length(clusters)){
  
  results <- read.delim(file=paste0("homer_annotate_subgroups_combined_output/",clusters[i],"_homer_annotate_GO/biological_process.txt"), header=TRUE)
  results$pval <- 10^(results$logP)
  results$q <- p.adjust(results$pval, method = "BH")
  results_sub1 <- subset(results, results$q < 0.05)
  results_sub2 <- subset(results_sub1, results_sub1$Target.Genes.in.Term >= 15)
  sig_pathways <- results_sub2$Term
  
  select <- rep("no", length(sig_pathways))
  select[which(grepl("Interleukin",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("TLR",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("MyD88",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Cytokine",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("STAT",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("Toll-Like Receptors",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("neutrophil",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("granulocyte",sig_pathways)==TRUE)] <- "yes"
  select[which(grepl("inflammasome",sig_pathways)==TRUE)] <- "yes"
  
  
  if(i == 1)
  {
    pathways_final <- sig_pathways[which(select=="yes")]
    fdr_final <- results_sub2$q[which(select=="yes")]
    p_final <- results_sub2$pval[which(select=="yes")]
    clust_final <- rep(paste0(clusters[i]), length(sig_pathways[which(select=="yes")]))
  }
  else
  {
    pathways_final <- append(pathways_final, sig_pathways[which(select=="yes")])
    fdr_final <- append(fdr_final, results_sub2$q[which(select=="yes")])
    p_final <- append(p_final, results_sub2$pval[which(select=="yes")])
    clust_final <- append(clust_final, rep(paste0(clusters[i]), length(sig_pathways[which(select=="yes")])))
  }
}

lab <- rep("none", length(pathways_final))
lab[which(grepl("neutrophil",pathways_final)==TRUE)] <- "yes"
lab[which(grepl("granulocyte",pathways_final)==TRUE)] <- "yes"


df_table <- data.frame(
  X <- pathways_final,
  Y <- clust_final,
  p <- -log10(p_final),
  l <- lab
)
colnames(df_table) <- c("X", "Y", "p", "l")

label_cols <- rep("black", length(unique(df_table$X)))
path_order <- rep(0, length(unique(df_table$X)))
order_count <- 1

df_table$X[which(df_table$X == "regulation of adaptive immune response based on somatic recombination of immune receptors built from immunoglobulin superfamily domains")] <- "reg adaptive immunity based on somatic recomb."


#tiff(paste0(OUT_DIR,"biological_process_neutrophil_bubble.tiff"), units="in", width=5, height=2.5, res=300)
p <- ggplot(df_table, aes(x=X, y=Y, size=p, fill=p, color=l)) + geom_point(shape=21, aes(color=l), stroke=1) + theme_classic() + labs(x="", y="") + scale_fill_gradient(high="blue", low="white") + 
  theme(legend.position="none") + theme(axis.text.x = element_text(angle = 90)) + coord_flip() + scale_color_manual(values=c("black", "#EE4C97")) 
#dev.off()

final_file <- file.path(OUT_DIR,paste0("biological_process_neutrophil_bubble.pdf"))
ggsave(filename = final_file, plot = p, height = 2.5, width = 5)

