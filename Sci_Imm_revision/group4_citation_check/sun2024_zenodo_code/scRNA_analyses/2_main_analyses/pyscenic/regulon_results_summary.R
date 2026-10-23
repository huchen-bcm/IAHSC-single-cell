setwd("/scRNA_analyses/2_main_analyses/pyscenic")
OUT_DIR <- "get_diff_regulons_outs/"

## p < 0.15
HSC_clusters <- c("HSC_a", "HSC_b")
sig_regulons_hsc <- vector()
for(i in 1:length(HSC_clusters))
{
  dat <- read.csv(file=paste0("get_diff_regulons_outs/",HSC_clusters[i],"_wilcox_data_no_outliers.csv"))
  sig_regulons_hsc <- append(dat$X[which(dat$pval_wilcox < 0.15)], sig_regulons_hsc)
}
sig_regulons_hsc <- unique(sig_regulons_hsc)



CMP_clusters <- c("CMP_a", "CMP_b", "CMP_c")
sig_regulons_cmp <- vector()
for(i in 1:length(CMP_clusters))
{
  dat <- read.csv(file=paste0("get_diff_regulons_outs/",CMP_clusters[i],"_wilcox_data_no_outliers.csv"))
  sig_regulons_cmp <- append(dat$X[which(dat$pval_wilcox < 0.15)], sig_regulons_cmp)
}
sig_regulons_cmp <- unique(sig_regulons_cmp)



GMP_clusters <- c("GMP_a", "GMP_b")
sig_regulons_gmp <- vector()
for(i in 1:length(GMP_clusters))
{
  dat <- read.csv(file=paste0("get_diff_regulons_outs/",GMP_clusters[i],"_wilcox_data_no_outliers.csv"))
  sig_regulons_gmp <- append(dat$X[which(dat$pval_wilcox <0.15)], sig_regulons_gmp)
}
sig_regulons_gmp <- unique(sig_regulons_gmp)


MEP_clusters <- c("MEP_a", "MEP_b", "MEP_c")
sig_regulons_mep <- vector()
for(i in 1:length(MEP_clusters))
{
  dat <- read.csv(file=paste0("get_diff_regulons_outs/",MEP_clusters[i],"_wilcox_data_no_outliers.csv"))
  sig_regulons_mep <- append(dat$X[which(dat$pval_wilcox < 0.15)], sig_regulons_mep)
}
sig_regulons_mep <- unique(sig_regulons_mep)



MLP_clusters <- c("MLP_a", "MLP_b")
sig_regulons_mlp <- vector()
for(i in 1:length(MLP_clusters))
{
  dat <- read.csv(file=paste0("get_diff_regulons_outs/",MLP_clusters[i],"_wilcox_data_no_outliers.csv"))
  sig_regulons_mlp <- append(dat$X[which(dat$pval_wilcox < 0.15)], sig_regulons_mlp)
}
sig_regulons_mlp <- unique(sig_regulons_mlp)


dat <- read.csv(file=paste0("get_diff_regulons_outs/PreBNK_wilcox_data_no_outliers.csv"))
sig_regulons_PreBNK <- dat$X[which(dat$pval_wilcox < 0.15)]


df_bar <- data.frame(
  count <- c(length(sig_regulons_hsc), length(sig_regulons_cmp), length(sig_regulons_gmp), length(sig_regulons_mep), length(sig_regulons_mlp), length(sig_regulons_PreBNK)),
  name <- c("HSC", "CMP", "GMP", "MEP", "MLP", "PreBNK")
)
df_bar$name <- factor(df_bar$name, c("HSC", "CMP", "GMP", "MEP", "MLP", "PreBNK"))

tiff(paste0(OUT_DIR, "numRegulons_p0.15_bargraph.tiff"), units="in", width=3, height=3, res=250)
ggplot(df_bar, aes(x=name, y=count)) + geom_bar(stat="identity", aes(fill=name), color="black") + theme_classic() + theme(legend.position = "none") + labs(x="", y="Regulons (p < 0.15)")+
  scale_fill_manual(values=c('#EF5350','#FF7F0E', '#33CC00', '#8C5782', '#46B8DA', '#1A0099')) + theme(axis.text.x = element_text(angle = 90))
dev.off()



### p <= 0.05

HSC_clusters <- c("HSC_a", "HSC_b")
sig_regulons_hsc <- vector()
for(i in 1:length(HSC_clusters))
{
  dat <- read.csv(file=paste0("get_diff_regulons_outs/",HSC_clusters[i],"_wilcox_data_no_outliers.csv"))
  sig_regulons_hsc <- append(dat$X[which(dat$pval_wilcox < 0.05)], sig_regulons_hsc)
}
sig_regulons_hsc <- unique(sig_regulons_hsc)



CMP_clusters <- c("CMP_a", "CMP_b", "CMP_c")
sig_regulons_cmp <- vector()
for(i in 1:length(CMP_clusters))
{
  dat <- read.csv(file=paste0("get_diff_regulons_outs/",CMP_clusters[i],"_wilcox_data_no_outliers.csv"))
  sig_regulons_cmp <- append(dat$X[which(dat$pval_wilcox < 0.05)], sig_regulons_cmp)
}
sig_regulons_cmp <- unique(sig_regulons_cmp)



GMP_clusters <- c("GMP_a", "GMP_b")
sig_regulons_gmp <- vector()
for(i in 1:length(GMP_clusters))
{
  dat <- read.csv(file=paste0("get_diff_regulons_outs/",GMP_clusters[i],"_wilcox_data_no_outliers.csv"))
  sig_regulons_gmp <- append(dat$X[which(dat$pval_wilcox < 0.05)], sig_regulons_gmp)
}
sig_regulons_gmp <- unique(sig_regulons_gmp)


MEP_clusters <- c("MEP_a", "MEP_b", "MEP_c")
sig_regulons_mep <- vector()
for(i in 1:length(MEP_clusters))
{
  dat <- read.csv(file=paste0("get_diff_regulons_outs/",MEP_clusters[i],"_wilcox_data_no_outliers.csv"))
  sig_regulons_mep <- append(dat$X[which(dat$pval_wilcox < 0.05)], sig_regulons_mep)
}
sig_regulons_mep <- unique(sig_regulons_mep)



MLP_clusters <- c("MLP_a", "MLP_b")
sig_regulons_mlp <- vector()
for(i in 1:length(MLP_clusters))
{
  dat <- read.csv(file=paste0("get_diff_regulons_outs/",MLP_clusters[i],"_wilcox_data_no_outliers.csv"))
  sig_regulons_mlp <- append(dat$X[which(dat$pval_wilcox < 0.05)], sig_regulons_mlp)
}
sig_regulons_mlp <- unique(sig_regulons_mlp)


dat <- read.csv(file=paste0("get_diff_regulons_outs/PreBNK_wilcox_data_no_outliers.csv"))
sig_regulons_PreBNK <- dat$X[which(dat$pval_wilcox < 0.05)]


df_bar <- data.frame(
  count <- c(length(sig_regulons_hsc), length(sig_regulons_cmp), length(sig_regulons_gmp), length(sig_regulons_mep), length(sig_regulons_mlp), length(sig_regulons_PreBNK)),
  name <- c("HSC", "CMP", "GMP", "MEP", "MLP", "PreBNK")
)
df_bar$name <- factor(df_bar$name, c("HSC", "CMP", "GMP", "MEP", "PreBNK", "MLP"))

colorScale <- c("CMP" = "#F47D2B","HSC" = "#F37B7D","MEP"="#C06CAB","PreBNK"="#272E6A","GMP"="#89C75F","MLP"="#0C727C")

tiff(paste0(OUT_DIR, "numRegulons_p0.05_bargraph.tiff"), units="in", width=5, height=3, res=400)
ggplot(df_bar, aes(x=name, y=count)) + geom_bar(stat="identity", aes(fill=name, color=name)) + theme_classic() + theme(legend.position = "none") + labs(x="", y="Regulons (p < 0.05)")+
  scale_fill_manual(values=colorScale) + theme(axis.text.x = element_text(angle = 45)) + scale_color_manual(values=colorScale)+
  theme(axis.title.y = element_text(size=15))+ theme(axis.text.x = element_text(size=13))+ theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()
