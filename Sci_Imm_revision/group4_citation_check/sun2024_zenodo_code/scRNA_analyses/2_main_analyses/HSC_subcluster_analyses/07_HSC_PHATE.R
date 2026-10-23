
# Here we make PHATE maps of the HSCs 

# setup ----------------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/HSC_subcluster_analyses/")
source("bcf_csf_utils.R")
out <- "02_HSC_PHATE/plots/"
dir.create(out,recursive = T)


myColorScale <- c('0'='#EF5350','1'='#CC0C00','2'='#FF7F0E','3'='#FFD147','4'='#CC9900',
                  '5'='#33CC00', '6' ='#99CC00', '7'='#8C5782', '8'='#990080', '9c'='#BA68C8')
names(myColorScale) <- paste0("c",names(myColorScale))

library(scales)
col = "#99CC00"
show_col(col)

#### 1 Create a a HSC specific UMAP to emphasize the HSC clusters in the PHATE clusters 
pdata_file <- file.path("../FINAL_phates_meta_knn7_gamma0.5_t10.rds")
pdata <- readRDS(pdata_file)

hsc_phate <- pdata %>% 
  mutate(hsc= ifelse(grepl("^HSC",as.character(clust_names)),as.character(clust_names),"none")) %>% 
  mutate(hsc = factor(hsc, levels = c("HSC_a","HSC_b","none"))) %>% 
  ggplot(aes(x= PHATE1, y=PHATE2)) +
  geom_point_rast(aes(size=hsc, alpha=hsc, color=hsc)) +
  scale_alpha_manual(values = c(0.8,0.8,0.2)) + 
  scale_size_manual(values = c(0.9,0.9,0.4)) + 
  scale_color_manual(values = main_myColorScale) +
  theme(axis.line = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.text = element_blank(),
        axis.ticks  = element_blank(),
        legend.position = "none")

hsc_phate_file <- file.path(out, 
                             paste0("FINAL_hsc_phate",".pdf"))
#paste0("phates_", "knn",i_params[1],"_gamma",i_params[2],".pdf"))
ggsave(filename = hsc_phate_file, plot = hsc_phate, height = 2 , width = 2)


##### HSC only PHATE
library(reticulate)
use_python("/opt/anaconda3/bin/python")

hsc_file <- file.path("1_UMAP_label_clusters/HSC_final_UNIQUE_IDS.rds")

hsc_sub <- readRDS(hsc_file)

exp <- t(hsc_sub@assays$integrated@scale.data)
hsc_PHATE <- phate(exp, ndim = 2, knn=6, t=10, gamma=1)

pdata_hsc <- pdata %>% filter(grepl("^HSC",as.character(clust_names)))

phate_mb <- as.data.frame(hsc_PHATE$embedding)
colnames(phate_mb) <- c("hsc_PHATE1","hsc_PHATE2")
#plot(seur_PHATE)

hsc_pdata <-  cbind(pdata_hsc,
  phate_mb[rownames(pdata_hsc),c("hsc_PHATE1","hsc_PHATE2")],
  hsc_clusters=hsc_sub@meta.data[rownames(pdata_hsc),"seurat_clusters"]
) %>% 
  mutate(hsc_clusters_lab=paste0("c",hsc_clusters)) 


##### save PHATE data
hsc_pdata_file <- file.path("sc_pdata_phate_ndim2_knn6_t10_gamma1.csv")
write.csv(hsc_pdata,file = hsc_pdata_file)


##### if PHATE data for HSC avail then load this 
hsc_pdata <- read.csv(hsc_pdata_file)


##### MPP score 
hsc_sub <- readRDS(file=paste0("1_UMAP_label_clusters/HSC_final_UNIQUE_IDS.rds"))

MPP_average_expression <- data.frame(AverageExpression(hsc_sub, 
                                                       assays='SCT',
                                                       features = c("THY1", "PTPRC", "ITGA6"), 
                                                       group.by="seurat_clusters")[[1]])

MPP_scaled_data <- data.frame(scale(t(MPP_average_expression), scale=TRUE, center=TRUE))
MPP_scaled_data$average <- rowMeans(MPP_scaled_data)
row.names(MPP_scaled_data) <- str_replace(row.names(MPP_scaled_data),"X","")
MPP_scaled_data_ordered <- MPP_scaled_data[order(MPP_scaled_data$average, decreasing = TRUE),]
write.csv(MPP_scaled_data_ordered, file=paste0("MPP_cluster_data.csv"))


##### GSEA results 
target_pathways <- c("HALLMARK_APOPTOSIS", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
                     "HALLMARK_IL2_STAT5_SIGNALING", "HALLMARK_COMPLEMENT", "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY", "HALLMARK_HYPOXIA",
                     "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_CHOLESTEROL_HOMEOSTASIS","HALLMARK_KRAS_SIGNALING_UP","HALLMARK_KRAS_SIGNALING_DN")


#gsea_res <- list.files("9_gsea",full.names = T)
gsea_res <- list.files("9_gsea",full.names = T)
gsea_res <- gsea_res[grep("^NOT",basename(gsea_res))]
gsea_res <- gsea_res[grep("no_correlations$",basename(gsea_res))]
names(gsea_res) <- gsub("NOT_SORTED_gsea_result_","",gsub("_no_correlations$","",basename(gsea_res)))

gsea_list <- lapply(gsea_res, function(x){
  t <- read.csv(x)
  t %>% 
    filter(pathways %in% target_pathways)
})

gsea_p_cluster <- bind_rows(gsea_list, .id = "hsc_clusters_lab") %>% 
  select(hsc_clusters_lab,pathways,padj) %>% 
  spread(pathways,padj)


##### Plotting all PHATES 
main_clusters <- hsc_pdata %>% 
  droplevels() %>% 
  ggplot(aes(x= hsc_PHATE1, y=hsc_PHATE2)) +
  geom_point_rast(aes(color=clust_names),
                  size=0.8, alpha=0.7) +
  scale_color_manual(values=main_myColorScale) +
  theme(axis.line = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.text = element_blank(),
        axis.ticks  = element_blank(), 
        legend.position = "bottom")

ggsave(filename = file.path(out, paste0("main_clusters_phate.pdf")), 
       plot = main_clusters, 
       height = 3.25 , width = 3)



clusters <- hsc_pdata %>% 
  droplevels() %>% 
  ggplot(aes(x= hsc_PHATE1, y=hsc_PHATE2)) +
  geom_point_rast(aes(color=hsc_clusters_lab),
                  size=0.8, alpha=0.7) +
  scale_color_manual(values=myColorScale) +
  theme(axis.line = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.text = element_blank(),
        axis.ticks  = element_blank(), 
        legend.position = "bottom")

ggsave(filename = file.path(out, paste0("hsc_clusters_phate.pdf")), 
       plot = clusters, 
       height = 3.25 , width = 3)

stem_score <- hsc_pdata %>% 
  droplevels() %>% 
  mutate(stem_score=MPP_scaled_data$average[match(hsc_clusters,rownames(MPP_scaled_data))]) %>% 
  ggplot(aes(x= hsc_PHATE1, y=hsc_PHATE2)) +
  geom_point_rast(aes(color=stem_score),
                  size=0.8, alpha=0.7) +
  scale_color_met_c(name = "Tam") +
  theme(axis.line = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.text = element_blank(),
        axis.ticks  = element_blank(), 
        legend.position = "bottom")

ggsave(filename = file.path(out, paste0("hsc_STEMscore_phate.pdf")), 
       plot = stem_score, 
       height = 3.25 , width = 3)

final_paths <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_HYPOXIA","HALLMARK_IL2_STAT5_SIGNALING")

final_path_test <- c("HALLMARK_TNFA_SIGNALING_VIA_NFKB")
gsea_phates <- hsc_pdata %>% 
  left_join(gsea_p_cluster[,c("hsc_clusters_lab",final_path_test)],by = "hsc_clusters_lab") %>% 
  select(one_of(final_path_test), hsc_clusters_lab, hsc_PHATE2, hsc_PHATE1) %>% 
  pivot_longer(one_of(final_path_test)) %>% 
  drop_na(value) %>% 
  mutate(path_sig = value < 0.05) %>% 
  mutate(name = gsub("HALLMARK_","",name)) %>% 
  ggplot(aes(x= hsc_PHATE1, y=hsc_PHATE2)) +
  geom_point_rast(aes(color=-log10(value), size=path_sig, alpha=path_sig)) +
  scale_color_met_c(name = "Hokusai3")  + 
  scale_alpha_manual(values = c(0.3, 0.7)) + 
  scale_size_manual(values = c(0.3,0.7)) +
  theme(axis.line = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.text = element_blank(),
        axis.ticks  = element_blank(), 
        legend.position = "bottom")

ggsave(filename = file.path(out, paste0("hsc_GSEA_phate.pdf")), 
       plot = gsea_phates, 
       height = 3.25 , width = 7)



phates_s <- clusters + stem_score + gsea_phates + 
  plot_layout(nrow = 1, widths = c(1,1,3)) & 
  theme(axis.line = element_blank(),
        plot.background = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(), 
        axis.text = element_blank(),
        axis.ticks  = element_blank(), 
        legend.position = "bottom")



i_phates_s_file <- file.path(out, 
                             paste0("FINAL_HSC_phates_", "knn","6","_gamma","1","_t10",".pdf"))
ggsave(filename = i_phates_s_file, plot = phates_s, height = 3.25 , width = 12.5 )
#paste0("phates_", "knn",i_params[1],"_gamma",i_params[2],".pdf"))

