
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
library(cowplot)

# Here we make an updated version of the gsea heatmap from "04_gsea_padj_heatmap"

# setup --------------------------------
setwd("/scRNA_analyses/2_main_analyses/4_mash_downstream")
OUT_DIR <- c("gsea_padj_plots/")
source("bcf_csf_utils.R")

theme_set(theme_minimal_grid() +
            theme(
              #base_size=2,
              legend.position = "bottom" , 
              text = element_text(family = "Helvetica"), 
              legend.text = element_text(size=8),
              legend.title =  element_text(size=8),
              axis.title = element_text(size=8),
              axis.text.y = element_text(size=8),
              strip.text.x = element_text(size=8, angle=0, vjust = 0.5, hjust = 0.5),
              strip.text.y = element_text(size=8, angle=0, vjust = 0.5, hjust = 0.5),
              axis.text.x = element_text(size=8))
)


# Read in the heatmap object from "04_gsea_padj_heatmap" -----------------------------------
gsea_file <- "gsea_padj_plots/Figure_2B_gsea.rds"
gsea <- readRDS(gsea_file)


### first generate a NES matrix to cluster the enrichments. 
top_gos <- gsea$data %>% 
  group_by(celltype) %>% 
  arrange(padj) 

top_gos <- top_gos %>% 
  #slice_head(n=10) %>% 
  ungroup() %>%
  pull(pathway_name) %>% 
  unique()

NES_mat <- gsea$data %>% 
  filter(pathway_name %in% top_gos) %>%
  mutate(ct_tr=paste0(celltype)) %>% 
  select(ct_tr,pathway_name,NES) %>% 
  spread(pathway_name,NES)

NES_mat <- as.data.frame(NES_mat)
rownames(NES_mat) <- NES_mat[,1]
NES_mat <- NES_mat[,-1]
NES_mat <- as.matrix(NES_mat)
NES_mat <- t(NES_mat)
NES_mat <- NES_mat[which(apply(NES_mat, 1, function(x)sum(is.na(x))) == 0),]
NES_hclust <- hclust(d = dist(NES_mat))

top_gos_final <- rownames(NES_mat)[NES_hclust$order]
top_gos_final_label <- gsub("HALLMARK_","",top_gos_final)
top_gos_final_label <- str_to_title(top_gos_final_label)


### re-generate plotting data adding the order of clusters 
go_bp_dat <- gsea$data %>% 
  filter(pathway_name %in% top_gos) %>% 
  mutate(label = str_to_title(gsub("HALLMARK_","",pathway_name))) %>% 
  mutate(label = factor(label, levels = rev(c(top_gos_final_label)))) %>% 
  #mutate(super_pop=celltype_dat$super_pop[match(cell_types, celltype_dat$work_pops)]) %>% 
  ungroup() 


#go_bp_dat_apop <- go_bp_dat %>% filter(label == "Apoptosis")
#go_bp_dat_gly <- go_bp_dat %>% filter(label == "Estrogen_response_late")

#View(gsea$data[gsea$data$pathway_name == "HALLMARK_APOPTOSIS",])
#View(gsea$data[gsea$data$pathway_name == "HALLMARK_GLYCOLYSIS",])

go_bp_plot <- go_bp_dat %>% 
  group_by(pathway_name) %>% 
  filter(any(padj < 0.1)) %>% 
  ungroup() %>% 
  drop_na() %>% 
  mutate(cluster = celltype ) %>% 
  mutate(cluster = factor(cluster, levels = names(myColorScale))) %>% 
  mutate(cluster_type = ct_info$type[match(cluster,ct_info$ct)]) %>% 
  mutate(cluster_type = factor(cluster_type, levels= names(ct_type_scale))) %>% 
  ggplot(aes(y=label, x=cluster))+
  geom_point(aes(size=-log10(padj), fill=NES, color=padj<0.1), shape=21, stroke=1) +
  scale_fill_met_c(name = "Benedictus",direction = -1) + 
  scale_color_manual(values = c("white","#53565A")) +
  scale_size(range = c(0.75,4.5)) +
  ylab("") + 
  facet_nested(~ cluster_type, scales = "free",space="free",nest_line = TRUE) +
  xlab("") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 7)) + 
  theme(text = element_text(size=9, family = "Helvetica"),   
        strip.text = element_text(size=9, family = "Helvetica",angle=0, hjust = 0.5),
        #axis.text = element_text(size=9, family = "Helvetica"),
        axis.text.x = element_text(size=7,angle=45, hjust = 1), 
        axis.text.y = element_text(size=8), 
        axis.title = element_text(size=8, family = "Helvetica"),
        legend.position = "right", 
        strip.background = element_blank(),
        legend.box="vertical") 


barplot_gsea <- barplot_new + go_bp_plot  + plot_layout(heights = c(0.3,1))

barplot_gsea_file <- file.path(OUT_DIR, 
                               paste0("Fig2_barplot_gsea",".pdf"))
ggsave(filename = barplot_gsea_file, plot = barplot_gsea, height = 5 , width = 5 )


