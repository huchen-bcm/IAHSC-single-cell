
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

# This code extends the code from "01_DE_genes_bargraphs" and creates a final version of the bargraph (the one seen in Fig 2B of this paper

# setup -------------------------------------------------------------------
setwd("/scRNA_analyses/2_main_analyses/4_mash_downstream")
OUT_DIR <- c("DE_genes_bargraphs/")
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



# Read in the bargraph that was created in "01_DE_genes_bargraphs" --------------------------------------------
barplot_file <- "/DE_genes_bargraphs/Figure_2A_bargraph.rds"
barplot <- readRDS(barplot_file)


barplot$data$cluster <- factor(barplot$data$cluster, c("HSC_a", "HSC_b", "CMP_a", "CMP_b", "CMP_c", "GMP_a", "GMP_b", "MEP_a", "MEP_b", "MEP_c", "MLP_a", "MLP_b", "PreBNK"))

barplot_new <- barplot$data %>% 
  mutate(fill_col= ifelse(group == "up","BCG","Placebo")) %>% 
  mutate(cluster_type = ct_info$type[match(cluster,ct_info$ct)]) %>% 
  mutate(cluster_type = factor(cluster_type, levels= names(ct_type_scale))) %>% 
  ggplot(aes(x = cluster, y = values)) + 
  geom_bar(aes(fill= fill_col), stat="identity", alpha=0.9) + 
  facet_grid(~cluster_type, scales = "free",space = "free") + 
  #scale_fill_met_d(name = "Benedictus") + 
  scale_fill_manual(values = vx_colors) + 
  ylab("Number of DR genes (lfsr < 0.01)") + 
  xlab("Cell type") + 
  theme(text = element_text(size=8, family = "Helvetica"),   
        strip.text = element_text(size=8, family = "Helvetica",angle=0, hjust = 0.5),
        #axis.text = element_text(size=9, family = "Helvetica"),
        axis.text.x = element_text(size=7,angle=45, hjust = 1), 
        axis.text.y = element_text(size=8), 
        strip.background = element_blank(),
        legend.position="none") 


#transposed 
t_barplot_new <- barplot$data %>% 
  mutate(fill_col= ifelse(group == "up","BCG","Placebo")) %>% 
  mutate(cluster_type = ct_info$type[match(cluster,ct_info$ct)]) %>% 
  mutate(cluster_type = factor(cluster_type, levels= names(ct_type_scale))) %>% 
  ggplot(aes(y = cluster, x = values)) + 
  geom_bar(aes(fill= fill_col), stat="identity", alpha=0.9) + 
  facet_grid(cluster_type~., scales = "free",space = "free") + 
  #scale_fill_met_d(name = "Benedictus") + 
  scale_fill_manual(values = vx_colors) + 
  xlab("N DRG (lfsr < 0.01)") + 
  ylab("Cell type") + 
  theme(text = element_text(size=8, family = "Helvetica"),   
        strip.text = element_text(size=8, family = "Helvetica",angle=0, hjust = 0.5),
        #axis.text = element_text(size=9, family = "Helvetica"),
        axis.text.x = element_text(size=7,angle=45, hjust = 1), 
        axis.text.y = element_text(size=8), 
        strip.background = element_blank(),
        legend.position="none") + scale_y_discrete(limits=rev)

t_barplot_file <- file.path(OUT_DIR, 
                            paste0("Fig2_t_barplot",".pdf"))
ggsave(filename = t_barplot_file, plot = t_barplot_new, height = 2.5 , width = 1.8)
