
#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

library("ggplot2")
library ("DESeq2")
library(statmod)
library("RColorBrewer")
library(edgeR)
library(devtools)
library(cluster)
library(ggrepel)
library(cowplot)
library(dplyr)


setwd("/scATAC_analyses/2_main_analyses/emmreml_downstream")
out_dir <- "DA_peaks_barplot/"

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


barplot_file <- "DA_peaks_barplot/Figure_4C_barplot.rds"
barplot <- readRDS(barplot_file)


myColorScale <- c("CMP1" = "#FEE500","CMP2"="#F47D2B", "HSC1" = "#D51F26", "HSC2"="#F37B7D","MEP1"="#89288F","PreBNK"="#272E6A","GMP1"="#208A42","MLP"="#0C727C","MEP2"="#C06CAB", 
                  "MEP3"="#E6C2DC", "CMP3"="#D24B27","MEP4"="#9983BD", "CMP4"="#D8A767","CMP5"="#7E1416","GMP2"="#89C75F","GMP3"="#3BBCA8")

ct_info <- data.frame(ct=names(myColorScale)) %>% 
  mutate(type= gsub('[0-9]+', '',ct)) 

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
        legend.position="none") 

t_barplot_file <- file.path(out_dir, 
                            paste0("Fig4_ATAC_barplot",".pdf"))
ggsave(filename = t_barplot_file, plot = t_barplot_new, height = 3.3 , width = 2.5)
