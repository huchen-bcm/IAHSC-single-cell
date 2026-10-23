#https://www.kisungyou.com/Rdimtools/reference/nonlinear_PHATE.html

#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

library(here)
library(future)
library(future.apply)
library(scales)
library(tidyverse)
library(broom)
library(data.table)
library(patchwork)
#library(ggpubr)
#library(jcolors)
library(gtools)
library(cowplot)
library(ggbeeswarm)
library(ggridges)
library(ggrepel)
library(ggh4x)
library(colorspace)
library(ggtext)
library(ggrastr)
library(normentR)
library(ggsci)
library(waffle)
library(rstatix)
library(broom)
library(org.Hs.eg.db)
#library(biomaRt)
library(clusterProfiler)
library(limma)
library(edgeR)
library(GSVA)
library(xCell)
library(fgsea)
library(EnrichmentBrowser)
library(phateR)
library(MetBrewer)

# This code is for graphics settings ----------------------------------------------------


dx_colors <- met.brewer(name = "Java")
dx_colors <- dx_colors[c(1:3,5)]
names(dx_colors) <- c("PSC","UC","CD","NC")

main_myColorScale <- c('HSC_a'='#EF5350','HSC_b'='#CC0C00','CMP_a'='#FF7F0E','CMP_b'='#FFD147','CMP_c'='#CC9900',
                       'GMP_a'='#33CC00', 'GMP_b' ='#99CC00', 'MEP_a'='#8C5782', 'MEP_b'='#990080', 'MEP_c'='#BA68C8',
                       'MLP_a'='#46B8DA', 'MLP_b'='#0288D1', 'PreBNK'='#1A0099', 'none'='lightgray')

ct_type_scale <- c('HSC'='#EF5350', 'CMP'='#FFD147','GMP'='#33CC00','MEP'='#990080','MLP'='#0288D1','none'='lightgray')


ct_info <- data.frame(ct=names(main_myColorScale)) %>% 
  mutate(type= gsub("_.*","",ct)) 
ct_info$type[ct_info$type == "PreBNK"] <- "MLP"


vx_colors <- met.brewer(name = "Benedictus", n = 12)[c(1,12)]
names(vx_colors) <- c("BCG","Placebo")



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

update_geom_defaults("point", list(shape = c(16),
                                   color= "darkgray"))
