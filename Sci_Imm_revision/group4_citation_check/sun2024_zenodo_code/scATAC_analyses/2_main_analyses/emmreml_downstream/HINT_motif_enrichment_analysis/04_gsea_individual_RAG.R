

#R version 4.1.0 loaded
#rstudio/1.4
#geos/3.9.1 loaded
#hdf5/1.12.0 loaded 
#cmake loaded

.libPaths("/project/lbarreiro/USERS/sarah/Rlibs_new")

library("ggplot2")
library ("DESeq2")
library(statmod)
library(goseq)
library("RColorBrewer")
library(edgeR)
library(devtools)
library(cluster)
library(fgsea)
library(data.table)
library(ggplot2)
library(cowplot)
library(msigdbr)

set.seed(831)


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


  plotEnrichment_custom <- function(pathway, stats,gseaParam=1, ticksSize = 0.25, point_size=0.25, line_color="navy", label=NULL){
    
    rnk <- rank(-stats)
    ord <- order(rnk)
    statsAdj <- stats[ord]
    statsAdj <- sign(statsAdj) * (abs(statsAdj)^gseaParam)
    statsAdj <- statsAdj/max(abs(statsAdj))
    pathway <- unname(as.vector(na.omit(match(pathway, names(statsAdj)))))
    pathway <- sort(pathway)
    
    gseaRes <- calcGseaStat(statsAdj, 
                            selectedStats = pathway, 
                            returnAllExtremes = TRUE)

    bottoms <- gseaRes$bottoms
    tops <- gseaRes$tops
    
    n <- length(statsAdj)
    
    xs <- as.vector(rbind(pathway - 1, pathway))
    ys <- as.vector(rbind(bottoms, tops))
    
    toPlot <- data.frame(x = c(0, xs, n + 1), 
                         y = c(0, ys, 0))
    
    diff <- (max(tops) - min(bottoms))/8
    
    x = y = NULL

    g <- ggplot(toPlot, aes(x = x, y = y)) + 
          geom_line(color = line_color, size=1.75, alpha=0.75) + 
          #geom_point(color = point_color, size = point_size) + 
          geom_hline(yintercept = max(tops), colour = "red", linetype = "dashed") + 
          geom_hline(yintercept = min(bottoms), colour = "red", linetype = "dashed") + 
          geom_hline(yintercept = 0, colour = "black") + 
          geom_segment(data = data.frame(x = pathway), 
            mapping = aes(x = x, y = -diff/2, xend = x, yend = diff/2), size = ticksSize) + 
          theme(panel.border = element_blank(), 
                panel.grid.minor = element_blank()) + 
          labs(x = "Ranked TFs", y = "Enrichment Score")

    if(!is.null(label)){
      
      g <- g + 

      labs(subtitle = label)
    }
    g
  } 
  


setwd("/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/")
IN_DIR <- "/project/lbarreiro/USERS/sarah/HUMAN_BM_PROJECT/BM_CD34_scATAC/Rprojects/ArchR/analysis_FINAL_ALTERNATE/emmreml_downstream/HINT_motif_enrichment_analysis/hint_motif_enrich_results/scatterplot_data_individual/"
OUT_DIR <- "gsea_individual_RAG/"
dir.create(OUT_DIR)

###############################################################
## Get list of TFs that are active in HSC p <= 0.1 criteria ###
###############################################################

clusters <- c("CMP", "GMP")
for(i in 1:length(clusters)){
  name <- clusters[i]
  
  cluster_data <- read.csv(file=paste0(IN_DIR, "0.1_scatterplot_data_", name , ".csv"))
  percent <- cluster_data[,"perc"]
  p_vals <- cluster_data[,"peak_fdr"]
  rank_vals <- p_vals*percent
  names(rank_vals) <- cluster_data[,"name"]
  sorted_ranked_list <- sort(rank_vals, decreasing=FALSE)
  assign(paste0(name,"_table"), sorted_ranked_list)
  print(i)
}

df_list <- list(CMP_table, GMP_table)


for(i in 1:length(df_list)){
  
  regulon_data <- read.csv(file=paste0(IN_DIR, "0.1_scatterplot_data_",clusters[i] ,".csv")) ##using CMP but same for GMP
  HSC_active <- regulon_data$name[which(regulon_data$l == "HSC active")]
  pathways_list <- list()
  pathways_list[[1]] <- as.character(HSC_active)
  names(pathways_list) <- c("HSC_active")
  
  name <- clusters[i]
  sorted_ranked_list <- df_list[[i]]
  
  fgseaRes <- fgsea(pathways = pathways_list, stats= sorted_ranked_list,
                    maxSize  = 500, nperm=100000)
  
  padjusted <- p.adjust(fgseaRes$pval, method='BH')
  
  gsea_result_df<- data.frame(
    pathways <- fgseaRes$pathway,
    pval <- fgseaRes$pval,
    NES <- fgseaRes$NES,
    padj <- fgseaRes$padj
  )
  colnames(gsea_result_df) <- c("pathways","pval","NES", "padj")
  write.csv(gsea_result_df, file=paste0(OUT_DIR,name,"_0.1_gsea_result"))
  print(i)
  
  file_name = paste0(OUT_DIR, paste0("GMP_gsea_0.1.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plotEnrichment(pathways_list[["HSC_active"]],
                       sorted_ranked_list))
  dev.off()
 
}


###############################################################
## Get list of TFs that are active in HSC p <= 0.05 criteria ##
###############################################################

clusters <- c("CMP", "GMP")
for(i in 1:length(clusters)){
  name <- clusters[i]
  
  cluster_data <- read.csv(file=paste0(IN_DIR, "0.05_scatterplot_data_", name , ".csv"))
  percent <- cluster_data[,"perc"]
  p_vals <- cluster_data[,"peak_fdr"]
  rank_vals <- p_vals*percent
  names(rank_vals) <- cluster_data[,"name"]
  sorted_ranked_list <- sort(rank_vals, decreasing=FALSE)
  assign(paste0(name,"_table"), sorted_ranked_list)
  print(i)
}

df_list <- list(CMP_table, GMP_table)

c_colors <- c('CMP'='#FFD147','GMP'='#33CC00')

clusters <- c("CMP", "GMP")
for(i in 1:length(df_list)){
  
  regulon_data <- read.csv(file=paste0(IN_DIR, "0.05_scatterplot_data_",clusters[i] ,".csv")) ##using CMP but same for GMP
  HSC_active <- regulon_data$name[which(regulon_data$l == "HSC active")]
  pathways_list <- list()
  pathways_list[[1]] <- as.character(HSC_active)
  names(pathways_list) <- c("HSC_active")
  
  name <- clusters[i]
  sorted_ranked_list <- df_list[[i]]
  
  fgseaRes <- fgsea(pathways = pathways_list, stats= sorted_ranked_list,
                    maxSize  = 500, nperm=100000)
  
  padjusted <- p.adjust(fgseaRes$pval, method='BH')
  
  gsea_result_df<- data.frame(
    pathways <- fgseaRes$pathway,
    pval <- fgseaRes$pval,
    NES <- fgseaRes$NES,
    padj <- fgseaRes$padj
  )
  colnames(gsea_result_df) <- c("pathways","pval","NES", "padj")
  write.csv(gsea_result_df, file=paste0(OUT_DIR,name,"_0.05_gsea_result"))
  print(i)
  
  file_name = paste0(OUT_DIR, paste0(name,"_gsea_0.05.tiff"))
  tiff(file_name, units="in", width=3, height=3, res=400)
  print(plotEnrichment(pathways_list[["HSC_active"]],
                       sorted_ranked_list))
  dev.off()
  
  p <- plotEnrichment(pathway = pathways_list[["HSC_active"]],
                      stats = sorted_ranked_list)

  plot_label <- paste0(name,"\n",
    "p=" ,format(gsea_result_df$pval, digits=3),"\n",
    "NES=", format(gsea_result_df$NES, digits=3))
  
  p_custom <- plotEnrichment_custom( 
                  pathway = pathways_list[["HSC_active"]],
                  stats = sorted_ranked_list, 
                  line_color=c_colors[name],
                  label=plot_label)

  p_custom_file <- file.path(OUT_DIR, paste0(name,"_custom_gsea_0.05.pdf"))
  ggsave(p_custom_file, p_custom, width=2.5, height=2)
  
}
