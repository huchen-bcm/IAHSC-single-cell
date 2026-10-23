##load libraries
library(glmnet)
library(dplyr)
library(stringi)
library(stringr)
library(readxl)
library(data.table)
library(ggplot2)



## Set directory structure
setwd("/scATAC_analyses/2_main_analyses/elastic_net/")
OUT_DIR <- "scatterplots/"
dir.create(OUT_DIR)

######################
### CMP5 IL1B 0.7 ####
######################


## read in predicted values 

read_predicted <- read.table(file=paste0("CMP5_EN_outputs/IL1B_FC_score_predicted_ATAC_CMP5_n15_a0.7.txt"), sep = ",", header=TRUE)

metadata <- read.table(paste0("CMP5_cyt_FC_final.txt"))

which(metadata$Donor != read_predicted$ind )

## plot il1b CMP5_a0.7

test <- cor.test(metadata$IL1B_FC, read_predicted$predicted, method = c("spearman")) 

real_plot_df <- data.frame(
  pred_val <- read_predicted$predicted,
  real_vals <- metadata$IL1B_FC
)
colnames(real_plot_df) <- c("pred_val", "real_vals")

p_real <- ggplot(real_plot_df, aes(x=real_vals, y=pred_val)) + geom_point(shape = 21, size = 3.5, color="grey", fill = "#00468B") + 
  labs(title=paste0("Elastic net IL1B prediction"), x = "Real FC", y=  "Predicted FC") +
  geom_smooth(method=lm, color="grey") + ylim(0.6, 1.8) + theme_linedraw()+theme(axis.title.x=element_text(size=11))+theme(axis.title.y=element_text(size=11))


final_file <- file.path(OUT_DIR,paste0("CMP5_IL1B_a0.7_scatter_predicted.pdf"))
ggsave(filename = final_file, plot = p_real, height = 3, width = 3)



