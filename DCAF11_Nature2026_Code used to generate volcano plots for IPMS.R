library(ggplot2)
library(dplyr)
library(tibble)
library(tidyverse)
library(GGally)
library(ggrepel)
library(dplyr)
  

#Set your working directory
setwd("/Users/hyoon3/Library/CloudStorage/Dropbox/DCAF11 DDX18 paper folder_shared/DCAF11_DDX18/4_Assay data/4_4_Proteomics/IP-MS/ip-esf_361/260326/")


#Import your file using pathname
Data=read.csv("/Users/hyoon3/Library/CloudStorage/Dropbox/DCAF11 DDX18 paper folder_shared/DCAF11_DDX18/4_Assay data/4_4_Proteomics/IP-MS/ip-esf_549/ip-esf_549_73-96-rerun_esf6_diann_directdia_20250707v3-code_globalnorm_3c4r_20250708_rmPCMT1-CRBN_Limma_output_M12_vs_DMSO.csv")

# Set your threshold for your hit
# Copy the column name and add below
Fold.change <- 3
Pval <- -3
Data$Significant <- ifelse(Data$logFC > abs(Fold.change) & log10(Data$P.Value) < Pval, "Hit", "Noise")

# Change your title here
Experiment.name <- "M12_Ip-esf_549"

# Modify the plot to map color to the new column
plot1 <- 
  ggplot(Data, aes(x = logFC, y = log10(P.Value))) +
  geom_text_repel(
    data = filter(Data, Significant == "Hit"), 
    aes(label = Gene.Symbol),
    size = 5, 
    box.padding = unit(0.5, "lines"), 
    point.padding = unit(0.005, "lines"),
    segment.color = "black",
    segment.size = 1,
    force = 2,
    nudge_y = 0.5,
    nudge_x = 0.5
  ) +
  theme_classic() +
  
  # Noise points with 50% transparency and no stroke
  geom_point(
    data = filter(Data, Significant == "Noise"),
    aes(color = Significant, fill = Significant, shape = Significant), 
    shape = 21, 
    size = 4, 
    stroke = 0,  
    fill = "#696969",
    alpha = 0.3  
  ) +
  
  # Hit points with full opacity and black stroke
  geom_point(
    data = filter(Data, Significant == "Hit"),
    aes(color = Significant, fill = Significant, shape = Significant), 
    shape = 21, 
    size = 4, 
    stroke = 1.2,  
    color = "black", 
    alpha = 1
  ) +
  
  # Customizing legend labels and order
  scale_y_reverse() +
  labs(
    x = "Log2 - Drug vs DMSO", 
    y = "Log10.(P-value)",
    title = Experiment.name
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 16),  
    axis.text = element_text(size = 14),   
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, linetype = "solid", colour = "black"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  geom_hline(
    yintercept = -3,
    linetype = "dashed", alpha = 0.5, color = "black"
  ) +
  geom_vline(
    xintercept = c(3, -3),
    linetype = "dashed", alpha = 0.5, color = "black"
  ) +
  xlim(-7, 7) +
  ylim(0, -12)

plot1

ggsave(paste0(Experiment.name, ".png"), plot1, width = 6, height = 5)
