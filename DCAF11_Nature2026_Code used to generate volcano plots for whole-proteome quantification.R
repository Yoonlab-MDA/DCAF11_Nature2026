library(ggplot2)
library(dplyr)
library(tibble)
library(tidyverse)
library(GGally)
library(ggrepel)
library(dplyr)


#Set your working directory
setwd("/Users/hyoon3/Library/CloudStorage/Dropbox/DCAF11 DDX18 paper folder_shared/DCAF11_DDX18/4_Assay data/4_4_Proteomics/wp-esf_637 (1)/260327/")


#Import your file using pathname
Data=read.csv("/Users/hyoon3/Library/CloudStorage/Dropbox/DCAF11 DDX18 paper folder_shared/DCAF11_DDX18/4_Assay data/4_4_Proteomics/wp-esf_637 (1)/wp-esf_637_esf4_diann_directdia_20250207v3-code_rs24k_3c3r_20250208_Limma_output_Pomalidomide_vs_DMSO.csv")

# Set your threshold for your hit
Fold.change <- -0.58
Pval <- -3
Data$Significant <- ifelse(Data$logFC < Fold.change & log10(Data$P.Value) < Pval, "Hit", "Noise")

# Add custom plotting groups
Data$PlotGroup <- Data$Significant
Data$PlotGroup[Data$Gene.Symbol %in% c("ELOC")] <- "BlueHit"
#Data$PlotGroup[Data$Gene.Symbol == "HMGB1"] <- "PurpleHit"

# Change your title here
Experiment.name <- "Pom Whole Proteome"

plot1 <- 
  ggplot(Data, aes(x = -log10(P.Value), y = logFC)) + 
  theme_classic() +
  
  # Noise points (back)
  geom_point(
    data = filter(Data, PlotGroup == "Noise"),
    aes(color = PlotGroup, fill = PlotGroup, shape = PlotGroup),
    size = 4, stroke = 1.2
  ) +
  
  # ELOC, NDUFB11: light blue with border (behind regular hits)
  geom_point(
    data = filter(Data, PlotGroup == "BlueHit"),
    shape = 21,
    size = 4,
    stroke = 1.2,
    color = "#5DADE2",
    fill = scales::alpha("#AED6F1", 0.5)
  ) +
  
  # Regular hit points (in front of BlueHit)
  geom_point(
    data = filter(Data, PlotGroup == "Hit"),
    aes(color = PlotGroup, fill = PlotGroup, shape = PlotGroup),
    size = 4, stroke = 1.2
  ) +
  
 
  # Labels for regular hits + BlueHit
  geom_text_repel(
    data = filter(Data, PlotGroup %in% c("Hit", "BlueHit")),
    aes(label = Gene.Symbol),
    size = 6,
    point.padding = 0.3,
    box.padding = 0.5,
    segment.color = "black",
    segment.size = 0.6,
    min.segment.length = 0,
    force = 2,
    max.overlaps = Inf
  ) +
  

  
  scale_color_manual(values = c(
    "Hit" = "black",
    "Noise" = "#001F1480"
  )) + 
  
  scale_fill_manual(values = c(
    "Hit" = "#F45053",
    "Noise" = "#001F1480"
  )) + 
  
  scale_shape_manual(values = c(
    "Hit" = 21,
    "Noise" = 16
  )) + 
  
  scale_y_reverse() +
  labs(
    x = "-Log10.(P-value)", 
    y = "Log2 - Pomalidomide vs DMSO",
    title = Experiment.name
  ) +
  theme(
    plot.title = element_text(hjust = 0.5),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(size = 1, linetype = "solid", colour = "black")
  ) +
  geom_hline(
    yintercept = c(0.58, -0.58),
    linetype = "dashed", alpha = 0.5, color = "black"
  ) +
  geom_vline(
    xintercept = 3,
    linetype = "dashed", alpha = 0.5, color = "black"
  ) +
  xlim(5, 0) +
  ylim(-3, 3)

plot1

ggsave(paste0(Experiment.name, ".png"), plot1, width = 6, height = 5)