library(tidyverse)
library(GGally)
library(ggrepel)
library(dplyr)

# Input file has to have the following structure: Col1 - Gene Name, Col2 - sgRNA sequence, remaining columns should containt read number data. 
# Colnames should correspond to the screening nomancluture
split_table_function <- function (main_path, master_input_file_name, screen_name){
  input_data <- read_csv(file.path(paste0(main_path, master_input_file_name)))
  colnames(input_data)[c(1,2)] <- c("Gene", "Sequence")
  colnames_for_single_comparison <- colnames(input_data) %>% str_subset(screen_name)
  selected_input_data <- input_data %>% select(Gene, Sequence, one_of(colnames_for_single_comparison))
  return(selected_input_data)
}
#split_output_table <- split_table_function(main_path, master_input_file_name, screen_name)
#write_csv(split_output_table, paste0(main_path, "split_table/", screen_name, ".csv"))

long_data_transformation_function <- function(wide_selected_input_data, screen_name){
  long_selected_input_data <- wide_selected_input_data %>% 
    gather(key="Experiment", value="Read_numbers", -Gene, -Sequence) %>% 
    separate(Experiment, into=c("Experiment",  "Cell_line", "Drug", "Gate", "Replicate"), sep="_")
  return(long_selected_input_data)}
normalization_function <- function(long_selected_input_data) {
  #substitute 0 and NA with 1
  long_selected_input_data [long_selected_input_data == 0 ] <- 1
  long_selected_input_data [is.na(long_selected_input_data)] <- 1
  
  #Read normalization (Every read number is divided by the mean of the reads within single sample followed by multiplication of average read number for the whole experiment)
  normalized_reads <- long_selected_input_data %>% 
    group_by(Experiment, Cell_line, Drug, Gate, Replicate) %>% 
    mutate(Norm_read_numbers = Read_numbers/mean(Read_numbers)) %>% 
    ungroup() %>% 
    mutate(Norm_read_numbers = Norm_read_numbers * mean(Read_numbers))
  
  return(normalized_reads)
}
wide_normalized_reads_function <- function(normalized_reads){
  normalized_reads_wide <- normalized_reads %>%
    select(-Read_numbers) %>% 
    unite(col = "Experiment", Experiment, Cell_line, Drug, Gate, Replicate, sep="_") %>% 
    spread(key = Experiment, value = Norm_read_numbers)
  return(normalized_reads_wide)
}
read_correlation_function <- function(normalized_reads_wide){
  read_count_log <- normalized_reads_wide %>% select_if(is.numeric)  %>% as.matrix() %>% log() %>% as.data.frame()
  correlation_key_word_plot <- 
    ggpairs(read_count_log, columns = 1:(ncol(read_count_log)),
            upper = list(continuous = wrap("density", alpha = 0.8, col = "steelblue"), combo = "box_no_facet"),
            lower = list(continuous = wrap("points", alpha = 0.2, col = "black", shape = 16, size = 0.8), combo = wrap("dot_no_facet", alpha = 0.4))) +
    theme_classic(base_size = 6) +
    theme(text = element_text(family = "Helvetica"),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(0.2, "cm"))
  return(correlation_key_word_plot)
}


# Calculates the ratio of gate 1 to gate 2 for each guide RNA and ranks them based on the ratio. Replicates are combined by summing up the ranks 
# and by calculating the median ratio for each guide. Finally, guides are merged to genes by calculating the medium ranks and ratios.
rank_guides = function(df_norm, gate_list, replicate_list, screen_name) {
  ranked_df <- df_norm %>%
    filter(Replicate %in% replicate_list,
           Gate %in% gate_list) %>%
    # Combine gates
    group_by(Gene, Sequence, Replicate) %>%
    summarise(Ratio = Norm_read_numbers[Gate == gate_list[1]]/Norm_read_numbers[Gate == gate_list[2]]) %>%
    ungroup() %>%
    group_by(Replicate) %>%
    mutate(Rank = rank(dplyr::desc(Ratio))) %>%
    ungroup() %>%
    # Combine replicates
    group_by(Gene, Sequence) %>%
    summarize(Rank = sum(Rank), Ratio = median(Ratio)) %>%
    ungroup() %>%
    # Combine sgRNAs
    group_by(Gene) %>%
    summarise(Rank = median(Rank), Ratio = median(Ratio))
  return(ranked_df)
}    

# Simulates a p value distribution
simulate_p_value <- function(total_sgRNAs = 2852, n_iterations = 500, replicate_list = c("Rep1", "Rep2", "Rep3"), n_sgRNAs = 4) {
  null_rank_matrix <- matrix(nrow = total_sgRNAs, ncol = n_iterations)
  
  for (i in 1:n_iterations) {
    null_rank_matrix[,i] <- rowSums(cbind(replicate(length(replicate_list), sample(1:total_sgRNAs, total_sgRNAs, replace=FALSE))))
    
  }
  
  null_rank_df <- data.frame(null_rank_matrix) 
  
  null_rank_df <- null_rank_df %>%
    mutate(ID = rep(1:(nrow(null_rank_df)/n_sgRNAs), length.out = nrow(null_rank_df))) %>%
    group_by(ID) %>%
    summarise_all(median, na.rm = T)
  
  null_rank_vector <- data.matrix(null_rank_df[-1]) %>% as.vector()
  return(null_rank_vector)
}

# Test the calculated ranks against the simulated distribution
calculate_p_value <- function(ranked_df, null_rank_vector, total_sgRNAs = 2852, n_iterations = 500, n_sgRNAs = 4, screen_name, FDR_thr = 0.1) {
  calculate_pVal <- function(i){
    sum(null_rank_vector <= as.numeric(ranked_df$Rank[i]))/(total_sgRNAs*(n_iterations/n_sgRNAs))
  }
  
  p_value <- map_dbl(1:(total_sgRNAs/n_sgRNAs), calculate_pVal)
  
  p_value_df <- data_frame(Gene = ranked_df$Gene, Ratio = ranked_df$Ratio, p_value) %>%
    mutate(p_value = ifelse(p_value < 0.5, p_value*2, (1 - p_value)*2),
           p_value = ifelse(p_value == 0.0, min(p_value[p_value > 0]), p_value)) %>%
    mutate(  FDR = p.adjust(p_value, method = "BH"),
             #q_value = qvalue(p_value)$qvalue,
             hits = ifelse((FDR < FDR_thr), T, F) ) %>%
    arrange(desc(p_value))
  return(p_value_df)
}

# Generate two volcano plots for a selected display of hits or every hit that makes the threshold
volcanoplot_screen <-  function(p_value_df, selected_hits, screen_name, gate_list = c("D", "A"), main_path = main_path) {
  experiment_path <- paste0(main_path, screen_name, "/")
  # Generates a volcano plot where x is the ratio of the sum of counts for each guide targeting a gene and y the p-value
  p = ggplot(p_value_df, aes(x = Ratio, y = -log10(p_value))) +
    geom_point(alpha = 0.4, shape = 16, size = 3, col = "darkgray") +
    scale_x_continuous(trans = "log2", breaks = c(0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32)) +
    scale_y_continuous(limits = c(0, -log10(min(p_value_df$p_value)) + 0.5)) +
    labs(title = screen_name, x = paste("Enrichment in", gate_list[1], "vs", gate_list[2]), y = expression('-log'[10]*'(p-value)')) +
    theme_classic(base_size = 10) +
    theme(
      text = element_text(family = "Helvetica"),
      legend.position = "bottom",
      panel.grid = element_blank(),
      title = element_text(size = 12),
      strip.background=element_blank())
  
  p1 = p + geom_point(data = filter(p_value_df, hits == T), alpha = 1, shape = 21, size = 1.5, fill = "#F45053") +
    geom_text_repel(
      data = filter(p_value_df, hits == T),
      aes(label = Gene),
      size = 2,
      box.padding = unit(0.3, "lines"),
      point.padding = unit(0.3, "lines"))
  
  p2 = p + geom_point(data = filter(p_value_df, Gene %in% selected_hits), alpha = 1, shape = 21, size = 3, fill = "#F45053") +
    geom_text_repel(
      data = filter(p_value_df, Gene %in% selected_hits),
      aes(label = Gene),
      size = 3,
      box.padding = unit(0.85, "lines"),
      point.padding = unit(0.4, "lines"))
  
  ggsave(paste0(experiment_path, "Volcano_", screen_name, ".png"), p1, width = 8, height = 8, units = "cm")
  ggsave(paste0(experiment_path, "Volcano_selected_", screen_name, ".png"), p2, width = 8, height = 8, units = "cm")
  
}

# Generates a plot for read numbers of the individual guide RNAs targeting genes that make the threshold
plot_individual_guides = function(df_norm, p_value_df, selected_hits, screen_name, gate_list = c("D", "A"), replicate_list = c("Rep1", "Rep2", "Rep3"), main_path = main_path) {
  experiment_path <- paste0(main_path, screen_name, "/")
  hit_list <- filter(p_value_df, hits == T)$Gene
  
  guides_df <- df_norm %>%
    filter(Replicate %in% replicate_list) %>% 
#           Gate %in% gate_list) %>% USE ALL GATES
    # Combine replicates
    group_by(Gene, Sequence, Gate) %>%
    summarise(Count = mean(Norm_read_numbers)) %>%
    ungroup()
  
  # Arrange guides by mean read count in the high gate
  levels = filter(guides_df, Gate == gate_list[1]) %>%
    group_by(Gene) %>%
    summarize(Mean = mean(Count)) %>%
    arrange(desc(Mean)) %>%
    `$`(Gene) %>%
    as.character()
  
  
  posn_jd = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.4)
  posn_d = position_dodge(width = 0.4)
  
  p1 = ggplot(filter(guides_df, Gene %in% hit_list), aes(x = factor(Gene, levels = levels), y = Count, fill = Gate, group = Gate)) +
    stat_summary(geom = "crossbar", 
                 fun.ymin = mean,
                 fun.ymax = mean,
                 fun.y = mean,
                 width = 0.2,
                 size = 0.3,
                 position = posn_d,
                 fatten = 1,
                 color = "gray",
                 show.legend = F) +
    stat_summary(geom = "errorbar", 
                 fun.ymin = function(x) {mean(x)-sd(x)}, 
                 fun.ymax = function(x) {mean(x)+sd(x)},
                 width = 0.1,
                 size = 0.3,
                 position = posn_d,
                 color = "gray",
                 show.legend = F) +
    geom_point(position = posn_jd, shape = 21, size = 1) +
    labs(y = "Normalized read count") +
    scale_fill_manual("Gate", values = c("green4", "green1", "orange", "red4")) +
    theme_classic(base_size = 6) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(0.2, "cm"),
          text = element_text(family = "Helvetica")
    )
  
  p2 = ggplot(filter(guides_df, Gene %in% selected_hits), aes(x = factor(Gene, levels = levels), y = Count, fill = Gate, group = Gate)) +
    stat_summary(geom = "crossbar", 
                 fun.ymin = mean,
                 fun.ymax = mean,
                 fun.y = mean,
                 width = 0.2,
                 size = 0.3,
                 position = posn_d,
                 fatten = 1,
                 color = "gray",
                 show.legend = F) +
    stat_summary(geom = "errorbar", 
                 fun.ymin = function(x) {mean(x)-sd(x)}, 
                 fun.ymax = function(x) {mean(x)+sd(x)},
                 width = 0.1,
                 size = 0.3,
                 position = posn_d,
                 color = "gray",
                 show.legend = F) +
    geom_point(position = posn_jd, shape = 21, size = 1) +
    labs(y = "Normalized read count") +
    scale_fill_manual("Gate", values = c("green4", "green1", "orange", "red4")) +
    theme_classic(base_size = 6) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(0.2, "cm"),
          text = element_text(family = "Helvetica")
          )
  
  
  
  ggsave(paste0(experiment_path, "Individual_guides_", screen_name, ".png"), p1, width = 8, height = 8, units = "cm")
  ggsave(paste0(experiment_path, "Individual_guides_selected_", screen_name, ".png"), p2, width = 8, height = 8, units = "cm")
}



# Wrapper function that runs all previous functions in sequence
run_pipeline <- function(main_path, input_file_name, screen_name, gate_list = c("D", "A"), replicate_list = c("Rep1", "Rep2", "Rep3", "Rep4"), 
                         n_iterations = 10, n_sgRNAs = 4, 
                         selected_hits, 
                         FDR_thr = 0.1) {
  
  experiment_path <- paste0(main_path, screen_name, "/")
  if (!dir.exists(experiment_path)){dir.create(experiment_path)}
  input_df <- read_csv(paste0(experiment_path, input_file_name))

  output_function_1 <- long_data_transformation_function( wide_selected_input_data = input_df, screen_name)
  total_sgRNAs <- unique(output_function_1$Gene) %>% length() * n_sgRNAs
  output_function_2 <- normalization_function (output_function_1)
  write_csv(output_function_2, paste0(experiment_path, "Normalized_reads_long", screen_name, ".csv"))
  output_function_3 <- wide_normalized_reads_function  (output_function_2)
  write_csv(output_function_3, paste0(experiment_path, "Normalized_reads_", screen_name, ".csv"))
  #correlation_plots <- read_correlation_function(output_function_3)
  #ggsave(paste0(experiment_path, "Correlation_", screen_name, ".png"), correlation_plots, width = 10, height = 10)
  
  output_function_4 <- rank_guides(df_norm = output_function_2, gate_list = gate_list, replicate_list = replicate_list, screen_name = screen_name)
  write_csv(output_function_4, paste0(experiment_path, "Ranked_", screen_name, ".csv"))
  output_function_5 <- simulate_p_value(total_sgRNAs = total_sgRNAs, n_iterations = n_iterations, replicate_list = replicate_list, n_sgRNAs = n_sgRNAs)
  output_function_6 <- calculate_p_value(ranked_df = output_function_4, null_rank_vector = output_function_5, total_sgRNAs = total_sgRNAs, n_iterations = n_iterations, n_sgRNAs = n_sgRNAs, screen_name = screen_name, FDR_thr = FDR_thr)
  write_csv(output_function_6, paste0(experiment_path, "Pval_", screen_name, ".csv"))
  
  volcanoplot_screen(p_value_df = output_function_6, selected_hits = selected_hits, screen_name = screen_name, gate_list = gate_list, main_path = main_path)
  plot_individual_guides(df_norm = output_function_2, p_value_df = output_function_6, selected_hits = selected_hits, screen_name = screen_name, gate_list = gate_list, replicate_list = replicate_list, main_path = main_path)
  
}

main_path <- "/Users/hyoon/Partners HealthCare Dropbox/Hojong Yoon/1_Publication/2025/DCAF11_DDX18/4_Assay data/4_6_CRISPR Screens/Bison/"
master_input_file_name <- "BISON.count.csv"

input_file <- read_csv(paste0(main_path , master_input_file_name))
input_file %>% colnames()

screen_name_list <- paste0(str_split(colnames(input_file)[-(1:2)], "_", simplify = TRUE)[,1], "_",
                           str_split(colnames(input_file)[-(1:2)], "_", simplify = TRUE)[,2], "_",
                           str_split(colnames(input_file)[-(1:2)], "_", simplify = TRUE)[,3]) %>% 
  unique()

#screen_name_list <- c("DepDeg.Set2.4.HEK293T.Cas9_H2.BAMBI_No.Treatment")
# screen_name_list <- c("DepDeg.Set1.5.U937.Cas9_PRDM1", "DepDeg.Set1.5.U937.Cas9_SPDEF_NoTreatment",
#                       "DepDeg.Set1.5.U937.Cas9_NEUROD1", "DepDeg.Set1.5.U937.Cas9_HRAS",
#                       "DepDeg.Set1.5.U937.Cas9_SQSTM1", "DepDeg.Set1.5.U937.Cas9_TP63",
#                       "DepDeg.Set1.5.U937.Cas9_XBP1", "DepDeg.Set1.6.U937.Cas9_HNF4A",
#                       "DepDeg.Set1.6.U937.Cas9_EBF1", "DepDeg.Set1.6.U937.Cas9_NFE2L2",
#                       "DepDeg.Set1.6.U937.Cas9_SREBF1", "DepDeg.Set1.6.U937.Cas9_RPP25L",
#                       "DepDeg.Set1.6.U937.Cas9_HERPUD1", "DepDeg.Set1.6.U937.Cas9_UBIAD1",
#                       "DepDeg.Set1.6.U937.Cas9_PIM2", "DepDeg.Set1.6.U937.Cas9_TTC7A",
#                       "DepDeg.Set1.7.U937.Cas9_FBXO7", "DepDeg.Set1.7.U937.Cas9_ATP6V0A2",
#                       "DepDeg.Set1.7.U937.Cas9_ITPK1", "DepDeg.Set1.7.U937.Cas9_WT1",
#                       "DepDeg.Set1.7.U937.Cas9_JUP", "DepDeg.Set1.7.U937.Cas9_FDFT1",
#                       "DepDeg.Set1.7.U937.Cas9_LMO2", "DepDeg.Set1.7.U937.Cas9_GRHL2",
#                       "DepDeg.Set1.7.U937.Cas9_SNAI2", "epDeg.Set1.7.U937.Cas9_ATP6V0E1",
#                       "DepDeg.Set1.8.U937.Cas9_MYB", "DepDeg.Set1.8.U937.Cas9_KCNK13",
#                       "DepDeg.Set1.8.U937.Cas9_FGFR3", "DepDeg.Set1.8.U937.Cas9_MARCHF5", 
#                       "DepDeg.Set1.8.U937.Cas9_ASNS", "DepDeg.Set1.8.U937.Cas9_FIS1", 
#                       "DepDeg.Set1.8.U937.Cas9_ATP1B3")

for (screen_name in screen_name_list){
  

split_output_table <- split_table_function(main_path, master_input_file_name, screen_name)# %>% 
  #filter(CCND1_Hep3b_DMSO_A_Rep1>50) %>% 
  #group_by(Gene) %>% 
  #mutate(count = n()) %>% 
  #filter(count>2)
ifelse(!dir.exists(paste0(main_path, screen_name)), dir.create(paste0(main_path, screen_name)), FALSE)
write_csv(split_output_table, paste0(main_path, screen_name, "/", screen_name, ".csv"))

#input_file %>% select(-Name, -Library.Sequence) %>% colnames() %>% str_extract(pattern = "[a-z]{2}.*_Rep.*")
#split_output_table %>% select(-Gene, -Sequence) %>% colnames() %>% str_extract(pattern = "(\\w){1}+(?=_Rep)")
replicate_list = str_split(colnames(split_output_table)[-(1:2)], "_", simplify = TRUE)[,5] %>% unique()


run_pipeline(main_path=paste0(main_path), 
                         input_file_name = paste0(screen_name, ".csv"),
                         screen_name = screen_name, 
                         gate_list = c("D","A"), 
                         replicate_list = replicate_list,  #c("Rep1", ) , "Rep2"
                         n_iterations = 100, n_sgRNAs = 4, 
                         selected_hits = c("DCAF11", "UBE2G1", "CUL4A", "COPS5"),
                         FDR_thr = 0.001)}
