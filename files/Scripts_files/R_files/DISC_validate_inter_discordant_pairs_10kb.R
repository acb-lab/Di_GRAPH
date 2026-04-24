#### R script to process inter-chromosomal discordant pairs (reads 75nt) for 10kb validation
### AFTER BLAST cross-validation
### Loop for each strain/sample/experiment in MYWD
### 24/04/2026 - Lydia

log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}

# -------------------------
# Read command line args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 7) {
  log_step("ERROR: Incorrect number of arguments")
  log_step("Usage: script.R <strain> <sample> <experiment> <file1> <fileoriginal> <filecontrol> <category_path>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

strain <- args[1]
sample   <- args[2]
experiment <- args[3]
path_to_file_1   <- args[4]

path_to_file_original <- args[5]
path_to_file_control <- args[6]
category_path   <- args[7]
strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))
log_step(paste("SAMPLE:", sample))
log_step(paste("EXPERIMENT:", experiment))
log_step(paste("file1:", path_to_file_1))
log_step(paste("fileoriginal:", path_to_file_original))
log_step(paste("filecontrol:", path_to_file_control))

# load libraries

library(readr)
library(stringr)
library(extrafont)
library(svglite)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)


path_to_categories_pairs_file <- file.path(category_path, "PMV_categories_pairs.tsv")
path_to_features_pairs_file <- file.path(category_path, "PMV_features_pairs.tsv")


#path_to_file_1
#path_to_file_2
#path_to_file_3
#path_to_file_original


# Functions

process_blast_filtered_files <- function (blast_filtered_file_path) {
  lines <- readLines(blast_filtered_file_path)
  split_lines <- strsplit(lines, "\t")
  parsed <- lapply(split_lines, function(x) {
    base <- x[1:9]
    
    extra <- x[-(1:9)]
    
    # Split extras into groups of 8
    if (length(extra) > 0) {
      groups <- split(extra, ceiling(seq_along(extra) / 8))
    } else {
      groups <- list()
    }
    
    list(
      base = base,
      extra = groups
    )
  })
  
  base_df <- do.call(rbind, lapply(parsed, function(x) x[["base"]]))
  base_df <- as.data.frame(base_df, stringsAsFactors = FALSE) %>% 
    rename("pair_group" = !!names(.[1]), "Feature_name_B" = !!names(.[2]), "hits" = !!names(.[3]),
           "mismatch" = !!names(.[4]), "length" = !!names(.[5]),
           "start_pos" = !!names(.[6]),  "end_pos" = !!names(.[7]),
           "seq_A" = !!names(.[8]),  "seq_chr" = !!names(.[9]))
  
  all_groups <- do.call(
    c,
    lapply(parsed, function(x) x$extra)
  )
  
  if (length(all_groups) > 0) {
    groups_df <- do.call(rbind, all_groups)
    groups_df <- as.data.frame(groups_df, stringsAsFactors = FALSE) %>%
      rename("pair_group" = !!names(.[1]), "Feature_name_B" = !!names(.[2]), 
             "mismatch" = !!names(.[3]), "length" = !!names(.[4]),
             "start_pos" = !!names(.[5]),  "end_pos" = !!names(.[6]),
             "seq_A" = !!names(.[7]),  "seq_chr" = !!names(.[8]))
    
    complete_df <- bind_rows(base_df, groups_df)
  } else {
    complete_df <- base_df
  }
  
  blast_filtered_df <- complete_df %>% group_by(pair_group) %>%
    mutate(hits = ifelse(is.na(hits), first(na.omit(hits)), hits)) %>%
    ungroup() %>%
    mutate(across(c(hits, mismatch, length, start_pos, end_pos), as.double)) %>% 
    mutate(valid = ifelse(hits == 0, "yes", "no"),
           start_pos_10kb = ifelse(start_pos < end_pos, start_pos-10000, end_pos-10000),
           end_pos_10kb = ifelse(end_pos > start_pos, end_pos+10000, start_pos+10000))
  
  return(blast_filtered_df)
  
}

get_blast_validated_0hits_complete_df <- function (chromosome_file) {
  
  complete_df <- chromosome_file %>% 
    separate(pair_group, c("pair_group_name", "number"), "-", remove = FALSE) %>% 
    mutate(valid_all = ifelse(valid == "yes", "yes", "no")) %>% 
    group_by(pair_group_name, pair_group) %>%
    mutate(all_yes_T_F = all(valid_all == "yes")) %>%
    group_by(pair_group_name) %>%
    mutate(group_all_yes_T_F = all(all_yes_T_F)) %>%
    ungroup() %>%
    mutate(valid_def = ifelse(group_all_yes_T_F, "yes", "no")) %>% 
    filter(valid_def == "yes") %>% select(pair_group) %>% 
    rename(Read_name_ID = pair_group) 
  
  
  return(complete_df)
  
}


get_blast_validated_10kb_complete_df <- function (chromosome_file, original_file) {
  
  complete_df <- chromosome_file %>% rename(Read_name_ID = pair_group) %>% 
    full_join(., original_file, by = "Read_name_ID") %>% 
    mutate(min_distance = pmin(
      abs(start_pos - Start_pos_B),
      abs(end_pos - Start_pos_B),
      na.rm = FALSE)) %>% 
    group_by(Read_name_ID) %>% arrange(min_distance) %>% 
    slice(1) %>%
    ungroup() %>% 
    select(!c("Feature_name_B.x")) %>% 
    mutate(valid_10kb = ifelse(Start_pos_B < start_pos_10kb | Start_pos_B > end_pos_10kb | hits == 0, "yes", "no")) %>% 
    group_by(pair_group_name) %>%
    mutate(valid_all = if_else(all(valid_10kb == "yes"), "yes", "no")) %>%
    ungroup() %>% 
    filter(valid_all == "yes") %>% select(Read_name_ID) 
  
  
  return(complete_df)
  
}         

chromosome_df <- process_blast_filtered_files(path_to_file_1)


log_step("Finding valid reads...") 
valid_read_name_0hits <- get_blast_validated_0hits_complete_df(chromosome_df)

log_step("Reading original df...") 
original_df <- read_tsv(path_to_file_original, col_names = TRUE) %>% mutate(Read_name_ID = with(., paste0(Read_name, "_", pair_group)))

log_step("Finding valid reads in 10kb area...") 
valid_read_name_10kb <- get_blast_validated_10kb_complete_df(chromosome_df, original_df) %>% unique()

valid_read_name <- bind_rows(valid_read_name_0hits, valid_read_name_10kb)
log_step("Comparing with original file...") 

valid_pairs_df <- semi_join(original_df, valid_read_name)

n_valid_reads <- nrow(valid_pairs_df)
n_original_reads <- nrow(original_df)

log_step("Calculating error rate...")

error_rate <- tibble(strain = strain_name, 
                     sample = sample,
                     experiment = experiment,
                     original_reads = n_original_reads,
                     valid_reads = n_valid_reads) %>% 
  mutate(valid_rate = (valid_reads/original_reads)*100)



# Write the new TSV
log_step("Saving processed dataframe...") 
log_step(paste0("Saving tsv files for: ", strain_name, " ",  sample, " ", experiment, "..."))

write_tsv(valid_pairs_df, file = paste0(strain, "/", sample, "_", experiment, "_inter_discordant_pairs_unique_processed_valid.tsv"))
write_tsv(error_rate, file = paste0(strain, "/", sample, "_", experiment, "_inter_discordant_pairs_unique_processed_valid_error_rate.tsv"))

# Calculate distribution of valid reads
df_categories_pairs <- read_tsv(path_to_categories_pairs_file,col_names = FALSE, show_col_types = FALSE) %>%
  rename("Category_A" = !!names(.[1]), "Category_B" = !!names(.[2]))

# Categories
categories <- unique(df_categories_pairs$Category_A)

# Process_categories
process_categories <- function(categories_i, valid_pairs_df_file, categories_df_file) {
  # Possible combinations
  combinations <- categories_df_file %>%
    filter(Category_A == categories_i) %>%
    distinct(Category_B) %>%
    pull()
  
  # Filter discordant pairs for each Category_A
  valid_pairs_filtered <- filter(valid_pairs_df_file, Category_A == categories_i)
  
  # Calculate number of reads 
  map_dfr(combinations, function(categoryB) {
    n <- valid_pairs_filtered %>% filter(Category_B == categoryB) %>% nrow()
    tibble(Category_A = categories_i, Category_B = categoryB, count = n)
  })
}

# Write the new TSV
log_step("Calculating discordant reads distribution along genomic categories...") 

# Get counts
count_df <- map_dfr(categories, process_categories, valid_pairs_df_file = valid_pairs_df, categories_df_file = df_categories_pairs) %>% 
  mutate(strain = strain_name,
         sample = sample,
         experiment = experiment)
write_tsv(count_df, file = paste0(strain, "/", sample, "_", experiment, "_inter_discordant_pairs_unique_processed_valid_counts.tsv"))


# Set genomic categories order
genomic_categories_order <- c(
  "ORF", "intergenic", "long_terminal_repeat", "transposable_element_gene",
  "LTR_retrotransposon", "tRNA_gene", "rRNA_gene", "ncRNA_gene",
  "snRNA_gene", "snoRNA_gene", "ARS", "centromere", "telomere"
)


log_step("Calculating global inter-chromosomal discordant read pairs distribution...")
# Calculate global inter-chromosomal discordant read pairs distribution 
global_distribution <- count_df %>% 
  group_by(strain, sample, experiment) %>% 
  mutate(global_percentage = (count / sum(count))*100,
         strain_sample_comb = paste0(strain, "_", sample)) %>% 
  replace(is.na(.), 0) %>% 
  group_by(strain, sample, experiment, strain_sample_comb,  Category_A, Category_B) %>%
  summarise(mean_global_percentage = mean(global_percentage, na.rm = TRUE),
            sd_global_percentage = sd(global_percentage, na.rm = TRUE)) %>%
  # Apply genomic_categories_order
  mutate(Category_A = factor(Category_A, levels = genomic_categories_order),
         Category_B = factor(Category_B, levels = genomic_categories_order)) %>%
  arrange(strain, sample, experiment, Category_A, Category_B)

write_tsv(global_distribution, file = paste0(strain, "/", sample, "_", experiment, "_inter_discordant_pairs_unique_processed_valid_global_distribution.tsv"))

log_step("Calculating category inter-chromosomal discordant read pairs distribution...")
# Calculate average per_category inter-chromosomal discordant read pairs distribution 
category_distribution <- count_df %>% 
  group_by(strain, sample, experiment, Category_A) %>% 
  mutate(category_percentage = (count / sum(count))*100,
         strain_sample_comb = paste0(strain, "_", sample)) %>% 
  replace(is.na(.), 0) %>%
  group_by(strain, sample, experiment, strain_sample_comb, Category_A, Category_B) %>%
  summarise(mean_category_percentage = mean(category_percentage, na.rm = TRUE),
            sd_category_percentage = sd(category_percentage, na.rm = TRUE)) %>%
  # Apply genomic_categories_order
  mutate(Category_A = factor(Category_A, levels = genomic_categories_order),
         Category_B = factor(Category_B, levels = genomic_categories_order)) %>%
  arrange(strain, sample, experiment,Category_A, Category_B)

write_tsv(category_distribution, file = paste0(strain, "/", sample, "_", experiment, "_inter_discordant_pairs_unique_processed_valid_category_distribution.tsv"))


# Inter-chromosomal discordant read pairs matrix
df_features_pairs <- read_tsv(path_to_features_pairs_file, col_names = TRUE, show_col_types = FALSE)
control_df <- read_tsv(path_to_file_control,col_names = TRUE) %>% filter(Category_A=="control_norm", Category_B=="control_norm")

#head(all_pairs)

get_discordant_matrix <- function(valid_pairs_df_file, all_pairs_df, control_df_file) {
  posA_info <- all_pairs_df %>%
    select(Feature_name_A, Position_A = Position, Essential_A = Essential) %>%
    distinct()
  
  posB_info <- all_pairs_df %>%
    select(Feature_name_B, Position_B = Position, Essential_B = Essential) %>%
    distinct()
  
  valid_pairs_df_posA <- valid_pairs_df_file %>%
    left_join(posA_info, by = "Feature_name_A")
  
  valid_pairs_df_posAB <- valid_pairs_df_posA %>%
    left_join(posB_info, by = "Feature_name_B")
  
  pair_counts <- valid_pairs_df_posAB %>%
    group_by(strain, sample, experiment, Feature_name_A, Feature_name_B) %>%
    summarise(count = n(), .groups = "drop")
  
  complete_matrix <- left_join(valid_pairs_df_posAB, pair_counts,
                               by = c("Feature_name_A", "Feature_name_B", "strain", "sample", "experiment")) %>%
    distinct()
  
  count_total <- nrow(control_df_file)
  complete_matrix <- complete_matrix %>%
    mutate(count_norm = (count / count_total) * 100) %>% 
    select(!c(pair_group, pair_group_name, number, Feature_name_B_prev, Feature_name_B_next))
  
  return(complete_matrix)
}

log_step("Generating discordant matrix...")
complete_matrix <- get_discordant_matrix(valid_pairs_df, df_features_pairs, control_df)
complete_matrix_noORF_nointergenic <- complete_matrix %>% filter(Category_A != "ORF", Category_B != "ORF", Category_A != "intergenic", Category_B != "intergenic")

log_step("Saving discordant matrix...")
write_tsv(complete_matrix, file = paste0(strain, "/", sample, "_", experiment, "_inter_discordant_pairs_unique_processed_valid_discordant_matrix.tsv"))
write_tsv(complete_matrix_noORF_nointergenic, file = paste0(strain, "/", sample, "_", experiment, "_inter_discordant_pairs_unique_processed_valid_discordant_matrix_noORF_nointergenic.tsv"))

log_step("Plotting discordant matrix...")
matrix_plot <-ggplot(complete_matrix, aes(y = Position_B, x = Position_A)) +
  #geom_hdr(xlim = c(0,14518), ylim = c(14518,0), method = "kde", fill = "brown4") + 
  geom_point(aes(colour = count_norm), alpha = 1, size = 0.8) +
  scale_colour_gradient2(low= "white", mid = "#4ae034", high = "#00641b",
                         midpoint = 1,
                         limits = c (0, 5),
                         oob = scales::squish) + 
  theme_bw(base_family = "Arial") + 
  theme(panel.background = element_blank()) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(legend.position="right") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  guides (colour = guide_colourbar(barwidth = 0.5, barheight = 5,
                                   frame.colour = "black", frame.linewidth = 0.25,
                                   ticks.colour = NA)) + 
  coord_cartesian(xlim = c(0,14518), ylim = c(14518, 0), expand=FALSE) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                                6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                                13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487), position = "top") +
  scale_y_continuous(breaks = c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                                6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                                13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487)) +
  geom_hline(yintercept=c(6570, 13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="solid", color = "black", linewidth=0.05) +
  geom_hline(yintercept=c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                          6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                          13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 144877),
             linetype="dashed", color = "black", linewidth=0.05) +
  geom_vline(xintercept=c(6570, 13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="solid", color = "black", linewidth=0.05) +
  geom_vline(xintercept=c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                          6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                          13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="dashed", color = "black", linewidth=0.05) +
  theme(axis.text.y=element_text(size=0)) +
  theme(axis.text.x=element_text(size=0)) +
  theme(axis.title.x = element_text(size=0)) +
  theme(axis.title.y = element_text(size=0)) +
  theme(axis.ticks = element_blank()) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "solid", linewidth = 0.05) +
  labs(
    title = paste0("Inter_chromosomal discordant matrix - ", strain_name, " - ", sample, " - ", experiment),
    colour = "Freq (%)")

#paste0(subdir, "/plot_", chr_name, "_75nt_", gsub("\\.tsv$", "", suffix), ".svg")
ggsave(
  filename = paste0(strain,"/", "Inter_chromosomal_discordant_matrix_plot_", strain_name, "_", sample, "_", experiment, ".svg"),
  plot = matrix_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)

log_step("Plotting discordant matrix - no ORF no intergenic...")
matrix_plot <-ggplot(complete_matrix_noORF_nointergenic, aes(y = Position_B, x = Position_A)) +
  #geom_hdr(xlim = c(0,14518), ylim = c(14518,0), method = "kde", fill = "brown4") + 
  geom_point(aes(colour = count_norm), alpha = 1, size = 2) +
  scale_colour_gradient2(low= "white", mid = "#4ae034", high = "#00641b",
                         midpoint = 1,
                         limits = c (0, 5),
                         oob = scales::squish) + 
  theme_bw(base_family = "Arial") + 
  theme(panel.background = element_blank()) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(legend.position="right") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  guides (colour = guide_colourbar(barwidth = 0.5, barheight = 5,
                                   frame.colour = "black", frame.linewidth = 0.25,
                                   ticks.colour = NA)) + 
  coord_cartesian(xlim = c(13195,14518), ylim = c(14518, 13195), expand=FALSE) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                                6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                                13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487), position = "top") +
  scale_y_continuous(breaks = c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                                6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                                13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487)) +
  geom_hline(yintercept=c(6570, 13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="solid", color = "black", linewidth=0.05) +
  geom_hline(yintercept=c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                          6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                          13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 144877),
             linetype="dashed", color = "black", linewidth=0.05) +
  geom_vline(xintercept=c(6570, 13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="solid", color = "black", linewidth=0.05) +
  geom_vline(xintercept=c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                          6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                          13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="dashed", color = "black", linewidth=0.05) +
  theme(axis.text.y=element_text(size=0)) +
  theme(axis.text.x=element_text(size=0)) +
  theme(axis.title.x = element_text(size=0)) +
  theme(axis.title.y = element_text(size=0)) +
  theme(axis.ticks = element_blank()) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "solid", linewidth = 0.05) +
  labs(
    title = paste0("Inter_chromosomal discordant matrix no ORF no intergenic - ", strain_name, " - ", sample, " - ", experiment),
    colour = "Freq (%)")

ggsave(
  filename = paste0(strain,"/", "Inter_chromosomal_discordant_matrix_reduced_plot_", strain_name, "_", sample, "_", experiment, ".svg"),
  plot = matrix_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)