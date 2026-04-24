#### R script to process inter-chromosomal discordant pairs (reads 75nt) for 10kb validation 
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

if (length(args) != 6) {
  log_step("ERROR: Incorrect number of arguments")
  log_step("Usage: script.R <strain> <sample> <experiment> <file1> <file2> <category_path>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

strain <- args[1]
sample   <- args[2]
experiment <- args[3]
path_to_file_1   <- args[4]
path_to_file_2 <- args[5]
category_path   <- args[6]
strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))
log_step(paste("SAMPLE:", sample))
log_step(paste("EXPERIMENT:", experiment))
log_step(paste("file1:", path_to_file_1))
log_step(paste("file2:", path_to_file_2))


# load libraries

library(readr)
library(stringr)
library(extrafont)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)



path_to_categories_file <- file.path(category_path, "PMV_categories.tsv")
path_to_categories_near_file <- file.path(category_path, "PMV_categories_near.tsv")

# Load categories dataframe 
# and categories_near_features dataframe
df_categories <- read_tsv(path_to_categories_file,col_names = FALSE, show_col_types = FALSE)
df_categories <- type.convert(df_categories, as.is = TRUE)
df_categories_near_features <- read_tsv(path_to_categories_near_file,col_names = TRUE, show_col_types = FALSE)


chromosomes <- c(
  "CHRI", "CHRII", "CHRIII", "CHRIV", "CHRV", "CHRVI", "CHRVII",
  "CHRVIII", "CHRIX", "CHRX", "CHRXI", "CHRXII", "CHRXIII", 
  "CHRXIV", "CHRXV", "CHRXVI", "CHRDISC")

log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}
process_df_categories_file <- function(df_categories_file, chromosomes_vector) {
  # Rename columns 
  all_df <- df_categories_file %>%
    rename(
      Category = !!names(.)[1],
      Feature_name = !!names(.)[2],
      Start_pos = !!names(.)[3],
      End_pos = !!names(.)[4],
      Chromosome = !!names(.)[5]
    )
  
  # Create a list to store results
  result_list <- list()
  
  for (chrom in chromosomes_vector) {
    # Filter by chromosome
    df_chrom <- filter(all_df, Chromosome == chrom)
    
    # If empty, skip
    if (nrow(df_chrom) == 0) next
    
    # For each unique Feature_name create a dataframe subset
    list_of_dfs <- lapply(unique(df_chrom$Feature_name), function(feat) {
      df_sub <- df_chrom %>%
        filter(Feature_name == feat) %>%
        select(Chromosome, Start_pos, End_pos, Category)
      df_sub$Feature_name <- feat
      df_sub
    })
    
    # Convert factor columns to proper types
    list_of_dfs <- lapply(list_of_dfs, function(x) type.convert(x, as.is = TRUE))
    
    names(list_of_dfs) <- unique(df_chrom$Feature_name)
    
    # Store in result list with chromosome name
    result_list[[chrom]] <- list_of_dfs
  }
  
  return(result_list)
}


filter_by_chromosome_A <- function(discordant_pairs_file, chromosomes_vector, df_categories_list_file) {
  
  chromosome_A_filter_results <- list()
  
  for (i in seq_along(chromosomes)) {
    chrom <- chromosomes[i]
    df_chr <- discordant_pairs_file %>% filter(Chr_name_A == chrom)
    
    list_of_dfs_chr <- df_categories_list_file[[i]]
    
    filtered_data_list <- list()
    
    for (j in seq_along(list_of_dfs_chr)) {
      current_df <- list_of_dfs_chr[[j]]
      
      filtered_data <- df_chr %>%
        filter(Start_pos_A >= current_df$Start_pos, Start_pos_A <= current_df$End_pos)
      
      filtered_data$Feature_name_A <- current_df$Feature_name
      filtered_data$Category_A <- current_df$Category
      filtered_data$Chromosome_A <- current_df$Chromosome
      
      filtered_data_list[[j]] <- filtered_data
    }
    
    chromosome_A_filter_results[[i]] <- bind_rows(filtered_data_list)
  }
  
  chromosome_A_filter_results_complete <- bind_rows(chromosome_A_filter_results) %>% unique()
  return(chromosome_A_filter_results)
}

filter_by_chromosome_B <- function(df_chr_A_filtered_file, chromosomes_vector, df_categories_list_file) {
  
  chromosome_B_filter_results <- list()
  
  for (i in seq_along(chromosomes)) {
    chrom <- chromosomes[i]
    df_chr <- df_chr_A_filtered_file %>% filter(Chr_name_B == chrom)
    
    list_of_dfs_chr <- df_categories_list_file[[i]]
    
    filtered_data_list <- list()
    
    for (j in seq_along(list_of_dfs_chr)) {
      current_df <- list_of_dfs_chr[[j]]
      
      filtered_data <- df_chr %>%
        filter(Start_pos_B >= current_df$Start_pos, Start_pos_B <= current_df$End_pos)
      
      filtered_data$Feature_name_B <- current_df$Feature_name
      filtered_data$Category_B <- current_df$Category
      filtered_data$Chromosome_B <- current_df$Chromosome
      
      filtered_data_list[[j]] <- filtered_data
    }
    
    chromosome_B_filter_results[[i]] <- bind_rows(filtered_data_list)
  }
  
  chromosome_B_filter_results_complete <- bind_rows(chromosome_B_filter_results) %>% unique()
  return(chromosome_B_filter_results)
}

prepare_for_blast <- function(df_chr_AB_filtered_file_X_ID) {
  blast_df <- df_chr_AB_filtered_file_X_ID %>% 
    filter (Category_A != "control") %>% 
    filter(Category_A != "control_norm") %>% 
    filter (Category_B != "control") %>% 
    filter(Category_B != "control_norm") %>% 
    mutate(Read_name_ID = with(., paste0(Read_name, "_", pair_group))) %>% 
    select(Read_name_ID, Sequence, Chromosome_B, Feature_name_B, Feature_name_B_prev, Feature_name_B_next)
  return(blast_df)
}

log_step("Loading categories list...")
# Get categories list
df_categories_list <- process_df_categories_file(df_categories, chromosomes)

log_step(paste0("Processing file: " , path_to_file_1,  "..."))

# Process File_1 : from T0
file_1 <- read_tsv(path_to_file_1, col_names = FALSE, show_col_types = FALSE)

file_1 <- file_1 %>%
  rename(
    Read_name = !!names(.)[1],
    Chr_name_A = !!names(.)[2],
    Start_pos_A = !!names(.)[3],
    Chr_name_B = !!names(.)[4],
    Start_pos_B = !!names(.)[5],
    Sequence = !!names(.)[6]
  )

log_step("Filtering by Feature_name_A...") 


df_chr_A_filtered_file_1 <- filter_by_chromosome_A(file_1, chromosomes, df_categories_list) %>% 
  bind_rows(.)

log_step("Filtering by Feature_name_B...")   

df_chr_AB_filtered_file_1 <- filter_by_chromosome_B(df_chr_A_filtered_file_1, chromosomes, df_categories_list) %>% 
  bind_rows(.) %>% 
  mutate(strain = rep(c(strain_name)), sample = rep(c("T0")), experiment = rep(c(experiment)))

df_chr_AB_filtered_file_1_nocontrols <- df_chr_AB_filtered_file_1 %>%
  filter (Category_A != "control", Category_B != "control" | Category_A != "control_norm", Category_B != "control_norm")

remove(df_chr_A_filtered_file_1)

log_step("Processing for BLAST analysis...") 

df_chr_AB_filtered_file_1_ID <- df_chr_AB_filtered_file_1_nocontrols %>%
  mutate(pair_id = pmin(paste(Feature_name_A, Feature_name_B, sep = "_"), paste(Feature_name_B, Feature_name_A, sep = "_"))) %>%
  group_by(Read_name, pair_id) %>%
  mutate(pair_group = paste0("id", cur_group_id(), "-", row_number())) %>%
  separate(pair_group, c("pair_group_name", "number"), "-", remove = FALSE)  %>% 
  mutate(Read_name_ID = paste0(Read_name, "_", pair_group)) %>%  
  ungroup() %>%
  select(-pair_id) %>% 
  left_join(., df_categories_near_features, by = "Feature_name_B")

control_file_1 <- df_chr_AB_filtered_file_1 %>% 
  filter (Category_A == "control" & Category_B == "control" | Category_A == "control_norm" & Category_B == "control_norm")

blast_df_file_1 <- prepare_for_blast(df_chr_AB_filtered_file_1_ID) %>% select(Read_name_ID, Sequence, Chromosome_B)



# Write the new TSV
log_step("Saving processed dataframes...") 
log_step(paste0("Saving tsv files for: ", strain_name,  " T0 ", experiment, "..."))


write_tsv(df_chr_AB_filtered_file_1_ID, file = paste0(strain, "/", "T0", "_", experiment, "_inter_discordant_pairs_unique_processed.tsv"))
write_tsv(control_file_1, file = paste0(strain, "/", "T0", "_", experiment, "_inter_discordant_pairs_unique_processed_control.tsv"))
write_tsv(blast_df_file_1, file = paste0(strain, "/", "T0", "_", experiment, "_inter_discordant_pairs_unique_processed_blast_chromosome.tsv"))


log_step(paste0("Processing file: " , path_to_file_2,  "..."))
# Process File_2 : from TSG/TLG/TLR

file_2 <- read_tsv(path_to_file_2, col_names = FALSE, show_col_types = FALSE)

file_2 <- file_2 %>%
  rename(
    Read_name = !!names(.)[1],
    Chr_name_A = !!names(.)[2],
    Start_pos_A = !!names(.)[3],
    Chr_name_B = !!names(.)[4],
    Start_pos_B = !!names(.)[5],
    Sequence = !!names(.)[6]
  )

log_step("Filtering by Feature_name_A...") 

df_chr_A_filtered_file_2 <- filter_by_chromosome_A(file_2, chromosomes, df_categories_list) %>% 
  bind_rows(.)

log_step("Filtering by Feature_name_B...") 

df_chr_AB_filtered_file_2 <- filter_by_chromosome_B(df_chr_A_filtered_file_2, chromosomes, df_categories_list) %>% 
  bind_rows(.) %>% 
  mutate(strain = rep(c(strain_name)), sample = rep(c(sample)), experiment = rep(c(experiment)))

df_chr_AB_filtered_file_2_nocontrols <- df_chr_AB_filtered_file_2 %>%
  filter (Category_A != "control", Category_B != "control" | Category_A != "control_norm", Category_B != "control_norm")


remove(df_chr_A_filtered_file_2)

log_step("Finding discordant_pairs not present in T0 sample...")

uncommon_pairs <- anti_join(df_chr_AB_filtered_file_2_nocontrols, df_chr_AB_filtered_file_1_nocontrols, by = c("Feature_name_A", "Feature_name_B")) 

log_step("Processing for BLAST analysis...") 

df_chr_AB_filtered_file_2_ID <- uncommon_pairs %>%
  mutate(pair_id = pmin(paste(Feature_name_A, Feature_name_B, sep = "_"), paste(Feature_name_B, Feature_name_A, sep = "_"))) %>%
  group_by(Read_name, pair_id) %>%
  mutate(pair_group = paste0("id", cur_group_id(), "-", row_number())) %>%
  separate(pair_group, c("pair_group_name", "number"), "-", remove = FALSE)  %>%  
  mutate(Read_name_ID = paste0(Read_name, "_", pair_group)) %>% 
  ungroup() %>%
  select(-pair_id) %>% 
  left_join(., df_categories_near_features, by = "Feature_name_B")

control_file_2 <- df_chr_AB_filtered_file_2 %>% filter (Category_A == "control" & Category_B == "control" | Category_A == "control_norm" & Category_B == "control_norm")

blast_df_file_2 <- prepare_for_blast(df_chr_AB_filtered_file_2_ID) %>% select(Read_name_ID, Sequence, Chromosome_B)


# Write the new TSV
log_step("Saving processed dataframes...") 
log_step(paste0("Saving tsv files for: ", strain_name, " ",  sample, " ", experiment, "..."))

write_tsv(df_chr_AB_filtered_file_2_ID, file = paste0(strain, "/", sample, "_", experiment, "_inter_discordant_pairs_unique_processed.tsv"))
write_tsv(control_file_2, file = paste0(strain, "/", sample, "_", experiment, "_inter_discordant_pairs_unique_processed_control.tsv"))
write_tsv(blast_df_file_2, file = paste0(strain, "/", sample, "_", experiment, "_inter_discordant_pairs_unique_processed_blast_chromosome.tsv"))