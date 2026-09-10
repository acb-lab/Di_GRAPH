#### R script to identify valid discordant pairs (5 different levels of BLAST validation) for 10kb validation
### For all experiments, BLAST validation +- 1 genomic feature
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
  log_step("Usage: script.R <strain> <sample> <experiment> <category_path> <file1> <fileorig> <filecontrol>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

strain <- args[1]
sample   <- args[2]
experiment   <- args[3]
category_path   <- args[4]
path_to_file_1   <- args[5]

path_to_file_original   <- args[6]
path_to_file_control   <- args[7]

strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))
log_step(paste("SAMPLE:", sample))
log_step(paste("EXPERIMENT:", experiment))




# load libraries

library(readr)
library(extrafont)
library(stringr)
library(svglite)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

# Variables from Bash

path_to_categories_pairs_file <- file.path(category_path, "PMV_categories_pairs.tsv")
path_to_features_pairs_file <- file.path(category_path, "PMV_features_pairs.tsv")


# Functions

#path_to_file_1
#path_to_file_2
#path_to_file_3
#path_to_file_original

# Functions
log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}

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



get_blast_minimum_distance_complete_df <- function (chromosome_file, original_file) {
  
  complete_df <- chromosome_file %>% rename(Read_name_ID = pair_group) %>% 
    full_join(., original_file, by = "Read_name_ID") %>% 
    mutate(min_distance = pmin(
      abs(start_pos - Start_pos_B),
      abs(end_pos - Start_pos_B),
      na.rm = FALSE)) %>% 
    group_by(Read_name_ID) %>% arrange(min_distance) %>% 
    slice(1) %>%
    ungroup() %>% 
    select(!c("Feature_name_B.x")) 
  
  
  return(complete_df)
  
}     




chromosome_df <- process_blast_filtered_files(path_to_file_1)

original_df <- read_tsv(path_to_file_original,col_names = TRUE) %>% 
  mutate(Read_name_ID = with(., paste0(Read_name, "_", pair_group))) %>%
  filter(Category_A != "control_norm", Category_B != "control_norm", Category_A != "control", Category_B != "control")


blast_processesd_minimum_distance_df <- get_blast_minimum_distance_complete_df(chromosome_df, original_df)

#  Valid: 0hits/+-10kb
process_blast_filtered_files_option1 <- function(blast_processesd_minimum_distance_file) {
  blast_filtered_df <- blast_processesd_minimum_distance_file %>% 
    mutate(valid_10kb = ifelse(Start_pos_B < start_pos_10kb | Start_pos_B > end_pos_10kb | hits == 0, "yes", "no")) %>%
    group_by(pair_group_name) %>%
    mutate(valid_all = if_else(all(valid_10kb == "yes"), "yes", "no")) %>%
    ungroup() %>%
    filter(valid_all == "yes")
  
  return(blast_filtered_df)
}



# Valid: >=1 hit & >= 3 mismatch & lenght = 75nt + 0hits/+-10kb
process_blast_filtered_files_option2 <- function(blast_processesd_minimum_distance_file) {
  blast_filtered_df <- blast_processesd_minimum_distance_file %>% 
    mutate(valid_10kb = ifelse(Start_pos_B < start_pos_10kb | Start_pos_B > end_pos_10kb | hits == 0, "yes", ifelse(hits !=0 & mismatch >2 & length == 75, "yes", "no"))) %>%
    group_by(pair_group_name) %>%
    mutate(valid_all = if_else(all(valid_10kb == "yes"), "yes", "no")) %>%
    ungroup() %>%
    filter(valid_all == "yes")
  
  return(blast_filtered_df)
}



# Valid: >=1 hit & = 2 mismatch & lenght = 75nt + >=1 hit & >= 3 mismatch & lenght = 75nt + 0hits/+-10kb
process_blast_filtered_files_option3 <- function(blast_processesd_minimum_distance_file) {
  blast_filtered_df <- blast_processesd_minimum_distance_file %>% 
    mutate(valid_10kb = ifelse(Start_pos_B < start_pos_10kb | Start_pos_B > end_pos_10kb | hits == 0, "yes", 
                               ifelse(hits !=0 & mismatch >2 & length == 75, "yes", 
                                      ifelse(hits !=0 & mismatch == 2 & length == 75, "yes", "no")))) %>%
    group_by(pair_group_name) %>%
    mutate(valid_all = if_else(all(valid_10kb == "yes"), "yes", "no")) %>%
    ungroup() %>%
    filter(valid_all == "yes")
  
  return(blast_filtered_df)
}

# Valid: >=1 hit & = 1 mismatch & lenght = 75nt + >=1 hit & = 2 mismatch & lenght = 75nt + >=1 hit & >= 3 mismatch & lenght = 75nt + 0hits/+-10kb 
process_blast_filtered_files_option4 <- function(blast_processesd_minimum_distance_file) {
  blast_filtered_df <- blast_processesd_minimum_distance_file %>% 
    mutate(valid_10kb = ifelse(Start_pos_B < start_pos_10kb | Start_pos_B > end_pos_10kb | hits == 0, "yes", 
                               ifelse(hits !=0 & mismatch >2 & length == 75, "yes", 
                                      ifelse(hits !=0 & mismatch == 2 & length == 75, "yes", 
                                             ifelse(hits != 0 & mismatch == 1 & length == 75, "yes", "no"))))) %>%
    group_by(pair_group_name) %>%
    mutate(valid_all = if_else(all(valid_10kb == "yes"), "yes", "no")) %>%
    ungroup() %>%
    filter(valid_all == "yes")
  
  return(blast_filtered_df)
}


# Valid: >=1 hit & = 1 mismatch & lenght = 75nt + >=1 hit & = 2 mismatch & lenght = 75nt + >=1 hit & >= 3 mismatch & lenght = 75nt + 0hits/+-10kb
process_blast_filtered_files_option5 <- function(blast_processesd_minimum_distance_file) {
  blast_filtered_df <- blast_processesd_minimum_distance_file %>% 
    mutate(valid_10kb = "yes") %>%
    group_by(pair_group_name) %>%
    mutate(valid_all = if_else(all(valid_10kb == "yes"), "yes", "no")) %>%
    ungroup() %>%
    filter(valid_all == "yes")
  
  return(blast_filtered_df)
}
###


chromosome_df_option1 <- process_blast_filtered_files_option1(blast_processesd_minimum_distance_df) 
chromosome_df_option2 <- process_blast_filtered_files_option2(blast_processesd_minimum_distance_df) 
chromosome_df_option3 <- process_blast_filtered_files_option3(blast_processesd_minimum_distance_df) 
chromosome_df_option4 <- process_blast_filtered_files_option4(blast_processesd_minimum_distance_df) 
chromosome_df_option5 <- process_blast_filtered_files_option5(blast_processesd_minimum_distance_df) 


log_step("Finding valid reads...") 

valid_read_name_option1 <- chromosome_df_option1 %>% select(Read_name_ID)
valid_read_name_option2 <- chromosome_df_option2 %>% select(Read_name_ID)
valid_read_name_option3 <- chromosome_df_option3 %>% select(Read_name_ID)
valid_read_name_option4 <- chromosome_df_option4 %>% select(Read_name_ID)
valid_read_name_option5 <- chromosome_df_option5 %>% select(Read_name_ID)




log_step("Comparing with original file...") 

valid_pairs_df_option1 <- semi_join(original_df, valid_read_name_option1)
valid_pairs_df_option2 <- semi_join(original_df, valid_read_name_option2)
valid_pairs_df_option3 <- semi_join(original_df, valid_read_name_option3)
valid_pairs_df_option4 <- semi_join(original_df, valid_read_name_option4)
valid_pairs_df_option5 <- semi_join(original_df, valid_read_name_option5)

n_valid_reads_option1 <- nrow(valid_pairs_df_option1)
n_valid_reads_option2 <- nrow(valid_pairs_df_option2)
n_valid_reads_option3 <- nrow(valid_pairs_df_option3)
n_valid_reads_option4 <- nrow(valid_pairs_df_option4)
n_valid_reads_option5 <- nrow(valid_pairs_df_option5)

n_original_reads <- nrow(original_df)

log_step("Calculating error rate...")

error_rate_option1 <- tibble(strain = strain_name, 
                             sample = sample,
                             experiment = experiment,
                             original_reads = n_original_reads,
                             valid_reads = n_valid_reads_option1) %>% 
  mutate(valid_rate = (valid_reads/original_reads)*100)

error_rate_option2 <- tibble(strain = strain_name, 
                             sample = sample,
                             experiment = experiment,
                             original_reads = n_original_reads,
                             valid_reads = n_valid_reads_option2) %>% 
  mutate(valid_rate = (valid_reads/original_reads)*100)

error_rate_option3 <- tibble(strain = strain_name, 
                             sample = sample,
                             experiment = experiment,
                             original_reads = n_original_reads,
                             valid_reads = n_valid_reads_option3) %>% 
  mutate(valid_rate = (valid_reads/original_reads)*100)

error_rate_option4 <- tibble(strain = strain_name, 
                             sample = sample,
                             experiment = experiment,
                             original_reads = n_original_reads,
                             valid_reads = n_valid_reads_option4) %>% 
  mutate(valid_rate = (valid_reads/original_reads)*100)

error_rate_option5 <- tibble(strain = strain_name, 
                             sample = sample,
                             experiment = experiment,
                             original_reads = n_original_reads,
                             valid_reads = n_valid_reads_option5) %>% 
  mutate(valid_rate = (valid_reads/original_reads)*100)




# Write the new TSV
log_step("Saving processed dataframe...") 
log_step(paste0("Saving tsv files for: ", strain_name, " ",  sample, " ", experiment, "..."))

write_tsv(valid_pairs_df_option1, file = paste0(strain, "/", sample, "_", experiment, "_option1","_inter_discordant_pairs_unique_processed_valid.tsv"))
write_tsv(valid_pairs_df_option2, file = paste0(strain, "/", sample, "_", experiment, "_option2","_inter_discordant_pairs_unique_processed_valid.tsv"))
write_tsv(valid_pairs_df_option3, file = paste0(strain, "/", sample, "_", experiment, "_option3","_inter_discordant_pairs_unique_processed_valid.tsv"))
write_tsv(valid_pairs_df_option4, file = paste0(strain, "/", sample, "_", experiment, "_option4","_inter_discordant_pairs_unique_processed_valid.tsv"))
write_tsv(valid_pairs_df_option5, file = paste0(strain, "/", sample, "_", experiment, "_option5","_inter_discordant_pairs_unique_processed_valid.tsv"))


write_tsv(error_rate_option1, file = paste0(strain, "/", sample, "_", experiment, "_option1", "_inter_discordant_pairs_unique_processed_valid_error_rate.tsv"))
write_tsv(error_rate_option2, file = paste0(strain, "/", sample, "_", experiment, "_option2", "_inter_discordant_pairs_unique_processed_valid_error_rate.tsv"))
write_tsv(error_rate_option3, file = paste0(strain, "/", sample, "_", experiment, "_option3", "_inter_discordant_pairs_unique_processed_valid_error_rate.tsv"))
write_tsv(error_rate_option4, file = paste0(strain, "/", sample, "_", experiment, "_option4", "_inter_discordant_pairs_unique_processed_valid_error_rate.tsv"))
write_tsv(error_rate_option5, file = paste0(strain, "/", sample, "_", experiment, "_option5", "_inter_discordant_pairs_unique_processed_valid_error_rate.tsv"))