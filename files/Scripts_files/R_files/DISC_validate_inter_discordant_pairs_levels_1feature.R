#### R script to identify valid discordant pairs (5 different levels of BLAST validation)
### For all experiments, BLAST validation +- 1 genomic feature
### Loop for each strain/sample/experiment in MYWD
### 14/04/2026 - Lydia


log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}

# -------------------------
# Read command line args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 9) {
  log_step("ERROR: Incorrect number of arguments")
  log_step("Usage: script.R <strain> <sample> <experiment> <category_path> <file1> <file2> <file3> <fileorig> <filecontrol>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

strain <- args[1]
sample   <- args[2]
experiment   <- args[3]
category_path   <- args[4]
path_to_file_1   <- args[5]
path_to_file_2   <- args[6]
path_to_file_3   <- args[7]
path_to_file_original   <- args[8]
path_to_file_control   <- args[9]

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

# 0 hits: valid
process_blast_filtered_files_option1 <- function (blast_filtered_file_path) {
  blast_filtered_df <- read.table(blast_filtered_file_path,
                                  header=FALSE, sep = "\t", 
                                  col.names = paste0("V",seq_len(5)), fill = TRUE) %>% 
    rename("pair_group" = !!names(.[1]), "Feature_name_B" = !!names(.[2]), "hits" = !!names(.[3]), 
           "mismatch" = !!names(.[4]), "length" = !!names(.[5])) %>% 
    mutate(valid = ifelse(hits == 0, "yes", "no"))
  
  return(blast_filtered_df)
  
}

# Valid: >=1 hit & >= 3 mismatch & lenght = 75nt + 0hits
process_blast_filtered_files_option2 <- function (blast_filtered_file_path) {
  blast_filtered_df <- read.table(blast_filtered_file_path,
                                  header=FALSE, sep = "\t", 
                                  col.names = paste0("V",seq_len(5)), fill = TRUE) %>% 
    rename("pair_group" = !!names(.[1]), "Feature_name_B" = !!names(.[2]), "hits" = !!names(.[3]), 
           "mismatch" = !!names(.[4]), "length" = !!names(.[5])) %>% 
    mutate(valid = ifelse(hits == 0, "yes", ifelse(hits !=0 & mismatch >2 & length == 75, "yes", "no")))
  
  return(blast_filtered_df)
  
}

# Valid: >=1 hit & = 2 mismatch & lenght = 75nt + >=1 hit & >= 3 mismatch & lenght = 75nt + 0hits 
process_blast_filtered_files_option3 <- function (blast_filtered_file_path) {
  blast_filtered_df <- read.table(blast_filtered_file_path,
                                  header=FALSE, sep = "\t", 
                                  col.names = paste0("V",seq_len(5)), fill = TRUE) %>% 
    rename("pair_group" = !!names(.[1]), "Feature_name_B" = !!names(.[2]), "hits" = !!names(.[3]), 
           "mismatch" = !!names(.[4]), "length" = !!names(.[5])) %>% 
    mutate(valid = ifelse(hits == 0, "yes", 
                          ifelse(hits !=0 & mismatch >2 & length == 75, "yes", 
                                 ifelse(hits !=0 & mismatch == 2 & length == 75, "yes", "no"))))
  
  return(blast_filtered_df)
  
}

# Valid: >=1 hit & = 1 mismatch & lenght = 75nt + >=1 hit & = 2 mismatch & lenght = 75nt + >=1 hit & >= 3 mismatch & lenght = 75nt + 0hits 
process_blast_filtered_files_option4 <- function (blast_filtered_file_path) {
  blast_filtered_df <- read.table(blast_filtered_file_path,
                                  header=FALSE, sep = "\t", 
                                  col.names = paste0("V",seq_len(5)), fill = TRUE) %>% 
    rename("pair_group" = !!names(.[1]), "Feature_name_B" = !!names(.[2]), "hits" = !!names(.[3]), 
           "mismatch" = !!names(.[4]), "length" = !!names(.[5])) %>% 
    mutate(valid = ifelse(hits == 0, "yes", 
                          ifelse(hits !=0 & mismatch >2 & length == 75, "yes", 
                                 ifelse(hits !=0 & mismatch == 2 & length == 75, "yes", 
                                        ifelse(hits != 0 & mismatch == 1 & length == 75, "yes", "no")))))
  
  return(blast_filtered_df)
  
}

# Valid: all + >=1 hit & = 1 mismatch & lenght = 75nt + >=1 hit & = 2 mismatch & lenght = 75nt + >=1 hit & >= 3 mismatch & lenght = 75nt + 0hits 
process_blast_filtered_files_option5 <- function (blast_filtered_file_path) {
  blast_filtered_df <- read.table(blast_filtered_file_path,
                                  header=FALSE, sep = "\t", 
                                  col.names = paste0("V",seq_len(5)), fill = TRUE) %>% 
    rename("pair_group" = !!names(.[1]), "Feature_name_B" = !!names(.[2]), "hits" = !!names(.[3]), 
           "mismatch" = !!names(.[4]), "length" = !!names(.[5])) %>% 
    mutate(valid = "yes")
  
  return(blast_filtered_df)
  
}


get_blast_validated_complete_df <- function (center_file, prev_file, next_file) {
  
  complete_df <- left_join(prev_file, next_file, by = "pair_group") %>% 
    left_join(center_file, ., by = "pair_group") %>% 
    rename("valid_prev" = !!names(.[11]), "valid_next" = !!names(.[16])) %>% 
    separate(pair_group, c("pair_group_name", "number"), "-", remove = FALSE) %>% 
    mutate_at(vars(valid_prev, valid_next), ~replace_na(., "yes")) %>% 
    mutate(valid_all = ifelse(valid == "yes" & valid_prev == "yes" & valid_next == "yes", "yes", "no")) %>% 
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


center_df_option1 <- process_blast_filtered_files_option1(path_to_file_1)
center_df_option2 <- process_blast_filtered_files_option2(path_to_file_1)
center_df_option3 <- process_blast_filtered_files_option3(path_to_file_1)
center_df_option4 <- process_blast_filtered_files_option4(path_to_file_1)
center_df_option5 <- process_blast_filtered_files_option5(path_to_file_1)

next_df_option1 <- process_blast_filtered_files_option1(path_to_file_1) ##
next_df_option2 <- process_blast_filtered_files_option2(path_to_file_1)
next_df_option3 <- process_blast_filtered_files_option3(path_to_file_1)
next_df_option4 <- process_blast_filtered_files_option4(path_to_file_1)
next_df_option5 <- process_blast_filtered_files_option5(path_to_file_1)

prev_df_option1 <- process_blast_filtered_files_option1(path_to_file_1) ##
prev_df_option2 <- process_blast_filtered_files_option2(path_to_file_1)
prev_df_option3 <- process_blast_filtered_files_option3(path_to_file_1)
prev_df_option4 <- process_blast_filtered_files_option4(path_to_file_1)
prev_df_option5 <- process_blast_filtered_files_option5(path_to_file_1)




log_step("Finding valid reads...") 

valid_read_name_option1 <- get_blast_validated_complete_df(center_df_option1, prev_df_option1, next_df_option1)
valid_read_name_option2 <- get_blast_validated_complete_df(center_df_option2, prev_df_option2, next_df_option2)
valid_read_name_option3 <- get_blast_validated_complete_df(center_df_option3, prev_df_option3, next_df_option3)
valid_read_name_option4 <- get_blast_validated_complete_df(center_df_option4, prev_df_option4, next_df_option4)
valid_read_name_option5 <- get_blast_validated_complete_df(center_df_option5, prev_df_option5, next_df_option5)
#head(valid_read_name)


original_df <- read_tsv(path_to_file_original,col_names = TRUE) %>% 
  mutate(Read_name_ID = with(., paste0(Read_name, "_", pair_group))) %>%
  filter(Category_A != "control_norm", Category_B != "control_norm", Category_A != "control", Category_B != "control")

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