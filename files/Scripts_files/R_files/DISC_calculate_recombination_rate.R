#### R script to calculate recombination rate (considers if there is a defined reference strain)
### Average from all experiments
### Loop for each strain/ in MYWD
### 13/04/2026 - Lydia



log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}

# -------------------------
# Read command line args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  log_step("ERROR: Incorrect number of arguments")
  log_step("Usage: script.R <root_dir> <reference_strain>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]

reference_strain   <- args[2]


log_step(paste("REFERENCE STRAIN:", reference_strain))

# load libraries

library(ggplot2)
library(extrafont)
library(svglite)
library(purrr)
library(stringr)
library(readr)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)


# Define a function to process tsv files
process_control_files <- function(file_path) {
  
  # Read file
  temp_file <- read_tsv(file_path, col_names = TRUE, show_col_types = FALSE)
  
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  
  # Prepare control dataframe
  control_count_file <- temp_file %>%
    select(strain, sample, experiment, Category_A, Category_B) %>% 
    filter(Category_A == "control_norm" & Category_B == "control_norm")
  
  ncontrol <- nrow(control_count_file)
  
  control_count_file_processed <- control_count_file %>% 
    mutate(control_count = ncontrol) %>% unique() %>%
    select(strain, sample, experiment, control_count)
  
  return(control_count_file_processed)
}


log_step("Finding control files...")
# Get all processed_control.tsv files recursively in root folder
control_files <- list.files(
  path = root_dir,
  pattern = "processed_control\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)



log_step("Processing control files...")
control_count_processed_df <- purrr::map_dfr(control_files, process_control_files)

write_tsv(control_count_processed_df, file.path(root_dir,"control_count.tsv"))

####### changed to adapt to -s flag #### 04/02/2026

if (reference_strain == "no") {
  control_count_processed_df_ratio <- control_count_processed_df %>%
    mutate(Ratio_vs_reference = 1)
}

if (reference_strain != "no") {
  reference_control_count_processed_df <- control_count_processed_df %>%
    filter(strain == reference_strain) %>%
    rename(reference_control_count = control_count) %>%
    select(-strain)
  
  control_count_processed_df_ratio <- control_count_processed_df %>%
    left_join(reference_control_count_processed_df, by = c("sample", "experiment")) %>%
    mutate(Ratio_vs_reference = control_count / reference_control_count)
}

write_tsv(control_count_processed_df_ratio, file.path(root_dir,"control_count_ratio.tsv"))