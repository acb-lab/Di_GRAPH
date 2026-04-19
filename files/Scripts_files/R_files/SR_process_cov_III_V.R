#### R script to process SR coverage profiles from chr III and V, 75nt and 18nt, average all experiments
### Loop for each subdir/timepoint file in MYWD
### 19/04/2026 - Lydia



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
  log_step("Usage: script.R <root_dir> <strain>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]

strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))

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


# root_dir <- ("/Users/lydia/GWSt/WD_test_v0.3.2/1_Wt")
# sample <- "TLR"

# Define a function to process tsv files
process_coverage_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  # directory name above the file
  
  # 
  parts <- str_split(file_base, "_", simplify = TRUE)
  
  # Validate and extract parts safely
  if (ncol(parts) >= 3) {
    sample_name <- parts[1]
    experiment_name <- parts[2]
  } else {
    warning(paste("Filename does not match expected format:", file_base))
    return(NULL)
  }
  
  # Read file
  temp_file <- read_tsv(file_path, col_names = FALSE, show_col_types = FALSE)
  
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  # Prepare nconcordant dataframe
  coverage_file <- temp_file %>%
    rename("chromosome" = !!names(.[1]),
           "coverage" = !!names(.[2]),
           "coordinate" = !!names(.[3])) %>%
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name
    )
  
  return(coverage_file)
}



log_step("Finding CHRIII coverage files...")
# Get all coverage.tsv files recursively in root folder
CHRIII_coverage_files <- list.files(
  path = root_dir,
  pattern = "75nt_CHRIII\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing CHRIII coverage files...")
CHRIII_avg_coverage <-  purrr::map_dfr(CHRIII_coverage_files, process_coverage_files) %>% 
  group_by(strain, sample, chromosome, coordinate) %>% 
  summarise(avg_coverage = mean(coverage),
            sd_coverage = sd(coverage)) %>% 
  ungroup()


log_step("Finding CHRV coverage files...")
# Get all coverage.tsv files recursively in root folder
CHRV_coverage_files <- list.files(
  path = root_dir,
  pattern = "75nt_CHRV\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


log_step("Processing CHRV coverage files...")
CHRV_avg_coverage <-  purrr::map_dfr(CHRV_coverage_files, process_coverage_files) %>% 
  group_by(strain, sample, chromosome, coordinate) %>% 
  summarise(avg_coverage = mean(coverage),
            sd_coverage = sd(coverage)) %>% 
  ungroup()



####
#Binning 100pb
bin_size <- 100

CHRIII_avg_coverage <- CHRIII_avg_coverage %>%
  group_by(strain, sample, chromosome) %>%
  mutate(bin = ceiling(row_number() / bin_size)) %>%
  ungroup()

CHRIII_avg_coverage_binned <- CHRIII_avg_coverage %>%
  group_by(strain, sample, chromosome, bin) %>%
  summarise(
    coordinate = last(coordinate),        
    sample = last(sample),          
    avg_coverage_binned = mean(avg_coverage),        
    .groups = "drop"
  )

write_tsv(CHRIII_avg_coverage_binned, file.path(root_dir, paste0("75nt_CHRIII_Coverage.tsv")))


CHRV_avg_coverage <- CHRV_avg_coverage %>%
  group_by(strain, sample, chromosome) %>%
  mutate(bin = ceiling(row_number() / bin_size)) %>%
  ungroup()

CHRV_avg_coverage_binned <- CHRV_avg_coverage %>%
  group_by(strain, sample, chromosome, bin) %>%
  summarise(
    coordinate = last(coordinate),        
    sample = last(sample),          
    avg_coverage_binned = mean(avg_coverage),        
    .groups = "drop"
  )

write_tsv(CHRV_avg_coverage_binned, file.path(root_dir, paste0("75nt_CHRV_Coverage.tsv")))

######

log_step("Finding CHRIII MATa coverage files...")
# Get all coverage.tsv files recursively in root folder
CHRIII_MATa_coverage_files <- list.files(
  path = root_dir,
  pattern = "75nt_CHRIII_MATa\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


log_step("Processing CHRIII MATa coverage files...")
CHRIII_MATa_avg_coverage <-  purrr::map_dfr(CHRIII_MATa_coverage_files, process_coverage_files) %>% 
  group_by(strain, sample, chromosome, coordinate) %>% 
  summarise(avg_coverage = mean(coverage),
            sd_coverage = sd(coverage)) %>% 
  ungroup()

write_tsv(CHRIII_MATa_avg_coverage, file.path(root_dir, paste0("75nt_CHRIII_MATa_Coverage.tsv")))


log_step("Finding CHRV MATa coverage files...")
# Get all coverage.tsv files recursively in root folder
CHRV_MATa_coverage_files <- list.files(
  path = root_dir,
  pattern = "75nt_CHRV_MATa\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


log_step("Processing CHRV MATa coverage files...")
CHRV_MATa_avg_coverage <-  purrr::map_dfr(CHRV_MATa_coverage_files, process_coverage_files) %>% 
  group_by(strain, sample, chromosome, coordinate) %>% 
  summarise(avg_coverage = mean(coverage),
            sd_coverage = sd(coverage)) %>% 
  ungroup()

write_tsv(CHRV_MATa_avg_coverage, file.path(root_dir, paste0("75nt_CHRV_MATa_Coverage.tsv")))

log_step("Finding CHRIII 18nt MATa coverage files...")
# Get all coverage.tsv files recursively in root folder
CHRIII_18nt_MATa_coverage_files <- list.files(
  path = root_dir,
  pattern = "CHRIII_18nt_ordered\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing CHRIII 18nt MATa coverage files...")
CHRIII_18nt_MATa_avg_coverage <-  purrr::map_dfr(CHRIII_18nt_MATa_coverage_files, process_coverage_files) %>% 
  group_by(strain, sample, chromosome, coordinate) %>% 
  summarise(avg_coverage = mean(coverage),
            sd_coverage = sd(coverage)) %>% 
  ungroup() %>% 
  mutate(coordinate_corrected = coordinate - 200753)

write_tsv(CHRIII_18nt_MATa_avg_coverage, file.path(root_dir, paste0("CHRIII_MATa_18nt_Coverage.tsv")))

log_step("Finding CHRV 18nt MATa coverage files...")
# Get all coverage.tsv files recursively in root folder
CHRV_18nt_MATa_coverage_files <- list.files(
  path = root_dir,
  pattern = "CHRV_18nt_ordered\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing CHRV 18nt MATa coverage files...")
CHRV_18nt_MATa_avg_coverage <-  purrr::map_dfr(CHRV_18nt_MATa_coverage_files, process_coverage_files) %>% 
  group_by(strain, sample, chromosome, coordinate) %>% 
  summarise(avg_coverage = mean(coverage),
            sd_coverage = sd(coverage)) %>% 
  ungroup() %>% 
  mutate(coordinate_corrected = coordinate - 289825)

write_tsv(CHRV_18nt_MATa_avg_coverage, file.path(root_dir, paste0("CHRV_MATa_18nt_Coverage.tsv")))
