#### R script to order datasets from genomic categories, 75nt, average all experiments
### Loop for each subdir/category file in MYWD
### 23/04/2026 - Lydia



log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}

# -------------------------
# Read command line args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  log_step("ERROR: Incorrect number of arguments")
  log_step("Usage: script.R <root_dir> <strain> <category_path> <analysis_suffix>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]

category_path <- args[3]
analysis_suffix <- args[4]


strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))
log_step(paste("ANALYSIS:", analysis_suffix))



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


# Define all categories
categories <- c("ORF", "LTR", "TEG", "Ty", "tRNA", "rRNA", "ncRNA", "snRNA", "snoRNA", "ARS", "Cen", "Tel", "Int")
combined <- list()

# Read the order file
orden_path <- file.path(category_path, "Nombre_ordenado_MATanoalpha_INTERGENICA.tsv")
orden <- read_tsv(orden_path, show_col_types = FALSE) %>% 
  rename(Feature_name_A = Nombre)

log_step("Processing categories files...")
# Loop through all categories
for (category in categories) {
  file_path <- file.path(root_dir, paste0(category, "_fingerprint_", basename(root_dir), "_", analysis_suffix, ".tsv"))
  
  if (file.exists(file_path)) {
    message("Reading: ", file_path)
    df <- read_tsv(file_path, show_col_types = FALSE)
    df_analysis <- data.frame(Analysis = rep(analysis_suffix, nrow(df)))
    completo <- bind_cols(df, df_analysis)
    completo_ordenado <- merge(completo, orden, by = "Feature_name_A")
    combined[[category]] <- completo_ordenado
  }
}


log_step("Combining files...")
# Combine and write final output
if (length(combined) > 0) {
  all_combined <- bind_rows(combined)
  
  # Sort by the first letter and number in 'Nombre_ordenado'
  df_sorted <- all_combined %>%
    mutate(
      Initial = substr(Nombre_ordenado, 1, 1),
      Number = as.numeric(gsub("_.*", "", substr(Nombre_ordenado, 2, nchar(Nombre_ordenado))))
    ) %>%
    arrange(Initial, Number)
  
  # df_sorted
  
  # Save the sorted file
  sorted_output <- file.path(root_dir, paste0("Genomic_sorted_", basename(root_dir), "_", gsub(" ", "", analysis_suffix), ".tsv"))
  write_tsv(df_sorted, sorted_output)
  message("Sorted file written to: ", sorted_output)
} else {
  warning("No data found for analysis: ", analysis_suffix)
}