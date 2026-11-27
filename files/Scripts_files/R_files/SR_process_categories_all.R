#### R script to process and plot SR coverage profiles for genomic categories, 75nt, average all experiments
### Single-session version of SR_process_categories.R (used by Di-GRAPH_piped.sh)
###
### SR_process_categories.R was started once per category (13 R sessions per strain), and
### every session re-read the same 12 genome-wide *_75nt.tsv per-base coverage files
### (~12M rows each). This script reads them ONCE and then runs the unchanged
### per-category computation, TSV output and plot for every category given.
### Results are identical to running SR_process_categories.R once per category.
###
### Usage: SR_process_categories_all.R <root_dir> <strain> <category_path> <status_file> <category>...
### <status_file> receives one line per category: "<category>\tOK" or "<category>\tERROR",
### so the calling shell can write the same log lines as with one session per category.

log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}

# -------------------------
# Read command line args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  log_step("ERROR: Incorrect number of arguments")
  log_step("Usage: script.R <root_dir> <strain> <category_path> <status_file> <category>...")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]
category_path <- args[3]
status_file <- args[4]
categories <- args[5:length(args)]

strain <- sub("/$", "", strain)
strain_name <- basename(strain)

log_step(paste("STRAIN:", strain))
log_step(paste("CATEGORIES:", paste(categories, collapse = " ")))

write_status <- function(status) {
  writeLines(paste(categories, status[categories], sep = "\t"), status_file)
}
status <- setNames(rep("ERROR", length(categories)), categories)
write_status(status)

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


# Chromosome positions
chr_positions <- list(
  CHRI = 1:230218, CHRII = 1:813184, CHRIII = 1:316513,
  CHRIV = 1:1531933, CHRV = 1:578179, CHRVI = 1:270161,
  CHRVII = 1:1090940, CHRVIII = 1:562643, CHRIX = 1:439888,
  CHRX = 1:745751, CHRXI = 1:666816, CHRXII = 1:1078177,
  CHRXIII = 1:924431, CHRXIV = 1:784333, CHRXV = 1:1091291,
  CHRXVI = 1:948066
)


# -------------------------------------------------------------------------
# Part 1 (once): read every coverage file and split it by chromosome.
# Same steps as the first half of process_coverage_files() in the original.
# Only the columns used later are kept (chromosome, coverage, coordinate);
# strain/sample/experiment are constant per file and are added back after
# filtering, before the same final select(), so the rows and values are the same.
# -------------------------------------------------------------------------
read_coverage_file <- function(file_path) {
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  # directory name above the file

  parts <- str_split(file_base, "_", simplify = TRUE)

  if (ncol(parts) >= 3) {
    sample_name <- parts[1]
    experiment_name <- parts[2]
  } else {
    warning(paste("Filename does not match expected format:", file_base))
    return(NULL)
  }

  temp_file <- read_tsv(file_path, col_names = FALSE, show_col_types = FALSE)

  if (nrow(temp_file) == 0) {
    return(NULL)
  }

  coverage_file <- temp_file %>%
    rename("chromosome" = !!names(.[1]),
           "coverage" = !!names(.[4])) %>%
    select(chromosome, coverage)
  rm(temp_file)

  chr_dfs <- lapply(names(chr_positions), function(chr) {
    coverage_file %>% filter(chromosome == chr) %>%
      mutate(coordinate = chr_positions[[chr]])
  })
  names(chr_dfs) <- names(chr_positions)

  list(strain = strain_name, sample = sample_name, experiment = experiment_name, chr_dfs = chr_dfs)
}


# -------------------------------------------------------------------------
# Part 2 (per category): second half of process_coverage_files(), unchanged.
# -------------------------------------------------------------------------
process_coverage_file <- function(cov, categories_df) {
  if (is.null(cov)) {
    return(NULL)
  }

  processed_df <- map_dfr(names(chr_positions), function(chr) {
    chr_df <- cov$chr_dfs[[chr]]

    cat_df <- categories_df %>% filter(Cromosoma == chr)

    map_dfr(unique(cat_df$Categoria), function(cat) {
      regions <- cat_df %>% filter(Categoria == cat)
      map_dfr(1:nrow(regions), function(i) {
        chr_df %>%
          filter(coordinate >= regions$Pos_inicio[i], coordinate <= regions$Pos_fin[i]) %>%
          mutate(Feature_name_A = cat, Category_A = regions$Tipo[i])
      })
    })
  })

  processed_df <- processed_df %>%
    mutate(strain = cov$strain, sample = cov$sample, experiment = cov$experiment) %>%
    select(strain, sample, experiment, chromosome, coverage, Feature_name_A, Category_A)
  return(processed_df)
}


process_category <- function(category, coverage_data) {
  log_step(paste("CATEGORY:", category))

  log_step("Finding category file...")
  category_file <- list.files(
    path = category_path,
    pattern = paste0(category, ".*\\.tsv$"),
    recursive = TRUE,
    full.names = TRUE
  )

  # Read category regions
  df_categories <- read_tsv(category_file, col_names = FALSE,
                            col_types = cols(X1 = col_character(), X2 = col_character(),
                                             X3 = col_double(), X4 = col_double(), X5 = col_character())) %>%
    rename(Tipo = X1, Categoria = X2, Pos_inicio = X3, Pos_fin = X4, Cromosoma = X5) %>%
    mutate(Categoria = as.factor(Categoria))

  log_step("Processing _75nt coverage files...")
  coverage_processed_df <- purrr::map_dfr(coverage_data, process_coverage_file, categories_df = df_categories)

  log_step("Processing _75nt coverage files summary...")
  coverage_processed_summary_df <- coverage_processed_df %>%
    group_by(strain, sample, experiment, chromosome, Feature_name_A, Category_A) %>%
    summarise(avg_coverage_feature = mean(coverage),
              sd_coverage_feature = sd(coverage)) %>%
    ungroup() %>%
    mutate(ratio_coverage_feature = avg_coverage_feature / avg_coverage_feature[sample == "T0"]) %>%
    ungroup() %>%
    group_by(strain, sample, chromosome, Feature_name_A, Category_A) %>%
    summarise(avg_ratio_coverage_feature = mean(ratio_coverage_feature),
              sd_ratio_coverage_feature = sd(ratio_coverage_feature)) %>%
    ungroup() %>%
    filter(sample != "T0")

  ## Save tsv file for each sample
  for (s in unique(coverage_processed_summary_df$sample)) {
    df <- subset(coverage_processed_summary_df, sample == s)
    write_tsv(df, paste0(root_dir,"/", category, "_fingerprint_", strain_name,"_T0vs", s,".tsv"))
  }

  log_step("Plotting...")
  ## Save svg plot for each sample
  for (s in unique(coverage_processed_summary_df$sample)) {
    df <- subset(coverage_processed_summary_df, sample == s)

    p <- ggplot(df) +
      geom_col(aes(y = avg_ratio_coverage_feature, x = Feature_name_A), fill = "#2F2C7E", width = 0.7,
               color = "black", linewidth = 0.1) +
      geom_errorbar(aes(x = Feature_name_A, ymin = avg_ratio_coverage_feature - sd_ratio_coverage_feature,
                        ymax = avg_ratio_coverage_feature + sd_ratio_coverage_feature),
                    width = 0.3, color = "grey8", alpha = 1, size = 0.3) +
      theme_classic() +
      theme(
        panel.grid = element_blank(),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA)
      ) +
      labs(title = paste0("Ratio T0vs", s), subtitle = paste0("Strain ", strain_name),
           y = "Ratio", x = "Feature")

    ggsave(
      filename = paste0(root_dir,"/", strain_name, "_", category, "_T0vs", s,"_plot.svg"),
      plot = p,
      width = 8,
      height = 6,
      dpi = 300,
      device = svglite,
      bg = "transparent"
    )
  }
  invisible(TRUE)
}


log_step("Finding _75nt coverage files...")
# Get all coverage.tsv files recursively in root folder (once: the category outputs
# written below do not match this pattern, so every category saw the same list)
coverage_files <- list.files(
  path = root_dir,
  pattern = "_75nt\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Reading _75nt coverage files (once for all categories)...")
coverage_data <- tryCatch(
  lapply(coverage_files, read_coverage_file),
  error = function(e) {
    # The original would have failed in every category on the same file
    log_step(paste("ERROR reading coverage files:", conditionMessage(e)))
    write_status(status)
    quit(status = 1)
  }
)

# Each category runs independently, as when each had its own R session
for (category in categories) {
  ok <- tryCatch(process_category(category, coverage_data),
                 error = function(e) {
                   log_step(paste("ERROR in category", category, ":", conditionMessage(e)))
                   FALSE
                 })
  status[category] <- if (isTRUE(ok)) "OK" else "ERROR"
  write_status(status)
}

quit(status = if (all(status == "OK")) 0 else 1)
