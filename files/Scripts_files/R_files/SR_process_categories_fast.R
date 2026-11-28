#### R script to process and plot SR coverage profiles for genomic categories, 75nt, average all experiments
### Fast single-session version of SR_process_categories.R (used by Di-GRAPH_fast.sh)
###
### Same inputs, outputs and arguments as SR_process_categories_all.R:
###   SR_process_categories_fast.R <root_dir> <strain> <category_path> <status_file> <category>...
###
### What changed (results are identical):
###  1. The 12 genome-wide *_75nt.tsv files are read once for all categories
###     (as in SR_process_categories_all.R), keeping only the coverage vector of
###     each chromosome.
###  2. The original filtered the whole chromosome table once per feature
###     (filter(coordinate >= start, coordinate <= end), ~6,600 dplyr calls per
###     chromosome set for ORF/Int) and bound all selected bases into one table
###     (~140M rows for ORF) only to average them per feature. Because the original
###     sets coordinate = row number (1..chromosome length), the bases of a feature
###     are simply rows ceiling(start)..floor(end). This version takes that slice and
###     calls the same mean() and sd() on the same values, in the same order, so every
###     number is bit-identical. The summary table is then built with the same
###     group_by() keys and types, so its rows come out in the same order.
###  3. Safety net: if a category file has the same feature twice on one chromosome
###     (never the case in PMV v9), that category runs the original per-feature code.

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
# Read every coverage file once; keep one coverage vector per chromosome.
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

  chromosome <- temp_file[[1]]
  coverage <- temp_file[[4]]
  rm(temp_file)

  chr_cov <- lapply(names(chr_positions), function(chr) {
    x <- coverage[which(chromosome == chr)]          # = filter(chromosome == chr)
    if (length(x) != length(chr_positions[[chr]])) {
      # the original mutate(coordinate = chr_positions[[chr]]) fails here too
      stop(sprintf("%s: %s has %d rows, expected %d", file_base, chr, length(x), length(chr_positions[[chr]])))
    }
    x
  })
  names(chr_cov) <- names(chr_positions)

  list(strain = strain_name, sample = sample_name, experiment = experiment_name, chr_cov = chr_cov)
}


# -------------------------------------------------------------------------
# Fast path: per-feature mean/sd straight from the coverage slice.
# Returns one row per (file, feature with >= 1 base), with the same key columns
# and types the original group_by() saw.
# -------------------------------------------------------------------------
feature_stats <- function(cov, categories_df) {
  if (is.null(cov)) return(NULL)
  lvls <- levels(categories_df$Categoria)
  out <- lapply(names(chr_positions), function(chr) {
    cat_df <- categories_df %>% filter(Cromosoma == chr)
    if (nrow(cat_df) == 0) return(NULL)
    x <- cov$chr_cov[[chr]]
    n <- length(x)
    # rows of the original filter: coordinate (= 1..n) within [Pos_inicio, Pos_fin]
    s <- pmax(1, ceiling(cat_df$Pos_inicio))
    e <- pmin(n, floor(cat_df$Pos_fin))
    keep <- !is.na(s) & !is.na(e) & s <= e
    if (!any(keep)) return(NULL)
    s <- s[keep]; e <- e[keep]
    avg <- vapply(seq_along(s), function(k) mean(x[s[k]:e[k]]), numeric(1))
    sdv <- vapply(seq_along(s), function(k) sd(x[s[k]:e[k]]), numeric(1))
    tibble(
      strain = cov$strain, sample = cov$sample, experiment = cov$experiment,
      chromosome = chr,
      Feature_name_A = factor(as.character(cat_df$Categoria[keep]), levels = lvls),
      Category_A = cat_df$Tipo[keep],
      avg = avg, sdv = sdv
    )
  })
  bind_rows(out)
}


# -------------------------------------------------------------------------
# Original per-feature code (safety net for category files with duplicated
# features on a chromosome). Same as process_coverage_files() in the original.
# -------------------------------------------------------------------------
process_coverage_file_slow <- function(cov, categories_df) {
  if (is.null(cov)) return(NULL)
  processed_df <- map_dfr(names(chr_positions), function(chr) {
    chr_df <- tibble(chromosome = chr, coverage = cov$chr_cov[[chr]]) %>%
      mutate(coordinate = chr_positions[[chr]])
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
  processed_df %>%
    mutate(strain = cov$strain, sample = cov$sample, experiment = cov$experiment) %>%
    select(strain, sample, experiment, chromosome, coverage, Feature_name_A, Category_A)
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

  fast_ok <- !any(duplicated(df_categories[, c("Cromosoma", "Categoria")]))

  if (fast_ok) {
    log_step("Per-feature coverage statistics (fast path)...")
    stats <- purrr::map_dfr(coverage_data, feature_stats, categories_df = df_categories)
    # Each group holds exactly one feature of one file, so taking its value
    # equals mean()/sd() over the group's bases, which were computed above.
    coverage_summary_1 <- stats %>%
      group_by(strain, sample, experiment, chromosome, Feature_name_A, Category_A) %>%
      summarise(avg_coverage_feature = avg[1],
                sd_coverage_feature = sdv[1],
                n_rows = n())
    stopifnot(all(coverage_summary_1$n_rows == 1))
    coverage_summary_1 <- coverage_summary_1 %>% select(-n_rows)
  } else {
    log_step("Duplicated features in category file: using the original per-feature code...")
    coverage_processed_df <- purrr::map_dfr(coverage_data, process_coverage_file_slow, categories_df = df_categories)
    coverage_summary_1 <- coverage_processed_df %>%
      group_by(strain, sample, experiment, chromosome, Feature_name_A, Category_A) %>%
      summarise(avg_coverage_feature = mean(coverage),
                sd_coverage_feature = sd(coverage))
  }

  log_step("Processing _75nt coverage files summary...")
  # From here on: unchanged from the original
  coverage_processed_summary_df <- coverage_summary_1 %>%
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
    log_step(paste("ERROR reading coverage files:", conditionMessage(e)))
    write_status(status)
    quit(status = 1)
  }
)

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
