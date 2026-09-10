#### R script to plot MATa-MATa' coverage (reads 75nt)
### Loop for each strain in MYWD
### 11/04/2026 - Lydia

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
strain   <- args[2]


strain <- sub("/$", "", strain)
strain_name <- basename(strain)

log_step(paste("ROOT_DIR:", root_dir))
log_step(paste("STRAIN:", strain))

# load libraries

library(ggplot2)
library(extrafont)
library(svglite)
library(ggdensity)
library(purrr)
library(stringr)
library(readr)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)



# Define a function to process tsv files
process_MATs_coverage_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  # directory name above the file
  
  # Expect filename like T4_E07_MAT_filtered.tsv
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
  
  # Prepare coverage dataframe
  coverage_file <- temp_file %>%
    rename("chromosome" = !!names(.[1]), "coverage" = !!names(.[4])) %>%
    select(chromosome, coverage) %>%
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name
    )
  
  return(coverage_file)
}


log_step("Finding MATa-MATa' coverage files...")
# Get all coverage.tsv files recursively in root folder
MAT_coverage_files <- list.files(
  path = root_dir,
  pattern = "MAT_filtered\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
MAT_coverage_files
log_step("Processing MATa-MATa' coverage files...")
MATs_processed_df <- purrr::map_dfr(MAT_coverage_files, process_MATs_coverage_files)


log_step("Processing chromosome III...")
III_df <- MATs_processed_df %>% filter(chromosome == "CHRIII") %>%
  group_by(strain, sample, experiment) %>% 
  mutate(position = row_number(), position_HO_centered = position - 200753) %>%
  ungroup() %>% 
  group_by(strain, sample, chromosome, position, position_HO_centered) %>%
  summarise(mean_coverage = mean(coverage, na.rm = TRUE), sd_coverage = sd(coverage, na.rm = TRUE))

log_step("Processing chromosome V...")
V_df <- MATs_processed_df %>% filter(chromosome == "CHRV") %>%
  group_by(strain, sample, experiment) %>% 
  mutate(position = row_number(), position_HO_centered = position - 289825) %>%
  ungroup() %>% 
  group_by(strain, sample, chromosome, position, position_HO_centered) %>%
  summarise(mean_coverage = mean(coverage, na.rm = TRUE), sd_coverage = sd(coverage, na.rm = TRUE))

MAT_df <- bind_rows(III_df, V_df)
log_step(paste("Number of samples found:", length(MAT_df)))


# Set polymorphisms positions
polymorphisms <- tibble(position_HO_centered = c(-634,-586,-541,-481,-427,-367,-304,-244,-211,-178,-118,-64,
                                                 0,64,129,194,259,324,395,454,519,584,649),
                        mean_coverage = rep(c(0)),
                        chromosome = rep(c("polymorphism")))

log_step("Splitting by sample...")
# Split by sample
split_MAT_df <- split(MAT_df, 
                      MAT_df$sample)

# Set output directory
#svg_output_dir <- file.path(root_dir, "MAT_plots")
#dir.create(svg_output_dir, showWarnings = FALSE, recursive = TRUE)
write_tsv(MAT_df, file.path(root_dir, paste0(strain_name, "_MAT_coverage_df.tsv")))


log_step("Plotting...")
# Generate plots
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  
  p <- ggplot(data_subset, aes(x = position_HO_centered, y = mean_coverage, color = chromosome)) +
    geom_point(data = polymorphisms, aes(x = position_HO_centered, y = mean_coverage),
               color = "chartreuse3", shape = 108, size = 5) +
    geom_line(linewidth = 1) + #1.2
    scale_color_manual(values = c("CHRIII" = "#2F2C7E", "CHRV" = "#A30000")) + #2F2C7E, A30000
    theme_classic(base_family = "Arial") +
    theme(
      panel.grid = element_line(color = "black", linewidth = 0.1),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none",
      aspect.ratio = 0.45,
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(hjust = 1)
    ) +
    theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 20),
          axis.title.y = element_text(vjust = 2, size = 20)) +
    theme(axis.text.x = element_text(vjust = 0, size = 15),
          axis.text.y = element_text(vjust = 0, size = 15)) +
    coord_cartesian(xlim = c(-1200, 1200), ylim = c(0, 30), expand = FALSE) +
    scale_x_continuous(limits = c(-1100, 1100), breaks = seq(-1100, 1100, 200)) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey3", linewidth = 0.5) +
    labs(
      title = paste0("Average coverage of MATa-MATa' discordant reads - ", strain_name, " - ", sample),
      x = "Relative Coordinate (bp)",
      y = "Average Coverage"
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("Average_MATa_MATa_coverage_", strain_name, "_", sample, ".svg")),
    plot = p,
    #width = 8,
    #height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})

# Generate plots for chromosome III (only coverage shape)
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  data_subset <- data_subset %>% filter(chromosome == "CHRIII")
  
  p <- ggplot(data_subset, aes(x = position_HO_centered, y = mean_coverage, color = chromosome)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("CHRIII" = "#2F2C7E")) +
    theme_classic(base_family = "Arial") +
    theme(
      panel.grid = element_line(color = "black", linewidth = 0.1),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none",
      aspect.ratio = 0.25,
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(hjust = 1)
    ) +
    theme(axis.line=element_blank(), axis.ticks=element_blank()) +
    theme(axis.text.x = element_text(vjust = 0, size = 0),
          axis.text.y = element_text(vjust = 0, size = 0))+
    coord_cartesian(xlim = c(-1200, 1200), ylim = c(0, 30), expand = FALSE) +
    scale_x_continuous(limits = c(-1100, 1100), breaks = seq(-1100, 1100, 200)) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
    labs(
      title = paste0("Average coverage of MATa-MATa' discordant reads chrIII - ", strain_name, " - ", sample),
      x = "",
      y = ""
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("Average_MATa_MATa_coverage_chrIII_", strain_name, "_", sample, ".svg")),
    plot = p,
    width = 8,
    height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})

# Generate plots for chromosome V (only coverage shape)
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  data_subset <- data_subset %>% filter(chromosome == "CHRV")
  
  p <- ggplot(data_subset, aes(x = position_HO_centered, y = mean_coverage, color = chromosome)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("CHRV" = "#A30000")) +
    theme_classic(base_family = "Arial") +
    theme(
      panel.grid = element_line(color = "black", linewidth = 0.1),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none",
      aspect.ratio = 0.25,
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(hjust = 1)
    ) +
    theme(axis.line=element_blank(), axis.ticks=element_blank()) +
    theme(axis.text.x = element_text(vjust = 0, size = 0),
          axis.text.y = element_text(vjust = 0, size = 0))+
    coord_cartesian(xlim = c(1200, -1200), ylim = c(0, 30), expand = FALSE) +
    scale_x_continuous(limits = c(-1100, 1100), breaks = seq(-1100, 1100, 200)) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
    labs(
      title = paste0("Average coverage of MATa-MATa' discordant reads chrV - ", strain_name, " - ", sample),
      x = "",
      y = ""
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("Average_MATa_MATa_coverage_chrV_", strain_name, "_", sample, ".svg")),
    plot = p,
    width = 8,
    height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})

# Define a function to process sam files
process_MATs_pairs_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  # directory name above the file
  
  
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
  temp_file <- read_tsv(file_path, skip = 19, col_names = FALSE, show_col_types = FALSE)
  
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  # Prepare coverage dataframe
  pairs_file <- temp_file %>%
    rename("chromosome_A" = !!names(.[3]), "position_A" = !!names(.[4]),
           "chromosome_B" = !!names(.[7]), "position_B" = !!names(.[8])) %>%
    select(chromosome_A, position_A, chromosome_B, position_B) %>%
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name
    )
  
  return(pairs_file)
}

log_step("Finding MATa-MATa' pairs files...")
# Get all pairs files recursively in root folder
MAT_pairs_files <- list.files(
  path = root_dir,
  pattern = "unique_MAT\\.sam$",
  recursive = TRUE,
  full.names = TRUE
)
MAT_pairs_files
log_step("Processing MATa-MATa' pairs files...")
MATs_pairs_processed_df <- purrr::map_dfr(MAT_pairs_files, process_MATs_pairs_files)

MAT_pairs_III <- MATs_pairs_processed_df %>% filter(chromosome_A == "CHRIII" & chromosome_B == "CHRV") %>% 
  rename(X_value_III = position_A, Y_value_V = position_B) %>% 
  select(X_value_III, Y_value_V, strain, sample, experiment)

MAT_pairs_V <- MATs_pairs_processed_df %>% filter(chromosome_A == "CHRV" & chromosome_B == "CHRIII") %>% 
  rename(X_value_III = position_B, Y_value_V = position_A) %>% 
  select(X_value_III, Y_value_V, strain, sample, experiment)

MAT_combined <- bind_rows(MAT_pairs_III, MAT_pairs_V)
write_tsv(MAT_combined, file.path(root_dir, paste0(strain_name, "_MAT_pairs_df.tsv")))

### Calculate proportion of GC events to the right or to the left of HO site
MAT_combined_bar_graph <- MAT_combined %>% mutate(right_left = ifelse(Y_value_V < 289825, "left", "right")) %>% 
  group_by(strain, sample, experiment, right_left) %>% 
  summarise(n_experiment = n()) %>%
  mutate(total_experiment = sum(n_experiment)) %>%                  
  mutate(percentage_experiment = 100 * n_experiment / total_experiment) %>%         # Percent for each right/left
  ungroup() %>% select(!c(n_experiment, total_experiment)) 
MAT_combined_summary <- MAT_combined_bar_graph %>% 
  group_by(strain, sample, right_left) %>% 
  summarise(average_percentage = mean(percentage_experiment),
            sd_percentage = sd(percentage_experiment))

write_tsv(MAT_combined_summary, file.path(root_dir, paste0(strain_name, "_MAT_pairs_summary_df.tsv")))


split_MAT_df <- split(MAT_combined, 
                      MAT_combined$sample)


log_step("Plotting...")
# Generate plots for MATa-MATa' pairs
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  
  p <- ggplot(data_subset, aes(x = X_value_III, y = Y_value_V)) + theme_light(base_family = "Arial") +
    theme(panel.background = element_blank(), plot.background = element_rect(fill = "transparent", colour = NA)) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    scale_x_continuous(breaks = seq(199653,201853, by = 1100), limits = c(199353,202153)) + 
    scale_y_continuous(breaks = seq(288725,290925, by = 1100), limits = c(288425,291225))+ 
    geom_hdr(xlim = c(199353,202153), ylim = c(288425,291225), method = "kde", fill = "brown4") + 
    geom_point (colour = "black", size = 0.5) + 
    theme(aspect.ratio = 1) + 
    geom_hline(yintercept=289825, linetype="dashed", color = "black", linewidth=0.2, alpha = 0.3) + 
    geom_vline(xintercept=200753, linetype="dashed", color = "black", linewidth=0.2, alpha = 0.3) +
    labs(
      title = paste0("MATa-MATa' discordant reads pairs - ", strain_name, " - ",sample),
      x = "Chr.III coordinates",
      y = "Chr.V coordinates"
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("MATa_MATa_pairs_", strain_name, "_", sample, ".svg")),
    plot = p,
    width = 8,
    height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})