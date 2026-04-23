#### R script to processand plot SR coverage profiles from chr III and V, 75nt and 18nt, average all experiments
### Loop for each subdir/timepoint file in MYWD
### 20/04/2026 - Lydia



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


# root_dir <- ("/Users/lab2.3/GWS3/WD_test_v0.3.2/1_Wt")
# strain <- "1_Wt"
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

# Define a function to prepare coverage files for plotting (after processing)
prepare_for_plotting <- function(coverage_df) {
  data <- coverage_df %>%
    mutate(
      Timepoint = factor(sample, levels = rev(c("T0", "TSG", "TLG", "TLR"))),
      Coverage = as.numeric(avg_coverage_binned),
      Position = coordinate
    )
  
  t0_values <- data %>%
    filter(Timepoint == "T0") %>%
    select(Position, T0_Coverage = Coverage)
  
  data_normalized <- data %>%
    left_join(t0_values, by = "Position") %>%
    mutate(
      Normalized_Coverage = ifelse(T0_Coverage == 0, NA, Coverage / T0_Coverage)
    )
  
  data_normalized$Position <- as.factor(data_normalized$Position)
  
  return(data_normalized)
}

# Define a function to prepare coverage files for plotting (after processing)
prepare_for_plotting_MATa <- function(coverage_df) {
  data <- coverage_df %>%
    mutate(
      Timepoint = factor(sample, levels = rev(c("T0", "TSG", "TLG", "TLR"))),
      Coverage = as.numeric(avg_coverage),
      Position = coordinate
    )
  
  t0_values <- data %>%
    filter(Timepoint == "T0") %>%
    select(Position, T0_Coverage = Coverage)
  
  data_normalized <- data %>%
    left_join(t0_values, by = "Position") %>%
    mutate(
      Normalized_Coverage = ifelse(T0_Coverage == 0, NA, Coverage / T0_Coverage)
    )
  
  data_normalized$Position <- as.factor(data_normalized$Position)
  
  return(data_normalized)
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

processed_75nt_CHRIII_Coverage <- prepare_for_plotting(CHRIII_avg_coverage_binned)
processed_75nt_CHRV_Coverage <- prepare_for_plotting(CHRV_avg_coverage_binned)

###
#Plotting

log_step("Plotting...")

chromosome_plot <- ggplot(
  processed_75nt_CHRIII_Coverage,
  aes(x = bin, y = Timepoint, fill = Normalized_Coverage)
) +
  geom_raster() +
  scale_fill_gradientn(
    colors = rev(c("#AF2418", "#E1AC40", "#EFD24D", "#5D8B27", "#4EACE9", "#4573A1", "#4C1F8E")),
    values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)),
    na.value = "gray90",
    name = "Normalized Coverage",
    limits = c(0, 2),
    oob = scales::squish
  ) +
  scale_x_continuous(
    breaks = c(0, 1000, 2000, 3000),
    labels = scales::comma(c(1, 100000, 200000, 300000))) +
  labs(
    title = "Normalized 75nt CHRIII Coverage",
    x = "Genomic Position",
    y = "Timepoint"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_line(color ="Black", linewidth =0,1),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
  )


ggsave(
  filename = paste0(root_dir,"/", "plot_CHRIII_75nt_Coverage",".svg"),
  plot = chromosome_plot,
  width = 10,
  height = 6,
  dpi = 300, 
  device = svglite,
  bg = "transparent"
)

chromosome_plot <- ggplot(
  processed_75nt_CHRV_Coverage,
  aes(x = bin, y = Timepoint, fill = Normalized_Coverage)
) +
  geom_raster() +
  scale_fill_gradientn(
    colors = rev(c("#AF2418", "#E1AC40", "#EFD24D", "#5D8B27", "#4EACE9", "#4573A1", "#4C1F8E")),
    values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)),
    na.value = "gray90",
    name = "Normalized Coverage",
    limits = c(0, 2),
    oob = scales::squish
  ) +
  scale_x_continuous(
    breaks = c(0, 2000, 4000, 6000),
    labels = scales::comma(c(1, 200000, 400000, 600000))) +
  labs(
    title = "Normalized 75nt CHRV Coverage",
    x = "Genomic Position",
    y = "Timepoint"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_line(color ="Black", linewidth =0,1),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
  )


ggsave(
  filename = paste0(root_dir,"/", "plot_CHRV_75nt_Coverage",".svg"),
  plot = chromosome_plot,
  width = 10,
  height = 6,
  dpi = 300, 
  device = svglite,
  bg = "transparent"
)

####

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

####

processed_CHRIII_MATa_avg_coverage <- prepare_for_plotting_MATa(CHRIII_MATa_avg_coverage)
processed_CHRV_MATa_avg_coverage <- prepare_for_plotting_MATa(CHRV_MATa_avg_coverage)

###
#Plotting

log_step("Plotting...")

chromosome_plot <- ggplot(
  processed_CHRIII_MATa_avg_coverage,
  aes(x = coordinate, y = Timepoint, fill = Normalized_Coverage)
) +
  geom_raster() +
  scale_fill_gradientn(
    colors = rev(c("#AF2418", "#E1AC40", "#EFD24D", "#5D8B27", "#4EACE9", "#4573A1", "#4C1F8E")),
    values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)),
    na.value = "gray90",
    name = "Normalized Coverage",
    limits = c(0, 2),
    oob = scales::squish
  ) +
  # scale_x_continuous(
  #   breaks = c(0, 1000, 2000, 3000),
  #   labels = scales::comma(c(1, 100000, 200000, 300000))) +
  labs(
    title = "Normalized 75nt CHRIII MATa_Coverage",
    x = "Genomic Position",
    y = "Timepoint"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_line(color ="Black", linewidth =0,1),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
  )


ggsave(
  filename = paste0(root_dir,"/", "plot_CHRIII_75nt_MATa_Coverage",".svg"),
  plot = chromosome_plot,
  width = 10,
  height = 6,
  dpi = 300, 
  device = svglite,
  bg = "transparent"
)

chromosome_plot <- ggplot(
  processed_CHRV_MATa_avg_coverage,
  aes(x = coordinate, y = Timepoint, fill = Normalized_Coverage)
) +
  geom_raster() +
  scale_fill_gradientn(
    colors = rev(c("#AF2418", "#E1AC40", "#EFD24D", "#5D8B27", "#4EACE9", "#4573A1", "#4C1F8E")),
    values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1, 1.25, 1.5, 1.75, 2)),
    na.value = "gray90",
    name = "Normalized Coverage",
    limits = c(0, 2),
    oob = scales::squish
  ) +
  # scale_x_continuous(
  #   breaks = c(0, 1000, 2000, 3000),
  #   labels = scales::comma(c(1, 100000, 200000, 300000))) +
  labs(
    title = "Normalized 75nt CHRV MATa_Coverage",
    x = "Genomic Position",
    y = "Timepoint"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_line(color ="Black", linewidth =0,1),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
  )


ggsave(
  filename = paste0(root_dir,"/", "plot_CHRV_75nt_MATa_Coverage",".svg"),
  plot = chromosome_plot,
  width = 10,
  height = 6,
  dpi = 300, 
  device = svglite,
  bg = "transparent"
)
####

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


####
log_step("Plotting...")
all_data <- bind_rows(CHRIII_18nt_MATa_avg_coverage, CHRV_18nt_MATa_avg_coverage)

for (s in unique(all_data$sample)) {
  df <- subset(all_data, sample == s)
  
  p <- ggplot(df, aes(x = coordinate_corrected, y = avg_coverage, color = chromosome)) +
    geom_line(linewidth = 0.8) +
    labs(
      title = paste("Coverage of ChrIII and ChrV - Timepoint", s),
      x = "Coordinate (relative to MAT locus)",
      y = "Coverage"
    ) +
    scale_color_manual(values = c("CHRIII" = "#2F2C7E", "CHRV" = "#A30000")) +
    scale_x_continuous(breaks = seq(-800, 800, by = 200)) +
    coord_cartesian(ylim = c(0, 3)) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "right",
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA)
    )
  
  
  ggsave(
    filename = paste0(root_dir,"/", "plot_", s,"_18nt_MATa_Coverage.svg"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300, 
    device = svglite,
    bg = "transparent"
  )
 
}
