#### R script to processand plot SR coverage profiles from chr III and V, 75nt and 18nt, average all experiments
### Loop for each subdir/timepoint file in MYWD
### 21/04/2026 - Lydia



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



# Define a function to process tsv files
process_diff_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename((dirname(dirname(file_path))))  # directory name above the file
  
  # 
  parts <- str_split(file_base, "_", simplify = TRUE)
  
  # Validate and extract parts safely
  if (ncol(parts) >= 5) {
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
  
  # Prepare dataframe
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




log_step("Finding CHRIII diff files...")
# Get all coverage.tsv files recursively in root folder
CHRIII_coverage_files <- list.files(
  path = root_dir,
  pattern = "CHRIII_diff.*\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing CHRIII diff files...")
CHRIII_avg_coverage <-  purrr::map_dfr(CHRIII_coverage_files, process_diff_files) %>% 
  group_by(strain, sample, chromosome, coordinate) %>% 
  summarise(avg_coverage = mean(coverage),
            sd_coverage = sd(coverage)) %>% 
  ungroup()


log_step("Finding CHRV coverage files...")
# Get all coverage.tsv files recursively in root folder
CHRV_coverage_files <- list.files(
  path = root_dir,
  pattern = "CHRV_diff.*\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


log_step("Processing CHRV coverage files...")
CHRV_avg_coverage <-  purrr::map_dfr(CHRV_coverage_files, process_diff_files) %>% 
  group_by(strain, sample, chromosome, coordinate) %>% 
  summarise(avg_coverage = mean(coverage),
            sd_coverage = sd(coverage)) %>% 
  ungroup()


coords_III <- tibble(
  coordinate = c(200119, 200167, 200212, 200272,
                 200326, 200386, 200449, 200509, 
                 200542, 200575, 200635, 200689,
                 200753,
                 200817, 200882, 200947, 201012,
                 201077, 201148, 201207, 201272, 
                 201337, 201402),
  rel_coordinate = c(-634, -586, -541, -481,
                 -427, -367, -304, -244,
                 -211, -178, -118, -64,
                 0,
                 64, 129, 194, 259,
                 324, 395, 454, 519,
                 584, 649))

coords_V <- tibble(
  coordinate = c(289191, 289239, 289284, 289344,
                 289398, 289458, 289521, 289581,
                 289614, 289647, 289707, 289761,
                 289825,
                 289889, 289954, 290019, 290084, 
                 290149, 290220, 290278, 290343,
                 290408, 290473),
  rel_coordinate = c(-634, -586, -541, -481,
                     -427, -367, -304, -244,
                     -211, -178, -118, -64,
                     0,
                     64, 129, 194, 259,
                     324, 395, 454, 519,
                     584, 649))


CHRIII_avg_coverage <- left_join(CHRIII_avg_coverage, coords_III, by = "coordinate")
CHRV_avg_coverage <- left_join(CHRV_avg_coverage, coords_V, by = "coordinate")

data <- bind_rows(CHRIII_avg_coverage, CHRV_avg_coverage)


###

for (s in unique(data$sample)) {
  df <- subset(data, sample == s)
  
  p <- ggplot(df, aes(x=rel_coordinate, y=avg_coverage, color=chromosome, fill=chromosome)) +
    geom_ribbon(aes(ymin=avg_coverage - sd_coverage, ymax=avg_coverage + sd_coverage),
                alpha=0.2, color=NA) +
    geom_line(size=1.2) +
    scale_x_continuous(breaks=unique(df$coordinate)) +  # Only show x-ticks for coordinates in column 4
    coord_cartesian(ylim = c(-125, 125)) +
    # scale_y_continuous(limits = c(-125, 125)) +
    scale_color_manual(values = c("CHRIII" = "#2F2C7E", "CHRV" = "#A30000")) +
    scale_fill_manual(values = c("CHRIII" = "#2F2C7E", "CHRV" = "#A30000")) +
    labs(
      title=paste0("Quant_Coverage of ChrIII and ChrV polymorphisms - ", s),
      x="Relative Coordinate (bp)",
      y="Average Coverage") +
    theme(
      axis.text.x = element_text(angle=90, vjust=0.5, hjust=1),
      panel.grid.minor.x = element_blank(),  # Remove vertical minor grid lines
      panel.grid.major.y = element_line(color = "gray", size = 0.5),  # Keep horizontal grid lines
      panel.grid = element_line(color ="Black", linewidth =0,1),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA),
    )
  

  
  ggsave(
    filename = paste0(root_dir,"/", "plot_", s,"_18nt_MATa_Coverage_Quant.svg"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300,
    device = svglite,
    bg = "transparent"
  )
  
  
}
###
