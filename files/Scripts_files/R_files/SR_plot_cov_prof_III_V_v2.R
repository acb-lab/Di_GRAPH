#### R script to plot SR coverage profiles from chr III and V, 75nt, average all experiments
### Loop for each strain/chr/suffix file in MYWD
### 19/04/2026 - Lydia



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
  log_step("Usage: script.R <file> <chromosome> <suffix> <strain> <root_dir>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

file_path <- args[1]
chr_name <- args[2]
suffix <-  args[3]
strain <-  args[4]
root_dir <- args[5]



strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))
log_step(paste("CHR:", chr_name))
log_step(paste("PLOT:", suffix))

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

file_path <- "/Users/lydia/GWSt/WD_test_v0.3.2/1_Wt/75nt_CHRV_Coverage.tsv"



title_text <- paste("Normalized 75nt", chr_name, suffix)
output_name <- paste0(root_dir, "/plot_", chr_name, "_75nt_", gsub("\\.tsv$", "", suffix), ".svg")

data <- read_tsv(file_path, col_names = TRUE)


data <- data %>%
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

p1 <- ggplot(
  data_normalized,
  aes(x = Position, y = Timepoint, fill = Normalized_Coverage)
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
    breaks = seq(
    min(data_normalized$Position),
    max(data_normalized$Position),
    by = 100000),
    labels = scales::comma) +
  # labs(
  #   title = title_text,
  #   x = "Genomic Position",
  #   y = "Timepoint"
  # ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_line(color ="Black", linewidth =0,1),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
  )

ggsave(output_name, plot = p1, width = 10, height = 6, dpi = 300, device = "svg")