#### R script to plot SR coverage profiles for CHRIII and CHRV
### Loop for each strain/chromosome/suffix file in MYWD
### 08/03/2026 - Lydia


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
  log_step("Usage: script.R <file_path> <chr> <suffix> <subdir>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

file_path <- args[1]
chr_name <- args[2]
suffix <- args[3]
subdir <- args[4]

log_step(paste("FILE:", file_path))
log_step(paste("CHR:", chr_name))
log_step(paste("SUFFIX:", suffix))
log_step(paste("SUBDIR:", subdir))

library(ggplot2)
library(extrafont)
library(dplyr)

title_text <- paste("Normalized 75nt", chr_name, suffix)
output_name <- paste0(subdir, "/plot_", chr_name, "_75nt_", gsub("\\.tsv$", "", suffix), ".svg")

data <- read.table(file_path, header = FALSE)
colnames(data) <- c("Position", "Timepoint", "Coverage")

data <- data %>%
  mutate(
    Timepoint = factor(Timepoint, levels = rev(c("0", "SG", "LG", "LR"))),
    Coverage = as.numeric(Coverage)
  )

t0_values <- data %>%
  filter(Timepoint == "0") %>%
  select(Position, T0_Coverage = Coverage)

data_normalized <- data %>%
  left_join(t0_values, by = "Position") %>%
  mutate(
    Normalized_Coverage = ifelse(T0_Coverage == 0, NA, Coverage / T0_Coverage)
  )

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
  scale_x_continuous(labels = scales::comma) +
  labs(
    title = title_text,
    x = "Genomic Position",
    y = "Timepoint"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_line(color ="Black", linewidth =0.1),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

ggsave(output_name, plot = p1, width = 10, height = 6, dpi = 300, device = "svg")