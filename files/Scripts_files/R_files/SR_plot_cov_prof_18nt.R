#### R script to plot SR coverage profiles for T0 and TLG, reads 18nt
### Loop for each subdir/timepoint file in MYWD
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
  log_step("Usage: script.R <file_chrIII> <file_chrV> <timepoint> <subdir>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

file_chrIII <- args[1]
file_chrV <- args[2]
timepoint <- args[3]
subdir <- args[4]

log_step(paste("file_chrIII:", file_chrIII))
log_step(paste("file_chrV:", file_chrV))
log_step(paste("timepoint:", timepoint))
log_step(paste("SUBDIR:", subdir))



# Load necessary libraries
library(ggplot2)
library(extrafont)
library(dplyr)

# # Assign file paths from the bash variables
# file_chrIII <- "$file_chrIII"
# file_chrV <- "$file_chrV"
# timepoint <- "$timepoint"  # Pass the timepoint variable to R

# Read the dataset for chrIII and chrV
data_chrIII <- read.table(file_chrIII, header=FALSE, col.names=c("Coordinate", "Timepoint", "Coverage", "Chromosome"))
data_chrV <- read.table(file_chrV, header=FALSE, col.names=c("Coordinate", "Timepoint", "Coverage", "Chromosome"))

# Combine the two datasets into one
all_data <- rbind(data_chrIII, data_chrV)

# Check if the necessary columns exist
if (!all(c("Coordinate", "Timepoint", "Coverage", "Chromosome") %in% colnames(all_data))) {
  stop("Error: Required columns are missing in the combined data")
}

# Plot the data using ggplot
p2 <- ggplot(all_data, aes(x=Coordinate, y=Coverage, color=Chromosome)) +
  geom_line(linewidth=0.8) +
  labs(
    title=paste("Coverage of ChrIII and ChrV - Timepoint", timepoint),
    x="Coordinate (relative to MAT locus)",
    y="Coverage"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust=0.5),  # Center the title
    legend.position="right",  # Position the legend
    panel.grid = element_line(color ="Black", linewidth =0,1),
    panel.background = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
  ) +
  scale_color_manual(values=c("CHRIII"="#2F2C7E", "CHRV"="#A30000")) +  # Color for each chromosome
  coord_cartesian(ylim=c(0, 3))  # Optional: Set y-axis limits
scale_x_continuous(breaks=seq(-800, 800, by=200))  # Set x-axis breaks to be at intervals of 200

# Define the output file path
output_name <- file.path("$subdir", paste0("plot_", timepoint, "_18nt_MATa_Coverage.svg"))

# Save the plot
ggsave(output_name, plot=p2, width=10, height=6, dpi=300, device="svg")