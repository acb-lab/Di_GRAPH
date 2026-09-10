#### R script to calculate error rate
### Average from all experiments
### Loop for each strain/ in MYWD
### 13/04/2026 - Lydia

log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}

# -------------------------
# Read command line args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  log_step("ERROR: Incorrect number of arguments")
  log_step("Usage: script.R <root_dir> <strain> <wd_dir>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]
wd_dir <- args[3]


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
process_error_rate_files <- function(file_path) {
  
  # Read file
  temp_file <- read_tsv(file_path, col_names = TRUE, show_col_types = FALSE)
  
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  
  # Prepare error_rate dataframe
  error_rate_file <- temp_file %>%
    select(strain, sample, experiment, valid_rate) 
  
  error_rate_file_processed <- error_rate_file %>% 
    mutate(error_rate = 100-valid_rate) 
  
  
  return(error_rate_file_processed)
}


log_step("Finding error_rate files...")
# Get all error_rate.tsv files recursively in root folder
error_rate_files <- list.files(
  path = root_dir,
  pattern = "processed_valid_error_rate\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)



log_step("Processing error_rate files...")
error_rate_processed_df <- purrr::map_dfr(error_rate_files, process_error_rate_files)

error_rate_processed_summary_df <- error_rate_processed_df %>%  group_by(strain, sample) %>% 
  summarise(mean_error_rate = mean(error_rate, na.rm = TRUE), 
            sd_error_rate = sd(error_rate, na.rm = TRUE)) %>% 
  filter(sample != "T0") %>% 
  ungroup

write_tsv(error_rate_processed_summary_df, file.path(root_dir, paste0(strain_name, "_error_rate_summary.tsv")))


samples_order <- c("TSG", "TLG", "TLR")

log_step("Plotting...")
# Generate plots

error_rate_processed_summary_df <- error_rate_processed_summary_df %>% 
  mutate(sample = factor(sample, levels = samples_order))





log_step("Plotting...")
# Generate plots

bar_plot <- ggplot(error_rate_processed_summary_df, aes(x = sample, y = mean_error_rate, fill = sample)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean_error_rate - sd_error_rate, ymax = mean_error_rate + sd_error_rate), 
                linewidth = 0.8, width = 0.2, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#252159", "#469CD7", "#8BE0FC")) +
  theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                               panel.background = element_blank(), 
                                               plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  coord_cartesian(ylim = c(0, 50), expand=FALSE) +
  theme(aspect.ratio = 0.75) + 
  scale_x_discrete(name = expression("Sample")) +
  scale_y_continuous(name = expression("Percentage"),
                     #limits = c(0, 1),
                     breaks = seq(0,50,10)) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 20, angle = 0), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  labs(
    title = paste0("Error rate - ", strain_name)
  )

ggsave(
  filename = paste0(strain,"/", "Error_rate_plot_", strain_name, ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)