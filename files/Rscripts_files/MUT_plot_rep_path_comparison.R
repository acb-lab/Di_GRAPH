#### R script to plot repair pathway comparison
### All strains in MYWD
### 07/03/2026 - Lydia

log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}

# -------------------------
# Read command line args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  log_step("ERROR: Incorrect number of arguments")
  log_step("Usage: script.R <root_dir>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]


log_step(paste("ROOT_DIR:", root_dir))


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



log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}

## Find all. repair_pathway_only_repair_summary.tsv files
# Define a function to process tsv files
process_repair_pathway_only_repair_summary_files <- function(file_path) {
  
  # Read file
  temp_file <- read_tsv(file_path, col_names = TRUE, show_col_types = FALSE)
  
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  # Prepare concordant dataframe
  repair_pathway_only_repair_summary_file <- temp_file 
  
  return(repair_pathway_only_repair_summary_file)
}


log_step("Finding repair_pathway_only_repair_summary files...")
# Get all repair_pathway_only_repair_summary.tsv files recursively in root folder
repair_pathway_only_repair_summary_files <- list.files(
  path = root_dir,
  pattern = "repair_pathway_only_repair_summary\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


repair_pathway_only_repair_summary_all_processed_df <- purrr::map_dfr(repair_pathway_only_repair_summary_files, 
                                                                      process_repair_pathway_only_repair_summary_files) 

samples_order <- c("T0", "TLG", "TLR")
repair_pathway_order <- c("NHEJ", "shHR", "lHR")



repair_pathway_only_repair_summary_all_processed_df_plots <- repair_pathway_only_repair_summary_all_processed_df %>% 
  mutate(repair_pathway = factor(repair_pathway, levels= repair_pathway_order))


bar_plot <- ggplot(repair_pathway_only_repair_summary_all_processed_df_plots, aes(x = strain, y = mean_ratio, fill = repair_pathway)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio), 
                linewidth = 0.8, width = 0.2, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("gold","dodgerblue", "purple3")) +
  theme_classic(base_family = "Helvetica") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                   panel.background = element_blank(), 
                                                   plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  coord_cartesian(ylim = c(0, 1), expand=FALSE) +
  theme(aspect.ratio = 0.7) + 
  scale_x_discrete(name = expression("Strain")) +
  scale_y_continuous(name = expression("Fraction"),
                     #limits = c(0, 1),
                     breaks = seq(0,1,0.2)) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 15, angle = 0), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  labs(
    title = paste0("Repair_pathway - comparison")
  )

ggsave(
  filename = paste0(root_dir, "/", "Repair_pathway_plot_comparison", ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)

log_step("Finding repair_pathway_summary files...")
# Get all repair_pathway_only_repair_summary.tsv files recursively in root folder
repair_pathway_summary_files <- list.files(
  path = root_dir,
  pattern = "repair_pathway_summary\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


repair_pathway_summary_all_processed_df <- purrr::map_dfr(repair_pathway_summary_files, 
                                                          process_repair_pathway_only_repair_summary_files) 




repair_pathway_summary_all_processed_df_plots <- repair_pathway_summary_all_processed_df %>% 
  mutate(repair_pathway = factor(repair_pathway, levels= repair_pathway_order),
         sample = factor(sample, levels = samples_order))

orden_plot <- tibble(position = c(1, 6, 11, 17, 22, 27, 33, 38, 43,
                                  2, 7, 12, 18, 23, 28, 34, 39, 44,
                                  3, 8, 13, 19, 24, 29, 35, 40, 45,
                                  4, 9, 14, 20, 25, 30, 36, 41, 46,
                                  5, 10, 15, 21, 26, 31, 37, 42, 47))

repair_pathway_summary_all_processed_df_plots <- bind_cols(repair_pathway_summary_all_processed_df_plots, orden_plot)


bar_plot <- ggplot(repair_pathway_summary_all_processed_df_plots, aes(x = position, y = mean_ratio, fill = repair_pathway)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio),
                linewidth = 0.8, width = 0.4, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("gold","dodgerblue", "purple3")) +
  theme_classic(base_family = "Helvetica") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                   panel.background = element_blank(), 
                                                   plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  coord_cartesian(ylim = c(0, 1), expand=FALSE) +
  theme(aspect.ratio = 1) + 
  scale_x_discrete(name = expression("Sample")) +
  scale_y_continuous(name = expression("Fraction"),
                     #limits = c(0, 1),
                     breaks = seq(0,1,0.2)) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 15, angle = 0), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  labs(
    title = paste0("Repair_pathway - comparison")
  )

ggsave(
  filename = paste0(root_dir, "/", "Repair_pathway_all_timepoints_plot_comparison", ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)


## Find all mutation_wig_bcf_mut_rate_summary.tsv files
# Define a function to process tsv files
process_mutation_wig_bcf_mut_rate_summary_files <- function(file_path) {
  
  # Read file
  temp_file <- read_tsv(file_path, col_names = TRUE, show_col_types = FALSE)
  
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  # Prepare concordant dataframe
  mutation_wig_bcf_mut_rate_summary_file <- temp_file 
  
  return(mutation_wig_bcf_mut_rate_summary_file)
}


log_step("Finding mutation_wig_bcf_mut_rate_summary files...")
# Get all mutation_wig_bcf_mut_rate_summary.tsv files recursively in root folder
mutation_wig_bcf_mut_rate_summary_files <- list.files(
  path = root_dir,
  pattern = "mutation_wig_bcf_mut_rate_summary\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


mutation_wig_bcf_mut_rate_summary_all_processed_df <- purrr::map_dfr(mutation_wig_bcf_mut_rate_summary_files, 
                                                                     process_mutation_wig_bcf_mut_rate_summary_files) %>% 
  filter(repair_pathway == "all")




bar_plot <- ggplot(mutation_wig_bcf_mut_rate_summary_all_processed_df, aes(x = strain, y = mean_mut_rate_global, fill = mutation)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean_mut_rate_global - sd_mut_rate_global, ymax = mean_mut_rate_global + sd_mut_rate_global), 
                linewidth = 0.8, width = 0.2, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("gold", "dodgerblue", "purple3")) +
  theme_classic(base_family = "Helvetica") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                   panel.background = element_blank(), 
                                                   plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  coord_cartesian(ylim = c(0, 4), expand=FALSE) +
  theme(aspect.ratio = 0.7) + 
  scale_x_discrete(name = expression("Strain")) +
  scale_y_continuous(name = expression("Pecentage"),
                     #limits = c(0, 1),
                     breaks = seq(0,6,2)) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 15, angle = 0), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  labs(
    title = paste0("Mutagenic rate - comparison")
  )

ggsave(
  filename = paste0(root_dir, "/", "Mutagenic_rate_plot_comparison", ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)