#### R script to calculate % of discordant reads in each Category_A and plot bar plot with SD
### Average from all experiments
### Loop for each strain/ in MYWD
### 12/04/2026 - Lydia

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
process_valid_counts_files <- function(file_path) {
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
  temp_file <- read_tsv(file_path, col_names = TRUE, show_col_types = FALSE)
  
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  # Prepare counts distribution dataframe
  counts_distribution_file <- temp_file %>%
    group_by(strain, sample, experiment, Category_A) %>%
    summarise(total_count = sum(count), .groups = "drop_last") %>%
    mutate(
      group_total = sum(total_count),
      percent = 100 * (total_count / group_total)
    ) %>%
    ungroup()
  
  return(counts_distribution_file)
}


log_step("Finding valid_counts files...")
# Get all valid counts.tsv files recursively in root folder
valid_counts_files <- list.files(
  path = root_dir,
  pattern = "valid_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
valid_counts_files
log_step("Processing valid_counts files...")
valid_counts_processed_df <- purrr::map_dfr(valid_counts_files, process_valid_counts_files)

genomic_categories_order <- c(
  "ORF", "intergenic", "long_terminal_repeat", "transposable_element_gene",
  "LTR_retrotransposon", "tRNA_gene", "rRNA_gene", "ncRNA_gene",
  "snRNA_gene", "snoRNA_gene", "ARS", "centromere", "telomere"
)

samples_order <- c("TSG", "TLG", "TLR")

valid_counts_processed_df_summary <- valid_counts_processed_df %>%
  group_by(strain, sample, Category_A) %>%
  summarise(
    mean_percent = mean(percent, na.rm = TRUE),
    sd_percent = sd(percent, na.rm = TRUE),
    .groups = "drop"
  )

valid_counts_processsed_df_summary_ordered <- valid_counts_processed_df_summary %>% mutate(Category_A = factor(Category_A, levels = genomic_categories_order)) %>%
  arrange(strain, sample, Category_A)

#write_tsv(valid_counts_processed_df_summary, file.path(root_dir, paste0(strain_name, "_valid_counts_summary.tsv")))
write_tsv(valid_counts_processsed_df_summary_ordered, file.path(root_dir, paste0(strain_name, "_valid_counts_summary_ordered.tsv")))
#write_tsv(valid_counts_processed_df, file.path(root_dir, paste0(strain_name, "_valid_counts.tsv")))

log_step("Plotting...")
# Generate plots

valid_counts_processsed_df_summary_ordered <-  valid_counts_processsed_df_summary_ordered %>% 
  mutate(Category_A = factor(Category_A, levels = genomic_categories_order)) %>% 
  mutate(sample = factor(sample, levels = samples_order))

bar_plot <- ggplot(valid_counts_processsed_df_summary_ordered, aes(x = Category_A, y = mean_percent, fill = sample)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean_percent - sd_percent, ymax = mean_percent + sd_percent), 
                linewidth = 0.8, width = 0.5, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#252159", "#469CD7", "#8BE0FC")) +
  theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                               panel.background = element_blank(), 
                                               plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  #coord_cartesian(ylim = c(0, 6), expand=FALSE) +
  coord_cartesian(ylim = c(0, 100), expand=FALSE) +
  theme(aspect.ratio = 0.75) + 
  scale_x_discrete(name = expression("Category"), 
                   labels = c("ORF", "Intergenic","LTR","TEG", "Ty", "tRNA", "rRNA", "ncRNA", "snRNA", "snoRNA", "ARS", "Centromere", "Telomere")) +
  scale_y_continuous(name = expression("Percentage"),
                     limits = c(0, 100),
                     breaks = seq(0,100,10)) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 20, angle = 90), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  labs(
    title = paste0("Discordant reads distribution - ", strain_name)
  )

ggsave(
  filename = paste0(strain,"/", "Discordant_reads_distribution_plot_", strain_name, ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)

bar_plot_reduced <- ggplot(valid_counts_processsed_df_summary_ordered, aes(x = Category_A, y = mean_percent, fill = sample)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean_percent - sd_percent, ymax = mean_percent + sd_percent), 
                linewidth = 0.8, width = 0.5, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#252159", "#469CD7", "#8BE0FC")) +
  theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                               panel.background = element_blank(), 
                                               plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  coord_cartesian(ylim = c(0, 6), expand=FALSE) +
  #coord_cartesian(ylim = c(0, 100), expand=FALSE) +
  theme(aspect.ratio = 0.75) + 
  scale_x_discrete(name = expression("Category"), 
                   labels = c("ORF", "Intergenic","LTR","TEG", "Ty", "tRNA", "rRNA", "ncRNA", "snRNA", "snoRNA", "ARS", "Centromere", "Telomere")) +
  scale_y_continuous(name = expression("Percentage"),
                     limits = c(0, 100),
                     breaks = seq(0,100,2)) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 20, angle = 90), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  labs(
    title = paste0("Discordant reads distribution - reduced - ", strain_name)
  )

ggsave(
  filename = paste0(strain,"/", "Discordant_reads_distribution_plot_reduced_", strain_name, ".svg"),
  plot = bar_plot_reduced,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)
