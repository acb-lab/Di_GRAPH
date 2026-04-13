#### R script to calculate read number of discordant reads in each Category_A and plot bar plot with SD
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
library(ggbreak)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)


# Define a function to process tsv files
process_valid_counts_files <- function(file_path) {
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
# Get all coverage.tsv files recursively in root folder
valid_counts_files <- list.files(
  path = root_dir,
  pattern = "valid_counts\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
valid_counts_files
log_step("Processing valid_counts files...")
valid_counts_processed_df <- purrr::map_dfr(valid_counts_files, process_valid_counts_files)



log_step("Finding ratio control file...")

ratio_control_file <- read_tsv(file.path(wd_dir, "control_count_ratio.tsv"), col_names = TRUE)



complete_valid_counts_df <- left_join(valid_counts_processed_df, ratio_control_file) %>%  
  mutate(total_count_norm = total_count / Ratio_vs_reference)

complete_valid_counts_summary_df <- complete_valid_counts_df %>% 
  group_by(strain, sample, Category_A) %>% 
  summarise(mean_total_count = mean(total_count),
            sd_total_count = sd(total_count),
            mean_total_count_norm = mean(total_count_norm),
            sd_total_count_norm = sd(total_count_norm))


genomic_categories_order <- c(
  "ORF", "intergenic", "long_terminal_repeat", "transposable_element_gene",
  "LTR_retrotransposon", "tRNA_gene", "rRNA_gene", "ncRNA_gene",
  "snRNA_gene", "snoRNA_gene", "ARS", "centromere", "telomere"
)

samples_order <- c("TSG", "TLG", "TLR")



complete_valid_counts_summary_df_ordered <- complete_valid_counts_summary_df %>% mutate(Category_A = factor(Category_A, levels = genomic_categories_order),
                                                                                        sample = factor(sample, levels = samples_order)) %>%
  arrange(strain, sample, Category_A)

write_tsv(complete_valid_counts_summary_df_ordered, file.path(root_dir, paste0(strain_name, "_category_A_valid_counts_summary_df.tsv")))



log_step("Plotting...")
# Generate plots

bar_plot <- ggplot(complete_valid_counts_summary_df_ordered, aes(x = Category_A, y = mean_total_count_norm, fill = sample)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean_total_count_norm - sd_total_count_norm, ymax = mean_total_count_norm + sd_total_count_norm),
                linewidth = 0.8, width = 0.5, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#252159", "#469CD7", "#8BE0FC")) +
  theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                               panel.background = element_blank(), 
                                               plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  scale_x_discrete(name = expression(""),
                   labels = c("ORF", "Intergenic","LTR","TEG", "Ty", "tRNA", "rRNA", "ncRNA", "snRNA", "snoRNA", "ARS", "Centromere", "Telomere")) +
  scale_y_break(c(125, 1000), scales = 0.5, ticklabels = c(seq(0, 125, 25), 1000, 3000, 5000, 7000)) +
  scale_y_continuous(
    limits = c(-10, 10000),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 20, angle = 90), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_blank()
  ) +
  labs(
    title = paste0("Category_A norm.read number  - ", strain_name)
  )

ggsave(
  filename = paste0(strain,"/", "Category_A_read_number_norm_plot_", strain_name, ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)

log_step("Plotting...")
# Generate plots

bar_plot <- ggplot(complete_valid_counts_summary_df_ordered, aes(x = Category_A, y = mean_total_count, fill = sample)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean_total_count - sd_total_count, ymax = mean_total_count + sd_total_count),
                linewidth = 0.8, width = 0.5, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("#252159", "#469CD7", "#8BE0FC")) +
  theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                               panel.background = element_blank(), 
                                               plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  scale_x_discrete(name = expression(""),
                   labels = c("ORF", "Intergenic","LTR","TEG", "Ty", "tRNA", "rRNA", "ncRNA", "snRNA", "snoRNA", "ARS", "Centromere", "Telomere")) +
  scale_y_break(c(125, 1000), scales = 0.5, ticklabels = c(seq(0, 125, 25), 1000, 3000, 5000, 7000)) +
  scale_y_continuous(
    limits = c(-10, 10000),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 20, angle = 90), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  theme(
    axis.line.y.right = element_blank(),
    axis.ticks.y.right = element_blank(),
    axis.text.y.right = element_blank(),
    axis.title.y.right = element_blank()
  ) +
  labs(
    title = paste0("Category_A read number  - ", strain_name)
  )

ggsave(
  filename = paste0(strain,"/", "Category_A_read_number_no_norm_plot_", strain_name, ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)