#### R script to calculate recombination rate (considers if there is a defined reference strain)
### Considering discordant and concordant reads
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
wd_dir   <- args[3]

strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))
log_step(paste("WD_DIR:", wd_dir))

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
process_concordant_row_count_files <- function(file_path) {
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
  concordant_count_file <- temp_file %>%
    rename("concordant_count" = !!names(.[1])) %>%
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name
    )
  
  return(concordant_count_file)
}

# Define a function to process tsv files
process_discordant_valid_files <- function(file_path) {
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
  
  # Prepare ndisc dataframe
  ndisc = nrow(temp_file)
  discordant_count_file <- tibble(discordant_count = ndisc) %>% 
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name
    )
  
  return(discordant_count_file)
}

log_step("Finding concordant row count files...")
# Get all concordant_row_count.tsv files recursively in root folder
concordant_row_count_files <- list.files(
  path = root_dir,
  pattern = "unique_row_count\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Finding discordant row count files...")

discordant_row_count_files <- list.files(
  path = root_dir,
  pattern = "unique_processed_valid\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Finding ratio control file...")

ratio_control_file <- read_tsv(file.path(wd_dir, "control_count_ratio.tsv"), col_names = TRUE)



log_step("Processing concordant row count files...")
concordant_processed_df <- purrr::map_dfr(concordant_row_count_files, process_concordant_row_count_files)

log_step("Processing discordant row count files...")
discordant_processed_df <- purrr::map_dfr(discordant_row_count_files, process_discordant_valid_files)

recombination_df <- left_join(concordant_processed_df, discordant_processed_df, by = c("strain", "sample", "experiment")) %>%
  select(strain, sample, experiment, concordant_count, discordant_count)

recombination_df$concordant_count <- as.numeric(recombination_df$concordant_count)
recombination_df$discordant_count <- as.numeric(recombination_df$discordant_count)

log_step("Calculation recombination rate...")


complete_recombination_df <- left_join(recombination_df, ratio_control_file) %>%  
  mutate(discordant_count_norm = discordant_count / Ratio_vs_reference,
         discordant_concordant_ratio = discordant_count / concordant_count,
         discordant_concordant_ratio_norm = discordant_concordant_ratio / Ratio_vs_reference) %>% 
  mutate(concordant_percentage = concordant_count / (concordant_count + discordant_count_norm) * 100,
         discordant_percentage = discordant_count_norm / (concordant_count + discordant_count_norm) * 100) %>%
  select(strain, sample, experiment, concordant_count, discordant_count, Ratio_vs_reference, concordant_percentage, discordant_percentage,
         discordant_count_norm, discordant_concordant_ratio, discordant_concordant_ratio_norm)


write_tsv(complete_recombination_df, file.path(root_dir, paste0(strain_name, "_recombination_df.tsv")))



log_step("Calculation recombination summary...")


recombination_summary <- complete_recombination_df %>%  group_by(strain, sample) %>% 
  summarise(mean_concordant_percentage = mean(concordant_percentage, na.rm = TRUE), 
            sd_concordant_percentage = sd(concordant_percentage, na.rm = TRUE),
            mean_discordant_percentage = mean(discordant_percentage, na.rm = TRUE),
            sd_discordant_percentage = sd(discordant_percentage, na.rm = TRUE),
            mean_discordant_concordant_ratio_norm = mean(discordant_concordant_ratio_norm, na.rm = TRUE),
            sd_discordant_concordant_ratio_norm = sd(discordant_concordant_ratio_norm, na.rm = TRUE)) %>%
  select(strain, sample, mean_discordant_percentage, sd_discordant_percentage,
         mean_discordant_concordant_ratio_norm, sd_discordant_concordant_ratio_norm)

write_tsv(recombination_summary, file.path(root_dir, paste0(strain_name, "_recombination_summary.tsv")))



samples_order <- c("TSG", "TLG", "TLR")

log_step("Plotting...")
# Generate plots



recombination_summary_graphs <- recombination_summary %>% ungroup() %>% filter(sample != "T0") %>% 
  pivot_longer(cols = c(
    mean_discordant_percentage, sd_discordant_percentage), 
    names_to = "category",
    values_to = "value") %>% 
  separate(category, into = c("metric", "type", "trash"), sep = "_") %>% 
  select(-trash) %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  ) %>% mutate(sample = factor(sample, levels = samples_order))





log_step("Plotting...")
# Generate plots

bar_plot <- ggplot(recombination_summary_graphs, aes(x = sample, y = mean, fill = type)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                linewidth = 0.8, width = 0.2, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("grey")) +
  theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                               panel.background = element_blank(), 
                                               plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  coord_cartesian(ylim = c(0, 0.2), expand=FALSE) +
  theme(aspect.ratio = 0.75) + 
  scale_x_discrete(name = expression("Sample")) +
  scale_y_continuous(name = expression("Percentage"),
                     #limits = c(0, 1),
                     breaks = seq(0,0.2,0.1)) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 20, angle = 0), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  labs(
    title = paste0("Recombination rate - ", strain_name)
  )

ggsave(
  filename = paste0(strain,"/", "Recombination_rate_plot_", strain_name, ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)


recombination_summary_ratio_graphs <- recombination_summary %>% ungroup() %>% filter(sample != "T0") %>% 
  pivot_longer(cols = c(
    mean_discordant_concordant_ratio_norm, sd_discordant_concordant_ratio_norm), 
    names_to = "category",
    values_to = "value") %>% 
  separate(category, into = c("metric", "type", "trash"), sep = "_") %>% 
  select(-trash) %>%
  pivot_wider(
    names_from = metric,
    values_from = value
  ) %>% mutate(sample = factor(sample, levels = samples_order))





log_step("Plotting...")
# Generate plots

bar_plot <- ggplot(recombination_summary_ratio_graphs, aes(x = sample, y = mean, fill = type)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), 
                linewidth = 0.8, width = 0.2, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("grey")) +
  theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                               panel.background = element_blank(), 
                                               plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  coord_cartesian(ylim = c(0, 0.002), expand=FALSE) +
  theme(aspect.ratio = 0.75) + 
  scale_x_discrete(name = expression("Sample")) +
  scale_y_continuous(name = expression("Ratio"),
                     #limits = c(0, 1),
                     breaks = seq(0,0.002,0.001)) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 20, angle = 0), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  labs(
    title = paste0("Recombination rate ratio - ", strain_name)
  )

ggsave(
  filename = paste0(strain,"/", "Recombination_rate_ratio_plot_", strain_name, ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)