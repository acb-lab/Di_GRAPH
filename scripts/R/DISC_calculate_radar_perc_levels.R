#### R script to plot radar plots (% CategoryA) (blast levels of validation)
### Loop for each strain/blast_option in MYWD
### 15/04/2026 - Lydia

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
  log_step("Usage: script.R <root_dir> <strain> <wd_dir> <blast_option>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]
wd_dir   <- args[3]
blast_option   <- args[4]


strain <- sub("/$", "", strain)
strain_name <- basename(strain)


log_step(paste("ROOT_DIR:", root_dir))
log_step(paste("STRAIN:", strain))

log_step(paste("BLAST_OPTION:", blast_option))

# load libraries

library(ggplot2)
library(svglite)
library(extrafont)
library(purrr)
library(stringr)
library(readr)
library(fmsb)
#library(ggbreak)
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
    mutate(blast_option = blast_option) %>% 
    group_by(strain, sample, experiment, blast_option,  Category_A) %>%
    summarise(total_count = sum(count), .groups = "drop_last") %>%
    mutate(
      group_total = sum(total_count),
      percent = 100 * (total_count / group_total)
    ) %>%
    ungroup()
  
  return(counts_distribution_file)
}


log_step("Finding valid_counts files...")
# Get all _inter_discordant_pairs_unique_processed_valid_counts files recursively in root folder

valid_counts_files <- list.files(
  path = root_dir,
  pattern = paste0(blast_option,"_inter_discordant_pairs_unique_processed_valid_counts\\.tsv$"),
  recursive = TRUE,
  full.names = TRUE
)


log_step("Processing valid_counts files...")
valid_counts_processed_df <- purrr::map_dfr(valid_counts_files, process_valid_counts_files)



log_step("Finding ratio control file...")

ratio_control_file <- read_tsv(file.path(wd_dir, "control_count_ratio.tsv"), col_names = TRUE)



complete_valid_counts_df <- left_join(valid_counts_processed_df, ratio_control_file) %>%  
  mutate(total_count_norm = total_count / Ratio_vs_reference)

genomic_categories_order <- c(
  "ORF", "intergenic", "long_terminal_repeat", "transposable_element_gene",
  "LTR_retrotransposon", "tRNA_gene", "rRNA_gene", "ncRNA_gene",
  "snRNA_gene", "snoRNA_gene", "ARS", "centromere", "telomere"
)

sample_order <- c("TSG" , "TLG", "TLR")


complete_valid_counts_summary_df <- complete_valid_counts_df %>% 
  group_by(strain, sample, Category_A) %>% 
  summarise(mean_total_count_norm = mean(total_count_norm)) %>% 
  mutate(Category_A = factor(Category_A, levels = genomic_categories_order),
         sample = factor(sample, levels = sample_order)) %>% 
  arrange(strain, sample, Category_A) %>% 
  ungroup()





######
# Define the transformation function
transform_y <- function(y) {
  ifelse(
    y <= 100,
    y * (0.75 / 100),                     # scale 0-100 to 0-0.75
    0.75 + ((y - 100) * (0.25 / (5000 - 100))) # scale 100-5000 to 0.25-1
  )
}

transform_y_50 <- function(y) {
  ifelse(
    y <= 50,
    y * (0.75 / 50),                     # scale 0-50 to 0-0.75
    0.75 + ((y - 50) * (0.25 / (5000 - 50))) # scale 50-5000 to 0.25-1
  )
}





# Add a transformed y column
complete_valid_counts_summary_df$mean_total_count_norm_trans <- transform_y(complete_valid_counts_summary_df$mean_total_count_norm)
complete_valid_counts_summary_df$mean_total_count_norm_trans_50 <- transform_y_50(complete_valid_counts_summary_df$mean_total_count_norm)


complete_valid_counts_summary_df_radar <- complete_valid_counts_summary_df %>% 
  select(Category_A, mean_total_count_norm_trans, sample) %>% 
  pivot_wider(names_from = "Category_A", values_from = "mean_total_count_norm_trans") %>% 
  as.data.frame() %>% 
  select(!c("sample"))

rownames(complete_valid_counts_summary_df_radar) <- c("TSG", "TLG", "TLR")

complete_valid_counts_summary_df_radar_50 <- complete_valid_counts_summary_df %>% 
  select(Category_A, mean_total_count_norm_trans_50, sample) %>% 
  pivot_wider(names_from = "Category_A", values_from = "mean_total_count_norm_trans_50") %>% 
  as.data.frame() %>% 
  select(!c("sample"))

rownames(complete_valid_counts_summary_df_radar_50) <- c("TSG", "TLG", "TLR")




max_min <- data.frame(
  ORF = c(1, 0), intergenic = c(1, 0), long_terminal_repeat = c(1, 0),
  transposable_element_gene = c(1, 0), LTR_retrotransposon = c(1, 0), tRNA_gene = c(1, 0),
  rRNA_gene = c(1, 0), ncRNA_gene = c(1, 0), snRNA_gene = c(1, 0),
  snoRNA_gene = c(1, 0), ARS = c(1, 0), centromere = c(1, 0),
  telomere = c(1, 0)
)



rownames(max_min) <- c("Max", "Min")


# Bind the variable ranges to the data
df_radar <- rbind(max_min, complete_valid_counts_summary_df_radar)
df_radar_50 <- rbind(max_min, complete_valid_counts_summary_df_radar_50)




log_step("Plotting...")
# Generate plots

svglite::svglite(file = paste0(strain,"/", "Radar_plot_", strain_name, "_", blast_option, "_TSG_TLG_TLR", ".svg"), width = 8, height = 8)
radarchartcirc(df_radar, axistype = 1,
               # Customize the polygon
               pcol = c("#252159", "#469CD7", "#8BE0FC"), 
               seg = 4,
               pty = 32, #32
               #pfcol = FALSE, 
               #pfcol = scales::alpha("black", 0.0), 
               plwd = 1.5, plty = 1,
               # Customize the grid
               cglcol = "grey", cglty = 2, cglwd = 0.8,
               # Customize the axis
               axislabcol = "grey9", calcex = 0.6,
               title = paste0("radar plot", "-", strain, "-", blast_option),
               # Variable labels
               vlcex = 0.7, vlabels = colnames(df_radar),
               caxislabels = c(0, 33.3, 66.6, 100, 5000)
)
dev.off()

svglite::svglite(file = paste0(strain,"/", "Radar_plot_", strain_name, "_", blast_option, "_50_TSG_TLG_TLR", ".svg"), width = 8, height = 8)
radarchartcirc(df_radar_50, axistype = 1,
               # Customize the polygon
               pcol = c("#252159", "#469CD7", "#8BE0FC"), 
               seg = 4,
               pty = 32, #32
               #pfcol = FALSE, 
               #pfcol = scales::alpha("black", 0.0), 
               plwd = 1.5, plty = 1,
               # Customize the grid
               cglcol = "grey", cglty = 2, cglwd = 0.8,
               # Customize the axis
               axislabcol = "grey9", calcex = 0.6,
               title = paste0("radar plot_50", "-", strain, "-", blast_option),
               # Variable labels
               vlcex = 0.7, vlabels = colnames(df_radar_50),
               caxislabels = c(0, 16.6, 33.3, 50, 5000)
)
dev.off()