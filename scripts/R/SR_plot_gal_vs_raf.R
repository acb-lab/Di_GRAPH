#### R script to plot gal_vs raf changes, 75nt, average all experiments
### Loop for each strain/ file in MYWD
### 23/04/2026 - Lydia



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
library(svglite)
library(purrr)
library(stringr)
library(readr)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)
library(ggrepel)





# Define a function to process tsv files
process_coverage_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  # directory name above the file
  
  # Read file
  temp_file <- read_tsv(file_path, col_names = TRUE, show_col_types = FALSE) 
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  # Prepare coverage dataframe
  coverage_file <- temp_file %>%  select(Feature_name_A, Category_A, avg_ratio_coverage_feature, sd_ratio_coverage_feature, Analysis) %>% 
    mutate(
      strain = strain_name
    )
  coverage_file_renamed <- coverage_file %>%  
    rename_with(~ paste0("avg_ratio_coverage_feature_", unique(coverage_file$Analysis)), .cols = "avg_ratio_coverage_feature") %>%
    rename_with(~ paste0("sd_ratio_coverage_feature_", unique(coverage_file$Analysis)), .cols = "sd_ratio_coverage_feature") %>%  select(!c(Analysis))
  
  return(coverage_file_renamed)
}


log_step("Finding coverage files...")
# Get all coverage.tsv files recursively in root folder
TLG_coverage_files <- list.files(
  path = root_dir,
  pattern = "Genomic_sorted_.*_T0vsTLG\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# Get all coverage.tsv files recursively in root folder
TLR_coverage_files <- list.files(
  path = root_dir,
  pattern = "Genomic_sorted_.*_T0vsTLR\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

# Get all coverage.tsv files recursively in root folder
TSG_coverage_files <- list.files(
  path = root_dir,
  pattern = "Genomic_sorted_.*_T0vsTSG\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing coverage files...")
TLG_processed_df <- purrr::map_dfr(TLG_coverage_files, process_coverage_files)
TLR_processed_df <- purrr::map_dfr(TLR_coverage_files, process_coverage_files)
TSG_processed_df <- purrr::map_dfr(TSG_coverage_files, process_coverage_files)

complete_df <- left_join(TLG_processed_df, TLR_processed_df) %>% 
  left_join(., TSG_processed_df)

# complete_df 
# write_tsv(complete_df, file.path(root_dir, paste0(strain_name, "_complete_df.tsv")))

complete_df_ratio <- complete_df %>% select(strain, Feature_name_A, Category_A, 
                                            avg_ratio_coverage_feature_T0vsTLG, avg_ratio_coverage_feature_T0vsTLR, avg_ratio_coverage_feature_T0vsTSG) %>% 
  filter(if_all(starts_with("avg_ratio_coverage_feature"), ~ !is.na(.) & is.finite(.))) %>%
  mutate(UP_DOWN_TLG = ifelse(avg_ratio_coverage_feature_T0vsTLG > 1.2 & (avg_ratio_coverage_feature_T0vsTLR > 0.8 & avg_ratio_coverage_feature_T0vsTLR <1.2), "UP", 
                              ifelse(avg_ratio_coverage_feature_T0vsTLG < 0.8 & (avg_ratio_coverage_feature_T0vsTLR > 0.8 & avg_ratio_coverage_feature_T0vsTLR <1.2), "DOWN", "NO_CHANGE"))) %>% 
  mutate(UP_DOWN_TSG = ifelse(avg_ratio_coverage_feature_T0vsTSG > 1.2 & (avg_ratio_coverage_feature_T0vsTLR > 0.8 & avg_ratio_coverage_feature_T0vsTLR <1.2), "UP", 
                              ifelse(avg_ratio_coverage_feature_T0vsTSG < 0.8 & (avg_ratio_coverage_feature_T0vsTLR > 0.8 & avg_ratio_coverage_feature_T0vsTLR <1.2), "DOWN", "NO_CHANGE")))

# write_tsv(complete_df_ratio, file.path(root_dir, paste0(strain_name, "_complete_df_ratio.tsv")))

n_features <- nrow(complete_df_ratio)

UP_DOWN_summary <- complete_df_ratio %>% group_by(strain) %>% 
  summarise(UP_TLG = sum(UP_DOWN_TLG == "UP"),
            DOWN_TLG = sum(UP_DOWN_TLG == "DOWN"),
            UP_TSG = sum(UP_DOWN_TSG == "UP"),
            DOWN_TSG = sum(UP_DOWN_TSG == "DOWN")) %>% 
  mutate(UP_TLG_perc = UP_TLG/n_features*100,
         DOWN_TLG_perc = DOWN_TLG/n_features*100,
         UP_TSG_perc = UP_TSG/n_features*100,
         DOWN_TSG_perc = DOWN_TSG/n_features*100) %>% 
  select(strain, ends_with("_perc"))



write_tsv(UP_DOWN_summary, file.path(root_dir, paste0(strain_name, "_UP_DOWN_summary.tsv")))

# 
categoria_colors <- c(
  "ORF" = "#4F71BE", 
  "intergenic" = "#FF051D", 
  "long_terminal_repeat" = "#DE8344", 
  "transposable_element_gene" = "#A5A5A5", 
  "LTR_retrotransposon" = "#F5C242", 
  "tRNA_gene" = "#6A99D0", 
  "rRNA_gene" = "#7EAB55", 
  "ncRNA_gene" = "#2D4374", 
  "snRNA_gene" = "#934D20", 
  "snoRNA_gene" = "#636363", 
  "ARS" = "#937424", 
  "centromere" = "#355D8D", 
  "telomere" = "#4B6733"
)
log_step("Plotting...")
plot_TLG<- ggplot(complete_df_ratio, aes(x = avg_ratio_coverage_feature_T0vsTLG, y =avg_ratio_coverage_feature_T0vsTLR)) +
  geom_point(aes(colour = Category_A, alpha = UP_DOWN_TLG),size = 2, shape = 16) + #0.5
  scale_colour_manual(values=categoria_colors) +
  scale_alpha_manual(values = c("NO_CHANGE" = 0.8, "UP" = 0.8, "DOWN" = 0.8)) +
  theme_classic(base_family = "Helvetica") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                   panel.background = element_blank(), 
                                                   plot.background = element_rect(fill = "transparent", colour = NA)) +
  # geom_text_repel(data = subset(complete_df_ratio, UP_DOWN_TLG == "UP" & ratio_medio_T0vsTLG > 2),
  #               aes(label = Nombre, color = categoria),
  #               size = 3, show.legend = FALSE, max.overlaps = Inf) +
  # geom_text_repel(data = subset(complete_df_ratio, UP_DOWN_TLG == "DOWN" & ratio_medio_T0vsTLG < 0.5),
  #                 aes(label = Nombre, color = categoria),
  #                 size = 3, show.legend = FALSE, max.overlaps = Inf) +
  theme(legend.position="none") + 
  coord_cartesian(xlim = c(0,5.3), ylim = c(0, 5.3), expand=FALSE) +
  theme(aspect.ratio = 0.8) +
  geom_vline(xintercept=c(0.75, 1.25),
             linetype="dashed", color = "black", linewidth=0.3) +
  geom_hline(yintercept=c(0.75, 1.25),
             linetype="dashed", color = "black", linewidth=0.3) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 20), #25
        axis.title.y = element_text(vjust = 1, size = 20)) + #25
  theme(axis.text.x = element_text(vjust = 0, size = 15), #20
        axis.text.y = element_text(vjust = 0, size = 15)) +
  labs(
    title = paste0("Gal vs Raf coverage - ", strain_name, " - ", "TLG"),
    x = "Gal_coverage",
    y = "Raf_coverage"
  )

# Save SVG
ggsave(
  filename = file.path(root_dir, paste0("Gal_vs_raf_coverage_", strain_name, "_","TLG", ".svg")),
  plot = plot_TLG,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)

plot_TSG <- ggplot(complete_df_ratio, aes(x = avg_ratio_coverage_feature_T0vsTSG, y = avg_ratio_coverage_feature_T0vsTLR)) +
  geom_point(aes(colour = Category_A, alpha = UP_DOWN_TSG),size = 2, shape = 16) + #0.5
  scale_colour_manual(values=categoria_colors) +
  scale_alpha_manual(values = c("NO_CHANGE" = 0.8, "UP" = 0.8, "DOWN" = 0.8)) +
  theme_classic(base_family = "Helvetica") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                   panel.background = element_blank(), 
                                                   plot.background = element_rect(fill = "transparent", colour = NA)) +
  # geom_text_repel(data = subset(complete_df_ratio, UP_DOWN_TSG == "UP" & ratio_medio_T0vsTSG > 2),
  #               aes(label = Nombre, color = categoria),
  #               size = 3, show.legend = FALSE, max.overlaps = Inf) +
  # geom_text_repel(data = subset(complete_df_ratio, UP_DOWN_TSG == "DOWN" & ratio_medio_T0vsTSG < 0.5),
  #                 aes(label = Nombre, color = categoria),
  #                 size = 3, show.legend = FALSE, max.overlaps = Inf) +
  theme(legend.position="none") + 
  coord_cartesian(xlim = c(0,5.3), ylim = c(0, 5.3), expand=FALSE) +
  theme(aspect.ratio = 0.8) +
  geom_vline(xintercept=c(0.75, 1.25),
             linetype="dashed", color = "black", linewidth=0.3) +
  geom_hline(yintercept=c(0.75, 1.25),
             linetype="dashed", color = "black", linewidth=0.3) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 20), #25
        axis.title.y = element_text(vjust = 1, size = 20)) + #25
  theme(axis.text.x = element_text(vjust = 0, size = 15), #20
        axis.text.y = element_text(vjust = 0, size = 15)) +
  labs(
    title = paste0("Gal vs Raf coverage - ", strain_name, " - ", "TSG"),
    x = "Gal_coverage",
    y = "Raf_coverage"
  )

# Save SVG
ggsave(
  filename = file.path(root_dir, paste0("Gal_vs_raf_coverage_", strain_name, "_TSG", ".svg")),
  plot = plot_TSG,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)