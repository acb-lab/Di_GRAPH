#### R script to calculate valid inter-chromosomal discordant pairs matrix (merging all experiments, blast levels of validation)
### Loop for each strain/sample in MYWD
### 14/04/2026 - Lydia

log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}

# -------------------------
# Read command line args
# -------------------------
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 5) {
  log_step("ERROR: Incorrect number of arguments")
  log_step("Usage: script.R <root_dir> <strain> <sample> <blast_option> <category_path>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]
sample   <- args[3]
blast_option   <- args[4]
category_path   <- args[5]

strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))
log_step(paste("SAMPLE:", sample))
log_step(paste("BLAST OPTION:", blast_option))


# load libraries

library(readr)
library(extrafont)
library(stringr)
library(svglite)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)


path_to_features_pairs_file <- file.path(category_path, "PMV_features_pairs.tsv")




log_step("Finding possible pairs...")

# Inter-chromosomal discordant read pairs matrix
df_features_pairs <- read_tsv(path_to_features_pairs_file, col_names = TRUE, show_col_types = FALSE)

log_step("Loading valid reads...")
valid_reads_files <- list.files(
  path = root_dir,
  pattern = paste0(sample, ".*_", blast_option, ".*inter_discordant_pairs_unique_processed_valid\\.tsv$"),
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing valid reads files...")
valid_pairs_df <- map_dfr(valid_reads_files, read_tsv)

log_step("Loading control reads...")
control_reads_files <- list.files(
  path = root_dir,
  pattern = paste0(sample, ".*inter_discordant_pairs_unique_processed_control\\.tsv$"),
  recursive = TRUE,
  full.names = TRUE
)
control_reads_files
log_step("Processing control reads files...")
control_df <- map_dfr(control_reads_files, read_tsv )%>% 
  filter(Category_A=="control_norm", Category_B=="control_norm")

get_discordant_matrix <- function(valid_pairs_df_file, all_pairs_df, control_df_file) {
  posA_info <- all_pairs_df %>%
    select(Feature_name_A, Position_A = Position, Essential_A = Essential) %>%
    distinct()
  
  posB_info <- all_pairs_df %>%
    select(Feature_name_B, Position_B = Position, Essential_B = Essential) %>%
    distinct()
  
  valid_pairs_df_posA <- valid_pairs_df_file %>%
    left_join(posA_info, by = "Feature_name_A")
  
  valid_pairs_df_posAB <- valid_pairs_df_posA %>%
    left_join(posB_info, by = "Feature_name_B")
  
  pair_counts <- valid_pairs_df_posAB %>%
    group_by(strain, sample, Feature_name_A, Feature_name_B) %>%
    summarise(count = n(), .groups = "drop")
  
  complete_matrix <- left_join(valid_pairs_df_posAB, pair_counts,
                               by = c("Feature_name_A", "Feature_name_B", "strain", "sample")) %>%
    distinct()
  
  count_total <- nrow(control_df_file)
  complete_matrix <- complete_matrix %>%
    mutate(count_norm = (count / count_total) * 100) %>% 
    select(!c(pair_group, pair_group_name, number, Feature_name_B_prev, Feature_name_B_next))
  
  return(complete_matrix)
}

log_step("Generating discordant matrix...")
complete_matrix <- get_discordant_matrix(valid_pairs_df, df_features_pairs, control_df)
complete_matrix_noORF_nointergenic <- complete_matrix %>% filter(Category_A != "ORF", Category_B != "ORF", Category_A != "intergenic", Category_B != "intergenic")

log_step("Saving discordant matrix...")
write_tsv(complete_matrix, file = paste0(strain, "/", sample, "_", blast_option, "_inter_discordant_pairs_unique_processed_valid_discordant_matrix_allexperiments.tsv"))
write_tsv(complete_matrix_noORF_nointergenic, file = paste0(strain, "/", sample, "_", blast_option, "_inter_discordant_pairs_unique_processed_valid_discordant_matrix_noORF_nointergenic_allexperiments.tsv"))

log_step("Plotting discordant matrix...")
matrix_plot <-ggplot(complete_matrix, aes(y = Position_B, x = Position_A)) +
  #geom_hdr(xlim = c(0,14518), ylim = c(14518,0), method = "kde", fill = "brown4") + 
  geom_point(aes(colour = count_norm), alpha = 1, size = 0.8) +
  scale_colour_gradient2(low= "white", mid = "#4ae034", high = "#00641b",
                         midpoint = 1,
                         limits = c (0, 5),
                         oob = scales::squish) + 
  theme_bw(base_family = "Arial") + 
  theme(panel.background = element_blank()) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(legend.position="right") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  guides (colour = guide_colourbar(barwidth = 0.5, barheight = 5,
                                   frame.colour = "black", frame.linewidth = 0.25,
                                   ticks.colour = NA)) + 
  coord_cartesian(xlim = c(0,14518), ylim = c(14518, 0), expand=FALSE) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                                6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                                13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487), position = "top") +
  scale_y_continuous(breaks = c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                                6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                                13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487)) +
  geom_hline(yintercept=c(6570, 13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="solid", color = "black", linewidth=0.05) +
  geom_hline(yintercept=c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                          6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                          13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 144877),
             linetype="dashed", color = "black", linewidth=0.05) +
  geom_vline(xintercept=c(6570, 13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="solid", color = "black", linewidth=0.05) +
  geom_vline(xintercept=c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                          6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                          13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="dashed", color = "black", linewidth=0.05) +
  theme(axis.text.y=element_text(size=0)) +
  theme(axis.text.x=element_text(size=0)) +
  theme(axis.title.x = element_text(size=0)) +
  theme(axis.title.y = element_text(size=0)) +
  theme(axis.ticks = element_blank()) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "solid", linewidth = 0.05) +
  labs(
    title = paste0("Inter_chromosomal discordant matrix - ", strain_name, " - ", sample, "-", blast_option),
    colour = "Freq (%)")


ggsave(
  filename = paste0(strain,"/", "Inter_chromosomal_discordant_matrix_plot_", strain_name, "_", sample, "_", blast_option, ".svg"),
  plot = matrix_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)

log_step("Plotting discordant matrix - no ORF no intergenic...")
matrix_plot <-ggplot(complete_matrix_noORF_nointergenic, aes(y = Position_B, x = Position_A)) +
  #geom_hdr(xlim = c(0,14518), ylim = c(14518,0), method = "kde", fill = "brown4") + 
  geom_point(aes(colour = count_norm), alpha = 1, size = 2) +
  scale_colour_gradient2(low= "white", mid = "#4ae034", high = "#00641b",
                         midpoint = 1,
                         limits = c (0, 5),
                         oob = scales::squish) + 
  theme_bw(base_family = "Arial") + 
  theme(panel.background = element_blank()) +
  theme(plot.background = element_rect(fill = "transparent", colour = NA))+
  theme(legend.position="right") +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  guides (colour = guide_colourbar(barwidth = 0.5, barheight = 5,
                                   frame.colour = "black", frame.linewidth = 0.25,
                                   ticks.colour = NA)) + 
  coord_cartesian(xlim = c(13195,14518), ylim = c(14518, 13195), expand=FALSE) +
  theme(aspect.ratio = 1) +
  scale_x_continuous(breaks = c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                                6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                                13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487), position = "top") +
  scale_y_continuous(breaks = c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                                6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                                13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487)) +
  geom_hline(yintercept=c(6570, 13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="solid", color = "black", linewidth=0.05) +
  geom_hline(yintercept=c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                          6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                          13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 144877),
             linetype="dashed", color = "black", linewidth=0.05) +
  geom_vline(xintercept=c(6570, 13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="solid", color = "black", linewidth=0.05) +
  geom_vline(xintercept=c(118, 574, 755, 1591, 1914, 2053, 2636, 2957, 3198, 3596, 3944, 4522, 5027, 5462, 6059,
                          6570, 6683, 7134, 7337, 8161, 8493, 8638, 9245, 9562, 9798, 10201, 10568, 11128, 11648, 12080, 12685,
                          13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487),
             linetype="dashed", color = "black", linewidth=0.05) +
  theme(axis.text.y=element_text(size=0)) +
  theme(axis.text.x=element_text(size=0)) +
  theme(axis.title.x = element_text(size=0)) +
  theme(axis.title.y = element_text(size=0)) +
  theme(axis.ticks = element_blank()) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "solid", linewidth = 0.05) +
  labs(
    title = paste0("Inter_chromosomal discordant matrix no ORF no intergenic - ", strain_name, " - ", sample, "-", blast_option),
    colour = "Freq (%)")


ggsave(
  filename = paste0(strain,"/", "Inter_chromosomal_discordant_matrix_reduced_plot_", strain_name, "_", sample, "_", blast_option, ".svg"),
  plot = matrix_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)
