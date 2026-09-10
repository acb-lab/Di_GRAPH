#### R script to identify recombination hotspots (in all experiments)
### Loop for each strain/sample in MYWD
### 12/04/2026 - Lydia

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
  log_step("Usage: script.R <root_dir> <strain> <sample> <path_to_freqs>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]
sample   <- args[3]
path_to_freqs   <- args[4]

strain <- sub("/$", "", strain)
strain_name <- basename(strain)


log_step(paste("ROOT_DIR:", root_dir))
log_step(paste("STRAIN:", strain))
log_step(paste("SAMPLE:", sample))


# load libraries

library(readr)
library(stringr)
library(svglite)
library(extrafont)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
# Suppress summarise info
options(dplyr.summarise.inform = FALSE)



log_step("Loading matrix files...")
matrix_files <- list.files(
  path = root_dir,
  pattern = paste0(sample, ".*inter_discordant_pairs_unique_processed_valid_discordant_matrix\\.tsv$"),
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing matrix files...")
matrix_files
matrix_df <- map_dfr(matrix_files, read_tsv) %>% 
  select(Feature_name_A, Category_A, Chromosome_A, Essential_A,
         Feature_name_B, Category_B, Chromosome_B, Essential_B,
         strain, sample, experiment, Read_name_ID)

log_step("Finding recombination hotspots...")
hotspots_feature_A <- matrix_df %>% group_by(Feature_name_A) %>%
  filter(n_distinct(experiment) == n_distinct(.$experiment)) %>%
  ungroup() %>% 
  select(1, 2, 3, 4) %>% rename("Category_A" = !!names(.[2]),
                                    "Chromosome_A" = !!names(.[3]),
                                    "Essential_A" = !!names(.[4])) %>% 
  unique()





log_step("Loading freq file...")
# Load control reads
freq_1_2_3 <- read_tsv(path_to_freqs, col_names = TRUE) %>% select("Feature_name_A", "Position_A", "count_norm") %>% 
  group_by(Feature_name_A) %>% 
  mutate(global_freq = sum(count_norm)) %>%  select(!c("count_norm")) %>%  unique()

hotspots_freq <- merge(hotspots_feature_A, freq_1_2_3, by = c("Feature_name_A")) %>% unique()

log_step("Saving hotspots dataframe...")
write_tsv(hotspots_freq, file = paste0(strain, "/", sample, "_inter_discordant_pairs_unique_processed_valid_discordant_hotspots.tsv"))

log_step("Plotting hotspots...")
hotspots_plot <- ggplot(hotspots_freq, aes(y = Position_A, x = 0)) +
  #geom_segment(aes(xend = 0.95, yend = Position_A, color = global_freq), alpha = 1)+
  geom_tile(aes(color = global_freq), alpha = 1)+
  scale_color_gradientn(
    colors = c("white", "#e31a1c", "#8b2500"),
    #values = scales::rescale(c(0, 25, 50)),
    na.value = "gray90",
    #limits = c(0, 50),
    limits = c(0, 25),
    oob = scales::squish) +
  guides (color = guide_colourbar(barwidth = 0.5, barheight = 10,
                                  frame.colour = "black", frame.linewidth = 0.25,
                                  ticks.colour = NA)) + 
  labs(color = "%") +
  theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                               panel.background = element_blank(), 
                                               plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.position="right") +
  #coord_cartesian(ylim = c(13195,14518), expand=FALSE) +
  #coord_cartesian(ylim = c(1, 14519), xlim = c(0,1), expand=FALSE) +
  coord_cartesian(ylim = c(1, 14519), xlim = c(0,0.1), expand=FALSE) +
  theme(aspect.ratio = 20) +
  scale_y_continuous(breaks = c(1, 6570,
                                13195, 13578, 13669, 13719, 13994, 14019, 14036, 14042, 14119, 14471, 14487)) +
  theme(axis.text.y=element_text(size=0)) +
  theme(axis.text.x=element_text(size=0)) +
  theme(axis.title.x = element_text(size=0)) +
  theme(axis.title.y = element_text(size=0)) +
  theme(axis.ticks.x = element_blank()) +
  labs(
    title = paste0("Hotspots - ", strain_name, " - ", sample)
  )


ggsave(
  filename = paste0(strain,"/", "Inter_chromosomal_discordant_hotspots_plot_", strain_name, "_", sample, ".svg"),
  plot = hotspots_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)