#### R script to calculate recombination hotspots distribution (in all strains, blast levels of validation)
###  MYWD
### 15/04/2026 - Lydia

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



# load libraries

library(ggplot2)
library(svglite)
library(purrr)
library(stringr)
library(readr)
library(scales)
library(extrafont)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)



TSG_hotspots_files <- list.files(
  path = root_dir,
  pattern ="TSG_option1_inter_discordant_pairs_unique_processed_valid_discordant_hotspots\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


TLG_hotspots_files <- list.files(
  path = root_dir,
  pattern ="TLG_option1_inter_discordant_pairs_unique_processed_valid_discordant_hotspots\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

TLR_hotspots_files <- list.files(
  path = root_dir,
  pattern ="TLR_option1_inter_discordant_pairs_unique_processed_valid_discordant_hotspots\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

all_TSG_valid_hotspots_files <- lapply(TSG_hotspots_files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  option <- str_extract(basename(f), "option[0-9]+")
  strain_name <- basename(dirname(f))
  df <- df %>% mutate(blast_option = option,
                      strain = strain_name,
                      sample = "TSG") %>% 
    filter(global_freq > 4)
  
  return(df)
}) %>%
  bind_rows()  

all_TLG_valid_hotspots_files <- lapply(TLG_hotspots_files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  option <- str_extract(basename(f), "option[0-9]+")
  strain_name <- basename(dirname(f))
  df <- df %>% mutate(blast_option = option,
                      strain = strain_name,
                      sample = "TLG") %>% 
    filter(global_freq > 4)
  
  return(df)
}) %>%
  bind_rows()  


all_TLR_valid_hotspots_files <- lapply(TLR_hotspots_files, function(f) {
  df <- read_tsv(f, show_col_types = FALSE)
  option <- str_extract(basename(f), "option[0-9]+")
  strain_name <- basename(dirname(f))
  df <- df %>% mutate(blast_option = option,
                      strain = strain_name,
                      sample = "TLR") %>% 
    filter(global_freq > 4)
  
  return(df)
}) %>%
  bind_rows()  



all_hotspots <- bind_rows(all_TSG_valid_hotspots_files,
                          all_TLG_valid_hotspots_files,
                          all_TLR_valid_hotspots_files)

sample_order <- c("TSG", "TLG", "TLR")

all_hotspots_ordered <- all_hotspots %>% mutate(sample = factor(sample, levels = sample_order)) 


p <- ggplot(all_hotspots_ordered, aes(x = global_freq, color = sample, fill = sample)) +
  geom_histogram(alpha = 0.25, position = "identity") +
  coord_cartesian(xlim = c(0,70), ylim = c(0, 60), expand=FALSE) +
  theme_classic(base_family = "Arial") +
  theme(panel.grid = element_line(color = "black", linewidth = 0.1),
        panel.background = element_blank(),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(panel.background = element_rect(fill = "white")) +
  theme(legend.position="right") +
  labs(title = paste0("Hotspots distribution")) +
  facet_wrap(sample~ strain, nrow = 3, ncol = 5) 

ggsave(
  filename = paste0(root_dir,"/", "Inter_chromosomal_discordant_hotspots_distribution.svg"),
  plot = p,
  width = 12,
  height = 8,
  device = svglite,
  bg = "transparent"
)