#### R script to calculate valid inter-chromosomal discordant pairs global distribution (blast levels of validation)
### Average from all experiments
### Loop for each strain/sample/blast_option in MYWD
### 14/04/2026 - Lydia

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
  log_step("Usage: script.R <root_dir> <strain> <sample> <blast_option>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]
sample   <- args[3]
blast_option   <- args[4]

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




# Functions
log_step <- function(message) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  message(sprintf("[%s] %s", timestamp, message))
}




log_step("Finding inter_discordant_pairs_unique_processed_valid_global_distribution files...")
# Get all coverage.tsv files recursively in root folder
discordant_reads_count_files <- list.files(
  path = root_dir,
  pattern = paste0(sample, ".*_", blast_option, ".*inter_discordant_pairs_unique_processed_valid_global_distribution\\.tsv$"),
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing inter_discordant_pairs_unique_processed_valid_global_distribution files...")
discordant_reads_count <- map_dfr(discordant_reads_count_files, read_tsv)


# Set genomic categories order
genomic_categories_order <- c(
  "ORF", "intergenic", "long_terminal_repeat", "transposable_element_gene",
  "LTR_retrotransposon", "tRNA_gene", "rRNA_gene", "ncRNA_gene",
  "snRNA_gene", "snoRNA_gene", "ARS", "centromere", "telomere"
)


# Calculate average global inter-chromosomal discordant read pairs distribution 
global_distribution <- discordant_reads_count %>% 
  group_by(strain, sample, strain_sample_comb,  Category_A, Category_B) %>%
  summarise(all_mean_global_percentage = mean(mean_global_percentage, na.rm = TRUE),
            all_sd_global_percentage = sd(mean_global_percentage, na.rm = TRUE)) %>%
  # Apply genomic_categories_order
  mutate(Category_A = factor(Category_A, levels = genomic_categories_order),
         Category_B = factor(Category_B, levels = genomic_categories_order)) %>%
  arrange(strain, sample, Category_A, Category_B) %>% 
  ungroup()

heatmap <- ggplot(global_distribution, aes(x = Category_B, y = Category_A, fill =all_mean_global_percentage)) +
  geom_tile(color = "black", linewidth = 0.2) +
  scale_fill_gradientn(
    colors = c("white", "#e31a1c", "#8b2500"),
    values = scales::rescale(c(0, 2, 100)),
    na.value = "gray90",
    limits = c(0, 100),
    oob = scales::squish) +
  guides (fill = guide_colourbar(barwidth = 0.5, barheight = 10,
                                 frame.colour = "black", frame.linewidth = 0.25,
                                 ticks.colour = NA)) + 
  labs(title = paste0("Inter_chromosomal discordant global distribution - ", strain_name, " - ", sample, " - ", blast_option),
       fill = "%") +
  scale_x_discrete(labels = c("ORF","Intergenic", "LTR","TEG", "Ty", "tRNA", "rRNA", "ncRNA", "snRNA", "snoRNA", "ARS", "Centromere", "Telomere")) +
  scale_y_discrete(limits = rev, labels = c("Telomere", "Centromere", "ARS", "snoRNA", "snRNA", "ncRNA", "rRNA", 
                                            "tRNA", "Ty", "TEG", "LTR", "Intergenic", "ORF")) +
  theme_minimal(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                               panel.background = element_blank(), 
                                               plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") + 
  theme(axis.text.y=element_text(size=8)) +
  theme(axis.text.x=element_text(size=8, angle = 90, hjust = 1)) +
  theme(axis.title.x = element_text(size=0)) +
  theme(axis.title.y = element_text(size=0)) + 
  theme(aspect.ratio = 1)


ggsave(
  filename = paste0(strain,"/", "Inter_chromosomal_discordant_pairs_global_distribution_heatmap_", strain_name, "_", sample, "_", blast_option, ".svg"),
  plot = heatmap,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)