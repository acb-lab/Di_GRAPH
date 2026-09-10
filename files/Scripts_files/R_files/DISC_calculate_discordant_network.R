#### R script to plot inter-chromosomal discordant network (blast level 1)
### Loop for each strain/sample in MYWD
### 15/04/2026 - Lydia

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
  log_step("Usage: script.R <root_dir> <strain> <sample>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]
sample   <- args[3]



strain <- sub("/$", "", strain)
strain_name <- basename(strain)


log_step(paste("ROOT_DIR:", root_dir))
log_step(paste("STRAIN:", strain))
log_step(paste("SAMPLE:", sample))

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



# Define a function to process tsv files
process_hotspots_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  # directory name above the file
  
  # Expect filename like T4_E07_MAT_filtered.tsv
  parts <- str_split(file_base, "_", simplify = TRUE)
  
  # Validate and extract parts safely
  if (ncol(parts) >= 3) {
    sample_name <- parts[1]
    option <- parts[2]
  } 
  else {
    warning(paste("Filename does not match expected format:", file_base))
    return(NULL)
  }
  
  # Read file
  temp_file <- read_tsv(file_path, col_names = TRUE, show_col_types = FALSE)
  
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  
  # Process file
  hotspots_file <- temp_file %>% 
    mutate(strain = strain_name, 
           sample = sample_name,
           blast_option = option) %>% 
    filter(global_freq > 6)
  
  return(hotspots_file)
}



# Define a function to process tsv files
process_matrix_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  # directory name above the file
  
  # Expect filename like T4_E07_MAT_filtered.tsv
  parts <- str_split(file_base, "_", simplify = TRUE)
  
  # Validate and extract parts safely
  if (ncol(parts) >= 3) {
    sample_name <- parts[1]
    option <- parts[2]
  } 
  else {
    warning(paste("Filename does not match expected format:", file_base))
    return(NULL)
  }
  
  # Read file
  temp_file <- read_tsv(file_path, col_names = TRUE, show_col_types = FALSE)
  
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  
  # Process file
  matrix_file <- temp_file %>% 
    mutate(strain = strain_name, 
           sample = sample_name,
           blast_option = option) %>% 
    select(Read_name_ID, 
           Chromosome_A, Feature_name_A, Category_A, Position_A, Essential_A,
           Chromosome_B, Feature_name_B, Category_B, Position_B, Essential_B,
           strain, sample, blast_option)
  
  return(matrix_file)
}


# Get all inter_discordant_pairs_unique_processed_valid_discordant_matrix_allexperiments.tsv files recursively in root folder
matrix_files <- list.files(
  path = root_dir,
  pattern = paste0(sample,".*option1_inter_discordant_pairs_unique_processed_valid_discordant_matrix_allexperiments\\.tsv$"),
  #pattern = "inter_discordant_pairs_unique_processed_valid_discordant_matrix_allexperiments\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing matrix files...")
matrix_processed_df <- purrr::map_dfr(matrix_files, process_matrix_files) 




log_step("Finding hotspots files...")

# Get all hotspots.tsv files recursively in root folder
hotspots_files <- list.files(
  path = root_dir,
  pattern = paste0(sample,".*option1_inter_discordant_pairs_unique_processed_valid_discordant_hotspots\\.tsv$"),
  #pattern = "inter_discordant_pairs_unique_processed_valid_discordant_hotspots\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing hotspots files...")
hotspots_processed_df <- purrr::map_dfr(hotspots_files, process_hotspots_files) 

write_tsv(hotspots_processed_df, file.path(root_dir, paste0(strain_name, "_", sample, "_hotspots_processed_df.tsv")))

hotspots_counts <- hotspots_processed_df %>%
  count(Feature_name_A, name = "n_occurrences")

#####

# Select hotspots from full matrix

hotspots_matrix <- left_join(hotspots_processed_df, matrix_processed_df, by = c("strain", "sample", "blast_option", "Feature_name_A")) %>% 
  select(Feature_name_A, Category_A.x, Chromosome_A.x, Essential_A.x, Position_A.x,
         Feature_name_B, Category_B, Chromosome_B, Essential_B, Position_B,
         strain, sample, blast_option, global_freq) %>% 
  rename("Category_A" = !!names(.[2]), "Chromosome_A" = !!names(.[3]),
         "Essential_A" = !!names(.[4]), "Position_A" = !!names(.[5])) %>% 
  unique() %>% 
  mutate(id = 1:n()) 

write_tsv(hotspots_matrix, file.path(root_dir, paste0(strain_name, "_", sample, "_hotspots_matrix.tsv")))



hotspots_matrix_A <- hotspots_matrix %>% 
  select(Feature_name_A, Category_A, strain, sample, id)

hotspots_matrix_B <- hotspots_matrix %>% 
  select(Feature_name_B, Category_B, strain, sample, id) %>% 
  rename(Category_A = Category_B)

hotspots_matrix_size <- hotspots_matrix %>% 
  group_by(Feature_name_A) %>% 
  count(Feature_name_A, name = "n_occurrences") %>% 
  ungroup()


hotspots_matrix_size_B <- hotspots_matrix %>% 
  group_by(Feature_name_B) %>% 
  count(Feature_name_B, name = "n_occurrences") %>% 
  ungroup()

hotspots_matrix_size_B_no_A <- hotspots_matrix_size_B %>% 
  rename(Feature_name_A = Feature_name_B) %>% 
  anti_join(., hotspots_matrix_size , by = "Feature_name_A") %>% 
  rename(Feature_name_B = Feature_name_A)

hotspots_matrix_A <- left_join(hotspots_matrix_A, hotspots_matrix_size) 

hotspots_matrix_B <- left_join(hotspots_matrix_B, hotspots_matrix_size_B_no_A) %>% 
  rename(Feature_name_A = Feature_name_B)

hotspots_matrix_A_B <- bind_rows(hotspots_matrix_A, hotspots_matrix_B)


edges <- hotspots_matrix %>%  select(Feature_name_A, Feature_name_B, Category_A, strain, id) %>% 
  rename("Source" = !!names(.[1]), "Target" = !!names(.[2])) %>% 
  as.data.frame()


########## NETWORK

library(tidyverse)
library(igraph)
library(ggraph)


node_sizes <- hotspots_matrix_size %>% 
  rename("name" = !!names(.[1]), "Size" = !!names(.[2]))


g <- graph_from_data_frame(edges, directed = FALSE)


# Derive node data
nodes <- data.frame(name = V(g)$name) %>%
  mutate(order = 1:n()) %>%
  left_join(hotspots_matrix_A_B %>% select(Feature_name_A, Category_A, strain, id, n_occurrences),
            by = c("name" = "Feature_name_A")) %>% 
  arrange(order)

nodes_def <- nodes %>% group_by(name) %>%
  slice(1) %>%
  ungroup() %>% 
  arrange(order)





# Attach node attributes to igraph
V(g)$Category_A <- nodes_def$Category_A
V(g)$Size <- nodes_def$n_occurrences
V(g)$strain <- nodes_def$strain
#V(g)$name

custom_palette <- c(
  "1_Wt" = "#00AFBB",   
  "2_exo1d" = "#293781",  
  "3_sgs1d" = "#AE2D2C",   
  "4_srs2d" = "#439645", 
  "5_rad51d" = "#5D297D",
  "1_Wt_2_exo1d" = "#15739e",
  "1_Wt_3_sgs1d" = "#576e74",
  "1_Wt_4_srs2d" = "#22a380",
  "1_Wt_5_rad51d" = "#2f6c9c",
  "2_exo1d_3_sgs1d " = "#6c3257",
  "2_exo1d_4_srs2d" = "#366763",
  "2_exo1d_5_rad51d" = "#43307f",
  "3_sgs1d_4_srs2d" = "#796239",
  "3_sgs1d_5_rad51d" = "#862b55",
  "4_srs2d_5_rad51d" = "#506061",
  "1_Wt_2_exo1d_3_sgs1d" = "#485c78",
  "1_Wt_2_exo1d_4_srs2d" = "#247f80",
  "1_Wt_2_exo1d_5_rad51d" = "#2d5a93",
  "1_Wt_3_sgs1d_4_srs2d" = "#507b64",
  "1_Wt_3_sgs1d_5_rad51d" = "#595777",
  "1_Wt_4_srs2d_5_rad51d" = "#357a7f",
  "2_exo1d_3_sgs1d_4_srs2d" = "#5e5351",
  "2_exo1d_3_sgs1d_5_rad51d" = "#672f63",
  "2_exo1d_4_srs2d_5_rad51d" = "#43526c",
  "3_sgs1d_4_srs2d_5_rad51d" = "#6f4f4f",
  "1_Wt_2_exo1d_3_sgs1d_4_srs2d" = "#476a6b",
  "1_Wt_2_exo1d_3_sgs1d_5_rad51d" = "#4d4f79",
  "1_Wt_2_exo1d_4_srs2d_5_rad51d" = "#326980",
  "1_Wt_3_sgs1d_4_srs2d_5_rad51d" = "#54676a",
  "2_exo1d_3_sgs1d_4_srs2d_5_rad51d" = "#5e495c",
  "1_Wt_2_exo1d_3_sgs1d_4_srs2d_5_rad51d" = "#4b5d6f"
    
)

custom_palette_categories = c("ORF" = "#4F71BE",
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
                              "telomere" = "#4B6733")

svglite::svglite(file = paste0(root_dir, "/", strain_name, "_", sample, "_Network_plot", ".svg"), width = 13, height = 13)
ggraph(g, layout = "fr") +  # "fr" = Fruchterman-Reingold layout
  geom_edge_link(alpha = 0.2, color = "black") +
  geom_node_point(aes(size = Size, color = Category_A), alpha = 1) +
  scale_color_manual(values = custom_palette_categories) +
  #geom_node_text(aes(label = name), repel = TRUE, size = 2) +
  geom_node_text(
    aes(label = ifelse(Size > 10 | name %in% hotspots_counts$Feature_name_A, name, "")),
    repel = TRUE,
    size = 3
  ) + 
  scale_size_continuous(range = c(3, 15)) +
  theme_void() +
  theme(legend.position = "none") +
  ggtitle(paste0(strain_name, " - ", sample, " - Network"))

dev.off()



hotspots_matrix_summary <- hotspots_matrix %>% 
  group_by(sample) %>% 
  mutate(total_sample = n()) %>% 
  ungroup() %>% 
  group_by(sample, strain) %>% 
  mutate(total_sample_strain = n()) %>% 
  ungroup() %>% 
  group_by(sample, Category_A) %>% 
  mutate(total_sample_category_A = n()) %>%
  ungroup() %>% 
  group_by(sample, Category_A, strain) %>% 
  mutate(total_sample_category_A_strain = n()) %>%
  ungroup() %>% 
  group_by(sample, Category_A, Category_B) %>% 
  mutate(total_sample_category_A_Category_B = n()) %>%
  ungroup() %>% 
  group_by(sample, Category_A, Category_B, strain) %>% 
  mutate(total_sample_category_A_Category_B_strain = n()) %>%
  ungroup() %>% 
  mutate(freq_category_A_100  = 100*(total_sample_category_A_Category_B_strain/total_sample_category_A),
         freq_sample_100  = 100*(total_sample_category_A_Category_B_strain/total_sample),
         freq_sample_100_soloA  = 100*(total_sample_category_A_strain/total_sample)) %>% 
  select(Category_A, Category_B,
         total_sample,
         total_sample_strain, 
         total_sample_category_A,
         total_sample_category_A_strain,
         total_sample_category_A_Category_B,
         total_sample_category_A_Category_B_strain,
         freq_category_A_100, freq_sample_100,
         freq_sample_100_soloA,
         strain, sample) %>% 
  unique()


hotspots_matrix_summary_positions_strain_sep <- hotspots_matrix_summary %>% 
  mutate(position_x = ifelse(Category_B == "ORF", 1,
                             ifelse(Category_B == "intergenic", 2,
                                    ifelse(Category_B == "long_terminal_repeat", 3,
                                           ifelse(Category_B == "transposable_element_gene", 4,
                                                  ifelse(Category_B == "LTR_retrotransposon", 5,
                                                         ifelse(Category_B == "tRNA_gene", 6,
                                                                ifelse(Category_B == "rRNA_gene", 7,
                                                                       ifelse(Category_B == "ncRNA_gene", 8,
                                                                              ifelse(Category_B == "snRNA_gene", 9,
                                                                                     ifelse(Category_B == "snoRNA_gene", 10,
                                                                                            ifelse(Category_B == "ARS", 11,
                                                                                                   ifelse(Category_B == "centromere", 12,
                                                                                                          ifelse(Category_B == "telomere", 13,"_")))))))))))))) %>% 
  mutate(position_y = ifelse(Category_A == "ORF", 13,
                             ifelse(Category_A == "intergenic", 12,
                                    ifelse(Category_A == "long_terminal_repeat", 11,
                                           ifelse(Category_A == "transposable_element_gene", 10,
                                                  ifelse(Category_A == "LTR_retrotransposon", 9,
                                                         ifelse(Category_A == "tRNA_gene", 8,
                                                                ifelse(Category_A == "rRNA_gene", 7,
                                                                       ifelse(Category_A == "ncRNA_gene", 6,
                                                                              ifelse(Category_A == "snRNA_gene", 5,
                                                                                     ifelse(Category_A == "snoRNA_gene", 4,
                                                                                            ifelse(Category_A == "ARS", 3,
                                                                                                   ifelse(Category_A == "centromere", 2,
                                                                                                          ifelse(Category_A == "telomere", 1,"_"))))))))))))))


hotspots_matrix_summary_positions_strain_sep$position_x <- as.numeric(hotspots_matrix_summary_positions_strain_sep$position_x)
hotspots_matrix_summary_positions_strain_sep$position_y <- as.numeric(hotspots_matrix_summary_positions_strain_sep$position_y)

write_tsv(hotspots_matrix_summary_positions_strain_sep, file.path(root_dir, paste0(strain_name, "_", sample, "_Network_summary_positions_strain_sep.tsv")))

hotspots_matrix_summary_positions_strain_sep_complete <- hotspots_matrix_summary_positions_strain_sep %>%
  complete(
    strain,
    sample,
    position_x = 1:13,
    position_y = 1:13,
    fill = list(freq_sample_100 = NA)
  )

p <- ggplot(hotspots_matrix_summary_positions_strain_sep_complete, aes(y = position_y, x = position_x)) +
  geom_tile(aes(fill = freq_sample_100), size = 1)+
  scale_fill_gradient2(low = "white",
                       mid = "#e31a1c",
                       high="blue2",
                       midpoint = 25,
                       na.value = "white",
                       limits = c(0, 50),
                       oob = scales::squish) +
  guides (fill = guide_colourbar(barwidth = 0.5, barheight = 10,
                                 frame.colour = "black", frame.linewidth = 0.25,
                                 ticks.colour = NA)) + 
  labs(fill = "%") +
  theme_classic(base_family = "Helvetica") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                   panel.background = element_blank(), 
                                                   plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +
  coord_cartesian(expand = FALSE) +
  geom_hline(yintercept = 0.5 + 0:13, colour = "gray7", size = 0.15) +
  geom_vline(xintercept = 0.5 + 0:13, colour = "gray7", size = 0.15) +
  scale_x_continuous(breaks = seq(1,13,1),
                     labels = c("ORF","Intergenic", "LTR","TEG", "Ty", "tRNA", "rRNA", 
                                "ncRNA", "snRNA", "snoRNA", "ARS", "Centromere", "Telomere")) +
  scale_y_continuous(breaks = seq(1,13,1),
                     labels = c("Telomere", "Centromere", "ARS", "snoRNA", "snRNA", "ncRNA", "rRNA", 
                                "tRNA", "Ty", "TEG", "LTR", "Intergenic", "ORF")) +
  theme(aspect.ratio = 1) +
  theme(axis.text.x=element_text(size=8, angle = 90, hjust = 1)) +
  theme(axis.title.x = element_text(size = 0)) +
  theme(axis.title.y = element_text(size = 0)) +
  ggtitle("Network quantification") +
  theme(axis.ticks.x = element_line()) +
  facet_grid(strain ~ sample) 

ggsave(
  filename = paste0(root_dir, "/", strain_name, "_", sample, "_Network_quantification_plot", ".svg"),
  plot = p,
  #width = 8,
  #height = 3,
  device = svglite,
  bg = "transparent"
)