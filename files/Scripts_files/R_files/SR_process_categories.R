#### R script to process and plot SR coverage profiles from for genomic categories, 75nt, average all experiments
### Loop for each subdir/category file in MYWD
### 23/04/2026 - Lydia



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
  log_step("Usage: script.R <root_dir> <strain> <category> <category_path>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]
category <- args[3]
category_path <- args[4]


strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))
log_step(paste("CATEGORY:", category))



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







# Chromosome positions
chr_positions <- list(
  CHRI = 1:230218, CHRII = 1:813184, CHRIII = 1:316513,
  CHRIV = 1:1531933, CHRV = 1:578179, CHRVI = 1:270161,
  CHRVII = 1:1090940, CHRVIII = 1:562643, CHRIX = 1:439888,
  CHRX = 1:745751, CHRXI = 1:666816, CHRXII = 1:1078177,
  CHRXIII = 1:924431, CHRXIV = 1:784333, CHRXV = 1:1091291,
  CHRXVI = 1:948066
)

log_step("Finding category file...")
category_file <- list.files(
  path = category_path,
  pattern = paste0(category, ".*\\.tsv$"),
  recursive = TRUE,
  full.names = TRUE
)

# Read category regions
df_categories <- read_tsv(category_file, col_names = FALSE,
                          col_types = cols(X1 = col_character(), X2 = col_character(),
                                           X3 = col_double(), X4 = col_double(), X5 = col_character())) %>%
  rename(Tipo = X1, Categoria = X2, Pos_inicio = X3, Pos_fin = X4, Cromosoma = X5) %>%
  mutate(Categoria = as.factor(Categoria))



process_coverage_files <- function(file_path, categories_df) {
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
  coverage_file <- temp_file %>%
    rename("chromosome" = !!names(.[1]),
           "coverage" = !!names(.[4])) %>%
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name
    )
  
  processed_df <- map_dfr(names(chr_positions), function(chr) {
    chr_df <- coverage_file %>% filter(chromosome == chr) %>%
      mutate(coordinate = chr_positions[[chr]])
    
    cat_df <- categories_df %>% filter(Cromosoma == chr)
    
    map_dfr(unique(cat_df$Categoria), function(cat) {
      regions <- cat_df %>% filter(Categoria == cat)
      map_dfr(1:nrow(regions), function(i) {
        chr_df %>%
          filter(coordinate >= regions$Pos_inicio[i], coordinate <= regions$Pos_fin[i]) %>%
          mutate(Feature_name_A = cat, Category_A = regions$Tipo[i])
      })
    })
  })
  
  processed_df <- processed_df %>% 
    select(strain, sample, experiment, chromosome, coverage, Feature_name_A, Category_A)
  return(processed_df)
}



log_step("Finding _75nt coverage files...")
# Get all coverage.tsv files recursively in root folder
coverage_files <- list.files(
  path = root_dir,
  pattern = "_75nt\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)

log_step("Processing _75nt coverage files...")
coverage_processed_df <-  purrr::map_dfr(coverage_files, process_coverage_files, categories_df = df_categories)

log_step("Processing _75nt coverage files summary...")
coverage_processed_summary_df <- coverage_processed_df %>% 
  group_by(strain, sample, experiment, chromosome, Feature_name_A, Category_A) %>% 
  summarise(avg_coverage_feature = mean(coverage),
            sd_coverage_feature = sd(coverage)) %>% 
  ungroup() %>% 
  mutate(ratio_coverage_feature = avg_coverage_feature / avg_coverage_feature[sample == "T0"]) %>%
  ungroup() %>% 
  group_by(strain, sample, chromosome, Feature_name_A, Category_A) %>% 
  summarise(avg_ratio_coverage_feature = mean(ratio_coverage_feature),
            sd_ratio_coverage_feature = sd(ratio_coverage_feature)) %>% 
  ungroup() %>% 
  filter(sample != "T0")

## Save tsv file for each sample

for (s in unique(coverage_processed_summary_df$sample)) {
  df <- subset(coverage_processed_summary_df, sample == s)
  
  write_tsv(df, paste0(root_dir,"/", category, "_fingerprint_", strain_name,"_T0vs", s,".tsv"))
  
}
####
log_step("Plotting...")
## Save svg plot for each sample
for (s in unique(coverage_processed_summary_df$sample)) {
  df <- subset(coverage_processed_summary_df, sample == s)
  
  p <- ggplot(df) + 
    geom_col(aes(y = avg_ratio_coverage_feature, x = Feature_name_A), fill = "#2F2C7E", width = 0.7, 
             color = "black", linewidth = 0.1) +
    geom_errorbar(aes(x = Feature_name_A, ymin = avg_ratio_coverage_feature - sd_ratio_coverage_feature, 
                      ymax = avg_ratio_coverage_feature + sd_ratio_coverage_feature), 
                  width = 0.3, color = "grey8", alpha = 1, size = 0.3) +
    theme_classic() +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", color = NA)
    ) +
    labs(title = paste0("Ratio T0vs", s), subtitle = paste0("Strain ", strain_name), 
         y = "Ratio", x = "Feature")
  

  
  
  ggsave(
    filename = paste0(root_dir,"/", strain_name, "_", category, "_T0vs", s,"_plot.svg"),
    plot = p,
    width = 8,
    height = 6,
    dpi = 300,
    device = svglite,
    bg = "transparent"
  )
  
}
###

