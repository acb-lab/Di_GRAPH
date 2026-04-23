#### R script to processand plot SR coverage profiles from chr III and V, 75nt and 18nt, average all experiments
### Loop for each subdir/timepoint file in MYWD
### 21/04/2026 - Lydia



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
category <- args[3]
category_file <- args[4]
sample <- args[5]

strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))



# Load packages
library(tidyverse)

# Variables from Bash
cepa <- "${strain}"
category <- "${category}"
category_file <- "${category_file}"
subdir <- "${subdir_escaped}"
time <- "${time}"



root_dir <- "/Users/lab2.3/GWS3/WD_test_v0.3.2/1_Wt"
strain <- "1_Wt"
category <- "Cen"
category_file <- "/Users/lab2.3/Documents/Di_GRAPH/files/Categories/11.PMV.Cen.tsv"
sample <-"TSG"


# File paths
plot_path <- file.path(subdir, paste0(cepa, "_", category, "_T0vs", time,"_plot.svg"))
tsv_path <- file.path(subdir, paste0(category, "_fingerprint_", cepa, "_T0vs", time,".tsv"))
#all_data_path <- file.path(subdir, paste0(category, "_all_data_", cepa, "_T0vs", time, ".tsv")) only to check raw data

# Chromosome positions
chr_positions <- list(
  CHRI = 1:230218, CHRII = 1:813184, CHRIII = 1:316513,
  CHRIV = 1:1531933, CHRV = 1:578179, CHRVI = 1:270161,
  CHRVII = 1:1090940, CHRVIII = 1:562643, CHRIX = 1:439888,
  CHRX = 1:745751, CHRXI = 1:666816, CHRXII = 1:1078177,
  CHRXIII = 1:924431, CHRXIV = 1:784333, CHRXV = 1:1091291,
  CHRXVI = 1:948066
)

# Read category regions
df_categories <- read_tsv(category_file, col_names = FALSE,
                          col_types = cols(X1 = col_character(), X2 = col_character(),
                                           X3 = col_double(), X4 = col_double(), X5 = col_character())) %>%
  rename(Tipo = X1, Categoria = X2, Pos_inicio = X3, Pos_fin = X4, Cromosoma = X5) %>%
  mutate(Categoria = as.factor(Categoria))

# # Function to process each experiment
# process_experiment <- function(exp_num, time) {
#   data_file <- file.path(root_dir, paste0(time, "_", "E", exp_num, "_75nt.tsv"))
#   df <- read_tsv(data_file, col_names = FALSE, col_types = cols_only(X1 = col_character(), X4 = col_double()))
#   
#   df_full <- tibble(
#     Cepa = cepa,
#     Experimento = exp_num,
#     Tiempo = time,
#     Cromosoma = df$X1,
#     Valor_real = df$X4
#   )
#   
#   map_dfr(names(chr_positions), function(chr) {
#     chr_df <- df_full %>% filter(Cromosoma == chr) %>%
#       mutate(Posicion = chr_positions[[chr]])
#     
#     cat_df <- df_categorias %>% filter(Cromosoma == chr)
#     
#     map_dfr(unique(cat_df\$Categoria), function(cat) {
#       regions <- cat_df %>% filter(Categoria == cat)
#       map_dfr(1:nrow(regions), function(i) {
#         chr_df %>%
#           filter(Posicion >= regions\$Pos_inicio[i], Posicion <= regions\$Pos_fin[i]) %>%
#           mutate(Nombre = cat, categoria = regions\$Tipo[i])
#       })
#     })
#   })
# }

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
  
  return(processed_df)
}

# # Function to process each experiment
# process_experiment <- function(exp_num, time) {
#   data_file <- file.path(root_dir, paste0(time, "_", "E", exp_num, "_75nt.tsv"))
#   df <- read_tsv(data_file, col_names = FALSE, col_types = cols_only(X1 = col_character(), X4 = col_double()))
#   
#   df_full <- tibble(
#     Cepa = cepa,
#     Experimento = exp_num,
#     Tiempo = time,
#     Cromosoma = df$X1,
#     Valor_real = df$X4
#   )
#   
#   map_dfr(names(chr_positions), function(chr) {
#     chr_df <- df_full %>% filter(Cromosoma == chr) %>%
#       mutate(Posicion = chr_positions[[chr]])
#     
#     cat_df <- df_categorias %>% filter(Cromosoma == chr)
#     
#     map_dfr(unique(cat_df\$Categoria), function(cat) {
#       regions <- cat_df %>% filter(Categoria == cat)
#       map_dfr(1:nrow(regions), function(i) {
#         chr_df %>%
#           filter(Posicion >= regions\$Pos_inicio[i], Posicion <= regions\$Pos_fin[i]) %>%
#           mutate(Nombre = cat, categoria = regions\$Tipo[i])
#       })
#     })
#   })
# }


log_step("Finding _75nt coverage files...")
# Get all coverage.tsv files recursively in root folder
coverage_files <- list.files(
  path = root_dir,
  pattern = "_75nt\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)


df <-  purrr::map_dfr(coverage_files, process_coverage_files, categories_df = df_categories)


df %>% 
  group_by(strain, sample, chromosome, coordinate) %>% 
  summarise(avg_coverage = mean(coverage),
            sd_coverage = sd(coverage)) %>% 
  ungroup()


# Process all experiments
experiments <- expand.grid(exp = 1:3, time = c("T0", time))
all_data <- pmap_dfr(experiments, ~process_experiment(..1, ..2))

# Export all_data(Only to check raw data)
# write_tsv(all_data, all_data_path)
# message("Raw data saved to: ", all_data_path)

# Compute ratios
ratio_data <- all_data %>%
  group_by(Cepa, Cromosoma, categoria, Nombre, Experimento, Tiempo) %>%
  summarise(Valor_real = sum(Valor_real), .groups = "drop") %>%
  group_by(Cepa, Cromosoma, categoria, Nombre, Experimento) %>%
  mutate(ratio = Valor_real / Valor_real[Tiempo == "T0"]) %>%
  filter(Tiempo == time) %>%
  group_by(Cepa, Cromosoma, categoria, Nombre) %>%
  summarise(ratio_medio = mean(ratio), desv_est = sd(ratio), .groups = "drop")

# Plot
p <- ggplot(ratio_data) + 
  geom_col(aes(y = ratio_medio, x = Nombre), fill = "#2F2C7E", width = 0.7, 
           color = "black", linewidth = 0.1) +
  geom_errorbar(aes(x = Nombre, ymin = ratio_medio - desv_est, 
                    ymax = ratio_medio + desv_est), 
                width = 0.3, color = "grey8", alpha = 1, size = 0.3) +
  theme_classic() +
  labs(title = paste0("Ratio T0vs", time), subtitle = paste0("Cepa ", cepa), 
       y = "Ratio", x = "Feature")

ggsave(
  filename = plot_path,
  plot = p,
  width = 8,
  height = 6,
  dpi = 300,
  device = "svg"
)
write_tsv(ratio_data, file = tsv_path)