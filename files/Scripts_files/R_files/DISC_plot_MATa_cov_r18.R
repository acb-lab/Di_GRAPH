#### R script to plot MATa-MATa' coverage (reads 18nt)
### Loop for each strain in MYWD
### 11/04/2026 - Lydia

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
strain   <- args[2]


strain <- sub("/$", "", strain)
strain_name <- basename(strain)

log_step(paste("ROOT_DIR:", root_dir))
log_step(paste("STRAIN:", strain))


# load libraries

library(ggplot2)
library(extrafont)
library(svglite)
library(ggdensity)
library(purrr)
library(stringr)
library(readr)
library(ggforce)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)


# Define a function to process tsv files
process_MATs_coverage_files <- function(file_path) {
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
  
  # Prepare coverage dataframe
  coverage_file <- temp_file %>%
    rename("chromosome" = !!names(.[1]), "coverage" = !!names(.[4])) %>%
    select(chromosome, coverage) %>%
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name
    )
  
  return(coverage_file)
}


log_step("Finding MATa-MATa' coverage files...")
# Get all coverage.tsv files recursively in root folder
MAT_coverage_files <- list.files(
  path = root_dir,
  pattern = "MAT_complete_r18_filtered\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)
MAT_coverage_files
log_step("Processing MATa-MATa' coverage files...")
MATs_processed_df <- purrr::map_dfr(MAT_coverage_files, process_MATs_coverage_files)


log_step("Processing chromosome III...")
III_df <- MATs_processed_df %>% filter(chromosome == "CHRIII") %>%
  group_by(strain, sample, experiment) %>% 
  mutate(position = row_number(), position_HO_centered = position - 200753) %>%
  ungroup() %>% 
  group_by(strain, sample, chromosome, position, position_HO_centered) %>%
  summarise(mean_coverage = mean(coverage, na.rm = TRUE), sd_coverage = sd(coverage, na.rm = TRUE))

log_step("Processing chromosome V...")
V_df <- MATs_processed_df %>% filter(chromosome == "CHRV") %>%
  group_by(strain, sample, experiment) %>% 
  mutate(position = row_number(), position_HO_centered = position - 289825) %>%
  ungroup() %>% 
  group_by(strain, sample, chromosome, position, position_HO_centered) %>%
  summarise(mean_coverage = mean(coverage, na.rm = TRUE), sd_coverage = sd(coverage, na.rm = TRUE))

MAT_df <- bind_rows(III_df, V_df)
log_step(paste("Number of samples found:", length(MAT_df)))


# Set polymorphisms positions
polymorphisms <- tibble(position_HO_centered = c(-634,-586,-541,-481,-427,-367,-304,-244,-211,-178,-118,-64,
                                                 0,64,129,194,259,324,395,454,519,584,649),
                        mean_coverage = rep(c(0)),
                        chromosome = rep(c("polymorphism")))

log_step("Splitting by sample...")
# Split by sample
split_MAT_df <- split(MAT_df, 
                      MAT_df$sample)

# Set output directory
#svg_output_dir <- file.path(root_dir, "MAT_plots_r18")
#dir.create(svg_output_dir, showWarnings = FALSE, recursive = TRUE)
write_tsv(MAT_df, file.path(root_dir, paste0(strain_name, "_MAT_coverage_r18_df.tsv")))


log_step("Plotting...")
# Generate plots
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  
  p <- ggplot(data_subset, aes(x = position_HO_centered, y = mean_coverage, color = chromosome)) +
    geom_point(data = polymorphisms, aes(x = position_HO_centered, y = mean_coverage),
               color = "chartreuse3", shape = 108, size = 5) +
    geom_line(linewidth = 1) + #1.2
    scale_color_manual(values = c("CHRIII" = "#2F2C7E", "CHRV" = "#A30000")) + #2F2C7E, A30000
    theme_classic(base_family = "Arial") +
    theme(
      panel.grid = element_line(color = "black", linewidth = 0.1),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none",
      aspect.ratio = 0.45,
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(hjust = 1)
    ) +
    theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 20),
          axis.title.y = element_text(vjust = 2, size = 20)) +
    theme(axis.text.x = element_text(vjust = 0, size = 15),
          axis.text.y = element_text(vjust = 0, size = 15)) +
    coord_cartesian(xlim = c(-1200, 1200), ylim = c(0, 30), expand = FALSE) +
    scale_x_continuous(limits = c(-1100, 1100), breaks = seq(-1100, 1100, 200)) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey3", linewidth = 0.5) +
    labs(
      title = paste0("Average coverage of MATa-MATa' discordant reads r18 - ", strain_name, " - ", sample),
      x = "Relative Coordinate (bp)", #Relative Coordinate (bp)",
      y = "Average Coverage" #Average Coverage"
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("Average_MATa_MATa_coverage_r18_", strain_name, "_", sample, ".svg")),
    plot = p,
    #width = 8,
    #height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})

# Generate plots for chromosome III (only coverage shape)
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  data_subset <- data_subset %>% filter(chromosome == "CHRIII")
  
  p <- ggplot(data_subset, aes(x = position_HO_centered, y = mean_coverage, color = chromosome)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("CHRIII" = "#2F2C7E")) +
    theme_classic(base_family = "Arial") +
    theme(
      panel.grid = element_line(color = "black", linewidth = 0.1),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none",
      aspect.ratio = 0.25,
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(hjust = 1)
    ) +
    theme(axis.line=element_blank(), axis.ticks=element_blank()) +
    theme(axis.text.x = element_text(vjust = 0, size = 0),
          axis.text.y = element_text(vjust = 0, size = 0))+
    coord_cartesian(xlim = c(-1200, 1200), ylim = c(0, 30), expand = FALSE) +
    scale_x_continuous(limits = c(-1100, 1100), breaks = seq(-1100, 1100, 200)) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
    labs(
      title = paste0("Average coverage of MATa-MATa' discordant reads chrIII r18 - ", strain_name, " - ", sample),
      x = "",
      y = ""
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("Average_MATa_MATa_coverage_r18_chrIII_", strain_name, "_", sample, ".svg")),
    plot = p,
    width = 8,
    height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})

# Generate plots for chromosome V (only coverage shape)
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  data_subset <- data_subset %>% filter(chromosome == "CHRV")
  
  p <- ggplot(data_subset, aes(x = position_HO_centered, y = mean_coverage, color = chromosome)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("CHRV" = "#A30000")) +
    theme_classic(base_family = "Arial") +
    theme(
      panel.grid = element_line(color = "black", linewidth = 0.1),
      panel.background = element_blank(),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none",
      aspect.ratio = 0.25,
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(hjust = 1)
    ) +
    theme(axis.line=element_blank(), axis.ticks=element_blank()) +
    theme(axis.text.x = element_text(vjust = 0, size = 0),
          axis.text.y = element_text(vjust = 0, size = 0))+
    coord_cartesian(xlim = c(1200, -1200), ylim = c(0, 30), expand = FALSE) +
    scale_x_continuous(limits = c(-1100, 1100), breaks = seq(-1100, 1100, 200)) +
    scale_y_continuous(limits = c(0, 30), breaks = seq(0, 30, 5)) +
    labs(
      title = paste0("Average coverage of MATa-MATa' discordant reads chrV r18 - ", strain_name, " - ", sample),
      x = "",
      y = ""
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("Average_MATa_MATa_coverage_r18chrV_", strain_name, "_", sample, ".svg")),
    plot = p,
    width = 8,
    height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})

# Define a function to process sam files
process_MATs_pairs_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  # directory name above the file
  
  
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
  temp_file <- read_tsv(file_path, skip = 19, col_names = FALSE, show_col_types = FALSE)
  
  # Skip if empty
  if (nrow(temp_file) == 0) {
    return(NULL)
  }
  
  # Prepare coverage dataframe
  pairs_file <- temp_file %>%
    rename("chromosome_A" = !!names(.[3]), "position_A" = !!names(.[4]),
           "chromosome_B" = !!names(.[7]), "position_B" = !!names(.[8])) %>%
    select(chromosome_A, position_A, chromosome_B, position_B) %>%
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name
    )
  
  return(pairs_file)
}

log_step("Finding MATa-MATa' pairs files...")
# Get all pairs files recursively in root folder
MAT_pairs_files <- list.files(
  path = root_dir,
  pattern = "unique_MAT_complete_r18\\.sam$",
  recursive = TRUE,
  full.names = TRUE
)
MAT_pairs_files
log_step("Processing MATa-MATa' pairs files...")
MATs_pairs_processed_df <- purrr::map_dfr(MAT_pairs_files, process_MATs_pairs_files)

MAT_pairs_III <- MATs_pairs_processed_df %>% filter(chromosome_A == "CHRIII" & chromosome_B == "CHRV") %>% 
  rename(X_value_III = position_A, Y_value_V = position_B) %>% 
  select(X_value_III, Y_value_V, strain, sample, experiment)

MAT_pairs_V <- MATs_pairs_processed_df %>% filter(chromosome_A == "CHRV" & chromosome_B == "CHRIII") %>% 
  rename(X_value_III = position_B, Y_value_V = position_A) %>% 
  select(X_value_III, Y_value_V, strain, sample, experiment)

MAT_combined <- bind_rows(MAT_pairs_III, MAT_pairs_V)
log_step("Saving MATa-MATa' pairs files...")
write_tsv(MAT_combined, file.path(root_dir, paste0(strain_name, "_MAT_pairs_r18_df.tsv")))

# ### Calculate proportion of GC and CO events

log_step("Calculating proportion of GC and CO events...")

poly_III <- tibble(Polymorphism = c("P_-12","P_-11","P_-10","P_-9","P_-8","P_-7","P_-6",
                                    "P_-5","P_-4","P_-3","P_-2","P_-1","P_0","P_1","P_2",
                                    "P_3","P_4","P_5","P_6","P_7","P_8","P_9","P_10"),
                   Start = c(200102,200150,200195,200255,200309,200369,200432,
                             200492,200525,200558,200618,200672,200736,200800,200865,
                             200930,200995,201060,201131,201190,201255,201320,201385),
                   End = c(200136,200184,200229,200289,200343,200403,200466,
                           200526,200559,200592,200652,200706,200770,200834,200899,
                           200964,201029,201094,201165,201224,201289,201354,201419))

poly_V <- tibble(Polymorphism = c("P_-12","P_-11","P_-10","P_-9","P_-8","P_-7","P_-6",
                                  "P_-5","P_-4","P_-3","P_-2","P_-1","P_0","P_1","P_2",
                                  "P_3","P_4","P_5","P_6","P_7","P_8","P_9","P_10"),
                 Start = c(289174,289222,289267,289327,289381,289441,289504,
                           289564,289597,289630,289690,289744,289808,289872,289937,
                           290002,290067,290132,290203,290262,290327,290392,290457),
                 End = c(289208,289256,289301,289361,289415,289475,289538,
                         289598,289631,289664,289724,289778,289842,289906,289971,
                         290036,290101,290166,290237,290296,290361,290426,290491))


process_MATs_18 <- function(MAT_combined_df, poly_III_df, poly_V_df){
  breaks_V <- c(poly_V_df$Start, tail(poly_V_df$End, 1))
  breaks_III <- c(poly_III_df$Start, tail(poly_III_df$End, 1))
  MAT_combined_df$Polymorphism_V <- cut(MAT_combined_df$Y_value_V, 
                                        breaks = breaks_V, 
                                        labels = poly_V_df$Polymorphism,
                                        include.lowest = FALSE,
                                        right = FALSE)
  
  MAT_combined_df$Polymorphism_III <- cut(MAT_combined_df$X_value_III, 
                                          breaks = breaks_III, 
                                          labels = poly_III_df$Polymorphism,
                                          include.lowest = FALSE,
                                          right = FALSE)
  
  MAT_combined <- MAT_combined_df %>% mutate(III_HO_centered = X_value_III - 200753,
                                             V_HO_centered = Y_value_V -289825,
                                             abs = abs(III_HO_centered) - abs(V_HO_centered),
                                             left_right = ifelse(X_value_III <= 200753, "left", "right"),
                                             group = ifelse(left_right == "left" & abs > 0, 1,
                                                            ifelse(left_right == "left" & abs < 0, 4,
                                                                   ifelse(left_right == "right" & abs > 0, 3,
                                                                          ifelse(left_right == "right" & abs < 0, 2, NA_real_)))),
                                             GC_CO = ifelse(group == 1 | group == 3, "GC", "CO"),
                                             group_GC = ifelse(group == 1, "GC_left",
                                                               ifelse(group == 3 & X_value_III < 200834, "GC_left",
                                                                      ifelse(group ==3 & X_value_III > 200834, "GC_right", 
                                                                             ifelse(group == 2, "d_GC_right",
                                                                                    ifelse(group == 4, "d_GC_left", NA_real_))))))
  
  MAT_combined_A <- MAT_combined %>% mutate(start = ifelse(group == 1, -1.5708,
                                                           ifelse(group == 2, -1.5708,
                                                                  ifelse(group == 3, 1.5708,
                                                                         ifelse(group == 4, 1.5708, NA_real_))))) %>% 
    mutate(AB = rep(c("A")))
  
  MAT_combined_B <- MAT_combined %>% mutate(start = ifelse(group == 1, 1.5708,
                                                           ifelse(group == 2, 1.5708,
                                                                  ifelse(group == 3, -1.5708,
                                                                         ifelse(group == 4, -1.5708, NA_real_))))) %>% 
    mutate(AB = rep(c("B")))
  
  MAT_combined_AB <- bind_rows(MAT_combined_A, MAT_combined_B)
  
  return(MAT_combined_AB)
  
}


MAT_combined_processed <- process_MATs_18(MAT_combined, poly_III, poly_V)



write_tsv(MAT_combined_processed, file.path(root_dir, paste0(strain_name, "_MAT_pairs_r18_processed_df.tsv")))

log_step("Calculating GC and CO distribution...")

GC_CO_summary <- MAT_combined_processed %>% filter(AB == "A") %>% 
  select(strain, sample, experiment, GC_CO) %>% 
  group_by(strain, sample, experiment, GC_CO) %>% 
  summarise(n_experiment = n()) %>%
  mutate(total_experiment = sum(n_experiment)) %>%                  
  mutate(percentage_experiment = 100 * n_experiment / total_experiment) %>%         
  ungroup() %>% select(!c(n_experiment, total_experiment)) 

# Define the fixed order of groups
GC_CO_levels <- c("GC", "CO")

GC_CO_combined_summary <- GC_CO_summary %>% 
  complete(strain, sample, experiment, GC_CO = GC_CO_levels, fill = list(percentage_experiment = 0)) %>% 
  group_by(strain, sample, GC_CO) %>% 
  summarise(average_percentage = mean(percentage_experiment, na.rm = TRUE),
            sd_percentage = sd(percentage_experiment, na.rm = TRUE)) %>%
  ungroup()

write_tsv(GC_CO_combined_summary, file.path(root_dir, paste0(strain_name, "_MAT_pairs_r18_GC_CO_summary_df.tsv")))

# Define the fixed order of groups
GC_CO_levels <- c("GC", "CO")

# Set the GC_CO factor with all levels
GC_CO_combined_summary <-GC_CO_combined_summary %>%
  mutate(GC_CO = factor(GC_CO, levels = GC_CO_levels))

GC_CO_combined_summary <-GC_CO_combined_summary %>%
  complete(sample, GC_CO = GC_CO_levels,
           fill = list(average_percentage = 0, sd_percentage = 0))


log_step("Calculating GC groups distribution...")


GC_groups_summary <- MAT_combined_processed %>% filter(AB == "A") %>% 
  select(strain, sample, experiment, group_GC) %>% 
  group_by(strain, sample, experiment, group_GC) %>% 
  summarise(n_experiment = n()) %>%
  mutate(total_experiment = sum(n_experiment)) %>%                  
  mutate(percentage_experiment = 100 * n_experiment / total_experiment) %>%         
  ungroup() %>% select(!c(n_experiment, total_experiment)) 

all_groups <- c("GC_left", "GC_right", "d_GC_left", "d_GC_right")
GC_groups_combined_summary <- GC_groups_summary %>% 
  complete(strain, sample, experiment, group_GC = all_groups, fill = list(percentage_experiment = 0)) %>% 
  group_by(strain, sample, group_GC) %>%
  summarise(average_percentage = mean(percentage_experiment, na.rm = TRUE),
            sd_percentage = sd(percentage_experiment, na.rm = TRUE)) %>%
  ungroup()

write_tsv(GC_groups_combined_summary, file.path(root_dir, paste0(strain_name, "_MAT_pairs_r18_groups_summary_df.tsv")))

# Define the fixed order of groups
GC_groups_levels <- c("GC_left", "GC_right", "d_GC_left", "d_GC_right")

# Set the group factor with all levels
GC_groups_combined_summary <-GC_groups_combined_summary %>%
  mutate(GC_group = factor(group_GC, levels = GC_groups_levels))

GC_groups_combined_summary <-GC_groups_combined_summary %>%
  complete(sample, GC_group = GC_groups_levels,
           fill = list(average_percentage = 0, sd_percentage = 0))


# Calculate proportion of polymorphism insertion for each polymorphism (III and V) - total-reads

MAT_r18_count_total_reads <- MAT_combined_processed %>% filter(AB == "A") %>% 
  select(X_value_III, Y_value_V, strain, sample, experiment, Polymorphism_V, Polymorphism_III) %>% 
  pivot_longer(cols= !c(X_value_III, Y_value_V, strain, sample, experiment)) %>% 
  rename(chromosome = name, polymorphism = value) %>% 
  group_by(strain, sample, experiment, chromosome, polymorphism) %>% 
  summarise(n_experiment = n()) %>%
  mutate(total_experiment = sum(n_experiment)) %>%                  
  mutate(percentage_experiment = 100 * n_experiment / total_experiment) %>%         
  ungroup() %>% select(!c(n_experiment, total_experiment)) 

# Define the fixed order of polymorphisms
polymorphism_levels <- c("P_-12","P_-11","P_-10","P_-9","P_-8","P_-7","P_-6",
                         "P_-5","P_-4","P_-3","P_-2","P_-1","P_0","P_1","P_2",
                         "P_3","P_4","P_5","P_6","P_7","P_8","P_9","P_10")

MAT_r18_count_total_reads_summary <- MAT_r18_count_total_reads %>% 
  complete(strain, sample, chromosome, experiment, polymorphism = polymorphism_levels, fill = list(percentage_experiment = 0)) %>%
  group_by(strain, sample, chromosome, polymorphism) %>% 
  summarise(average_percentage = mean(percentage_experiment, na.rm = TRUE),
            sd_percentage = sd(percentage_experiment, na.rm = TRUE)) %>% 
  ungroup()

#   # Define the fixed order of polymorphisms
# polymorphism_levels <- c("P_-12","P_-11","P_-10","P_-9","P_-8","P_-7","P_-6",
#                         "P_-5","P_-4","P_-3","P_-2","P_-1","P_0","P_1","P_2",
#                         "P_3","P_4","P_5","P_6","P_7","P_8","P_9","P_10")

# Set the polymorphism factor with all levels
MAT_r18_count_total_reads_summary <- MAT_r18_count_total_reads_summary %>%
  mutate(polymorphism = factor(polymorphism, levels = polymorphism_levels)) 

MAT_r18_count_total_reads_summary <- MAT_r18_count_total_reads_summary %>%
  complete(sample, polymorphism = polymorphism_levels, chromosome,
           fill = list(average_percentage = 0, sd_percentage = 0))



write_tsv(MAT_r18_count_total_reads_summary, file.path(root_dir, paste0(strain_name, "_MAT_pairs_r18_count_total_reads_summary_df.tsv")))

# Calculate proportion of polymorphism insertion for each polymorphism (III and V) - HOinc-reads

MAT_r18_count_HOinc_reads <- MAT_combined_processed %>% filter(AB == "A") %>% 
  select(X_value_III, Y_value_V, strain, sample, experiment, Polymorphism_V, Polymorphism_III) %>% 
  pivot_longer(cols= !c(X_value_III, Y_value_V, strain, sample, experiment)) %>% 
  rename(chromosome = name, polymorphism = value) %>% 
  group_by(strain, sample, experiment, chromosome, polymorphism) %>% 
  summarise(n_HOinc = n(), .groups = "drop") %>%
  mutate(n_P0_V = if_else(polymorphism == "P_0", n_HOinc, 0L)) %>% 
  group_by(strain, sample, experiment) %>% 
  mutate(n_P0_V_total = sum(n_P0_V, na.rm = TRUE)) %>%
  rename(count = n_HOinc) %>% 
  select(!c(n_P0_V)) %>% 
  mutate(percentage_experiment = 100 * count / n_P0_V_total) %>%         
  ungroup() %>% select(!c(count, n_P0_V_total)) 


MAT_r18_count_HOinc_reads_summary <- MAT_r18_count_HOinc_reads %>% 
  complete(strain, sample, chromosome, experiment, polymorphism = polymorphism_levels, fill = list(percentage_experiment = 0)) %>%
  group_by(strain, sample, chromosome, polymorphism) %>% 
  summarise(average_percentage = mean(percentage_experiment, na.rm = TRUE),
            sd_percentage = sd(percentage_experiment, na.rm = TRUE)) %>% 
  ungroup()

# Define the fixed order of polymorphisms
polymorphism_levels <- c("P_-12","P_-11","P_-10","P_-9","P_-8","P_-7","P_-6",
                         "P_-5","P_-4","P_-3","P_-2","P_-1","P_0","P_1","P_2",
                         "P_3","P_4","P_5","P_6","P_7","P_8","P_9","P_10")

# Set the polymorphism factor with all levels
MAT_r18_count_HOinc_reads_summary <- MAT_r18_count_HOinc_reads_summary %>%
  mutate(polymorphism = factor(polymorphism, levels = polymorphism_levels))

MAT_r18_count_HOinc_reads_summary <- MAT_r18_count_HOinc_reads_summary %>%
  complete(sample, polymorphism = polymorphism_levels, chromosome,
           fill = list(average_percentage = 0, sd_percentage = 0))

write_tsv(MAT_r18_count_HOinc_reads_summary, file.path(root_dir, paste0(strain_name, "_MAT_pairs_r18_count_HOinc_reads_summary_df.tsv")))



# GC_CO plotting

split_MAT_df <- split(MAT_combined_processed, 
                      MAT_combined_processed$sample)


log_step("Plotting...")
# Generate plots for MATa-MATa' pairs
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  
  p <- ggplot(data_subset) + theme_light(base_family = "Arial") +
    geom_arc_bar(aes(x0 = X_value_III, y0 = Y_value_V, r0 = 0, r = 35,
                     start = (start)+(pi/2), end = (start)+(pi/2) + pi, fill = AB),
                 color = NA, alpha = 0.15) +
    scale_fill_manual(values = c("A" = "#A30000", "B" = "#2F2C7E")) +
    #scale_fill_manual(values = c("A" = "#B86029", "B" = "#081D5C")) +
    coord_cartesian(xlim = c(199653,201853), ylim = c(288725, 290925), expand=FALSE) +
    scale_x_continuous(limits = c(199653, 201853),
                       breaks = seq(199653,201853,1100),
                       expand = FALSE) +
    scale_y_continuous(limits = c(288725, 290925),
                       breaks = seq(288725, 290925, 1100)) +
    geom_vline(xintercept=c(200753),
               linetype="dashed", color = "grey3", linewidth=0.1) +
    geom_hline(yintercept=c(289825),
               linetype="dashed", color = "grey3", linewidth=0.1) +
    theme_bw(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                            panel.background = element_blank()) +
    theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) +
    theme(legend.position="none") + 
    theme(aspect.ratio = 1) + 
    labs(
      title = paste0("MATa-MATa' discordant reads pairs r18 - ", strain_name, " - ",sample),
      x = "Chr.III coordinates",
      y = "Chr.V coordinates"
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("MATa_MATa_pairs_r18_", strain_name, "_", sample, ".svg")),
    plot = p,
    width = 8,
    height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})




# Bar plot distribution of polymorphisms vs HOinc reads
MAT_r18_count_HOinc_reads_summary<- MAT_r18_count_HOinc_reads_summary %>%
  mutate(polymorphism = factor(polymorphism, levels = polymorphism_levels)) %>%
  filter(polymorphism != "NA")

split_MAT_df <- split(MAT_r18_count_HOinc_reads_summary, 
                      MAT_r18_count_HOinc_reads_summary$sample)


log_step("Plotting...")
# Generate plots for MATa-MATa' pairs distribution
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  
  p <- ggplot(data_subset, aes(x = polymorphism, y = average_percentage, fill = chromosome)) +
    geom_col(position = "dodge2") +
    geom_errorbar(aes(ymin = average_percentage - sd_percentage, ymax = average_percentage + sd_percentage), 
                  linewidth = 0.8, width = 0.5, colour = "gray10", position = position_dodge(width = 0.9)) +
    scale_fill_manual(values=c("#2F2C7E", "#A30000")) +
    #scale_fill_manual(values=c("#081D5C", "#B86029")) +
    theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                 panel.background = element_blank(), 
                                                 plot.background = element_rect(fill = "transparent", colour = NA)) +
    theme(legend.position="right") +  
    coord_cartesian(ylim = c(0, 100),expand=FALSE) +
    theme(aspect.ratio = 1) + 
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(name = expression("Percentage"),
                       #limits = c(0, 100),
                       breaks = seq(0,100,10)) +
    theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 15),  #25
          axis.title.y = element_text(vjust = 1, size = 15)) + 
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 10, angle = 90), #20
          axis.text.y = element_text(vjust = 0, size = 10)) +
    labs(
      title = paste0("MATa-MATa' poly distribution vs HOinc r18 - ", strain_name, " - ",sample),
      x = "Polymorphism",
      y = "Percentage"
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("MATa_MATa_poly_distrib_vs_HOinc_r18_", strain_name, "_", sample, ".svg")),
    plot = p,
    #width = 8,
    #height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})

# Bar plot distribution of polymorphisms vs total reads for III or V
MAT_r18_count_HOinc_reads_summary<- MAT_r18_count_total_reads_summary %>%
  mutate(polymorphism = factor(polymorphism, levels = polymorphism_levels)) %>%
  filter(polymorphism != "NA")

split_MAT_df <- split(MAT_r18_count_HOinc_reads_summary, 
                      MAT_r18_count_HOinc_reads_summary$sample)


log_step("Plotting...")
# Generate plots for MATa-MATa' pairs distribution
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  
  p <- ggplot(data_subset, aes(x = polymorphism, y = average_percentage, fill = chromosome)) +
    geom_col(position = "dodge2") +
    geom_errorbar(aes(ymin = average_percentage - sd_percentage, ymax = average_percentage + sd_percentage), 
                  linewidth = 0.8, width = 0.5, colour = "gray10", position = position_dodge(width = 0.9)) +
    scale_fill_manual(values=c("#2F2C7E", "#A30000")) +
    #scale_fill_manual(values=c("#081D5C", "#B86029")) +
    theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                 panel.background = element_blank(), 
                                                 plot.background = element_rect(fill = "transparent", colour = NA)) +
    theme(legend.position="right") +  
    coord_cartesian(ylim = c(0, 100), expand=FALSE) +
    theme(aspect.ratio = 0.45) + 
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(name = expression("Percentage"),
                       #limits = c(0, 100),
                       breaks = seq(0,100,20)) +
    theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 20), 
          axis.title.y = element_text(vjust = 1, size = 20)) + 
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 15, angle = 90), 
          axis.text.y = element_text(vjust = 0, size = 15)) +
    labs(
      title = paste0("MATa-MATa' poly distribution vs total r18 - ", strain_name, " - ",sample),
      x = "Polymorphism",
      y = "Percentage"
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("MATa_MATa_poly_distrib_vs_total_r18_", strain_name, "_", sample, ".svg")),
    plot = p,
    #width = 8,
    #height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})

# Bar plot of GC or CO

split_MAT_df <- split(GC_CO_combined_summary, 
                      GC_CO_combined_summary$sample)


log_step("Plotting...")
# Generate plots for GC-CO
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  
  p <- ggplot(data_subset, aes(x = GC_CO, y = average_percentage, fill = GC_CO)) +
    geom_col(position = "dodge2") +
    geom_errorbar(aes(ymin = average_percentage - sd_percentage, ymax = average_percentage + sd_percentage), 
                  linewidth = 0.8, width = 0.5, colour = "gray10", position = position_dodge(width = 0.9)) +
    #scale_fill_manual(values=c("black", "chartreuse3")) +
    scale_fill_manual(values = c("CO" = "black", "GC" = "chartreuse3")) + 
    theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                 panel.background = element_blank(), 
                                                 plot.background = element_rect(fill = "transparent", colour = NA)) +
    theme(legend.position="right") +  
    coord_cartesian(ylim = c(0, 100),expand=FALSE) +
    theme(aspect.ratio = 1) + 
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(name = expression("Percentage"),
                       #limits = c(0, 100),
                       breaks = seq(0,100,10)) +
    theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
          axis.title.y = element_text(vjust = 1, size = 25)) + 
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 20, angle = 90), 
          axis.text.y = element_text(vjust = 0, size = 20)) +
    labs(
      title = paste0("GC vs CO percentage - ", strain_name, " - ",sample),
      x = "Repair product",
      y = "Percentage"
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("MATa_MATa_GC_CO_r18_", strain_name, "_", sample, ".svg")),
    plot = p,
    #width = 8,
    #height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})

# Bar plot of GC groups

split_MAT_df <- split(GC_groups_combined_summary, 
                      GC_groups_combined_summary$sample)


log_step("Plotting...")
# Generate plots for GC-groups
MAT_plot_list <- lapply(names(split_MAT_df), function(sample) {
  data_subset <- split_MAT_df[[sample]]
  
  p <- ggplot(data_subset, aes(x = GC_group, y = average_percentage, fill = GC_group)) +
    geom_col(position = "dodge2") +
    geom_errorbar(aes(ymin = average_percentage - sd_percentage, ymax = average_percentage + sd_percentage), 
                  linewidth = 0.8, width = 0.5, colour = "gray10", position = position_dodge(width = 0.9)) +
    #scale_fill_manual(values=c("black", "chartreuse3")) +
    #scale_fill_manual(values = c("CO" = "black", "GC" = "chartreuse3")) + 
    theme_classic(base_family = "Arial") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                 panel.background = element_blank(), 
                                                 plot.background = element_rect(fill = "transparent", colour = NA)) +
    theme(legend.position="right") +  
    coord_cartesian(ylim = c(0, 100), expand=FALSE) +
    theme(aspect.ratio = 1) + 
    scale_x_discrete(drop = FALSE) +
    scale_y_continuous(name = expression("Percentage"),
                       #limits = c(0, 100),
                       breaks = seq(0,100,10)) +
    theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
          axis.title.y = element_text(vjust = 1, size = 25)) + 
    theme(axis.text.x = element_text(hjust = 1, vjust = 0.5, size = 20, angle = 90), 
          axis.text.y = element_text(vjust = 0, size = 20)) +
    labs(
      title = paste0("GC group percentage - ", strain_name, " - ",sample),
      x = "Repair product",
      y = "Percentage"
    )
  
  # Save SVG
  ggsave(
    filename = file.path(root_dir, paste0("MATa_MATa_GC_group_r18_", strain_name, "_", sample, ".svg")),
    plot = p,
    #width = 8,
    #height = 3.6,
    device = svglite,
    bg = "transparent"
  )
  
  return(p)
  
})