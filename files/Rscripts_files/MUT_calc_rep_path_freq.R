#### R script to calculate repair pathway frequency
### Loop for each strain in MYWD
### 07/03/2026 - Lydia

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
  log_step("Usage: script.R <root_dir> <strain> <wd_dir>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain   <- args[2]
wd_dir   <- args[3]

strain <- sub("/$", "", strain)
strain_name <- basename(strain)

log_step(paste("ROOT_DIR:", root_dir))
log_step(paste("STRAIN:", strain))
log_step(paste("WD_DIR:", wd_dir))


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


# Define a function to process tsv files
process_nucleotide_row_count_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  # directory name above the file
  
  # 
  parts <- str_split(file_base, "_", simplify = TRUE)
  
  # Validate and extract parts safely
  if (ncol(parts) >= 3) {
    sample_name <- parts[1]
    experiment_name <- parts[2]
    nucleotide_name <- parts[3]
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
  
  # Prepare nucleotide_count df
  nucleotide_count_file <- temp_file %>%
    rename("nucleotide_count" = !!names(.[1])) %>%
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name,
      nucleotide = nucleotide_name
    )
  
  return(nucleotide_count_file)
}



log_step("Finding nucleotide row count files...")
# Get all nucleotide row count files recursively in root folder
nucleotide_row_count_files <- list.files(
  path = root_dir,
  pattern = "row_count\\.tsv$",
  recursive = TRUE,
  full.names = TRUE
)





log_step("Processing nucleotide row count files...")
nucleotide_processed_df <- purrr::map_dfr(nucleotide_row_count_files, process_nucleotide_row_count_files)
nucleotide_processed_df <- nucleotide_processed_df %>% 
  mutate(repair_pathway = ifelse(nucleotide == "200689AC", "NHEJ", 
                                 ifelse(nucleotide == "200689ACT", "shHR+NHEJ", 
                                        ifelse(nucleotide == "200689ACTG", "all", 
                                               ifelse(nucleotide == "200689AT", "shHR", 
                                                      ifelse(nucleotide == "200689G", "lHR", "-"))))))

nucleotide_processed_df_all <- nucleotide_processed_df %>% filter(repair_pathway == "all") %>% 
  mutate(nucleotide_count_all = nucleotide_count) %>% 
  select(strain, sample, experiment, nucleotide_count_all)


nucleotide_processed_df_ratio <- left_join(nucleotide_processed_df, nucleotide_processed_df_all) %>% 
  mutate(nucleotide_count_ratio = nucleotide_count / nucleotide_count_all) %>% 
  filter(repair_pathway %in% c("NHEJ", "shHR", "lHR")) 




nucleotide_processed_df_ratio_T0 <- nucleotide_processed_df_ratio %>% 
  filter(sample == "T0" & repair_pathway == "lHR") %>% 
  mutate(nucleotide_count_ratio_T0 = nucleotide_count_ratio) %>% 
  select(strain, experiment, nucleotide_count_ratio_T0)


nucleotide_processed_df_ratio_lHR <- nucleotide_processed_df_ratio %>% 
  filter(sample == "TLG" & repair_pathway == "lHR") %>% 
  left_join(., nucleotide_processed_df_ratio_T0) %>% 
  mutate(nucleotide_count_ratio_repair_V = 2*(nucleotide_count_ratio - nucleotide_count_ratio_T0)) %>% 
  select(strain, sample, experiment, nucleotide, repair_pathway, nucleotide_count_ratio_repair_V) %>% 
  rename("nucleotide_count_ratio" = !!names(.[6]))

nucleotide_processed_df_ratio_NHEJ <- nucleotide_processed_df_ratio %>% 
  filter(sample == "TLG" & repair_pathway == "NHEJ") %>% 
  select(strain, sample, experiment, nucleotide, repair_pathway, nucleotide_count_ratio) 

nucleotide_processed_df_ratio_shHR <- nucleotide_processed_df_ratio %>% 
  filter(sample == "TLG" & repair_pathway == "shHR") %>% 
  select(strain, sample, experiment, nucleotide, repair_pathway, nucleotide_count_ratio) 

nucleotide_processed_df_ratio_repair <- bind_rows(nucleotide_processed_df_ratio_lHR, nucleotide_processed_df_ratio_shHR) %>% 
  bind_rows(., nucleotide_processed_df_ratio_NHEJ)




write_tsv(nucleotide_processed_df_ratio, file.path(root_dir, paste0(strain_name, "_repair_pathway_df.tsv")))
write_tsv(nucleotide_processed_df_ratio_repair, file.path(root_dir, paste0(strain_name, "_repair_pathway_only_repair_df.tsv")))


nucleotide_processed_df_ratio_summary <- nucleotide_processed_df_ratio %>% 
  group_by(strain, sample, nucleotide, repair_pathway) %>% 
  summarise(mean_ratio = mean(nucleotide_count_ratio, na.rm = TRUE), 
            sd_ratio = sd(nucleotide_count_ratio, na.rm = TRUE))


nucleotide_processed_df_ratio_repair_summary <- nucleotide_processed_df_ratio_repair %>% 
  group_by(strain, sample, nucleotide, repair_pathway) %>% 
  summarise(mean_ratio = mean(nucleotide_count_ratio, na.rm = TRUE), 
            sd_ratio = sd(nucleotide_count_ratio, na.rm = TRUE))


nucleotide_processed_df_ratio_repair_fraction_summary <- nucleotide_processed_df_ratio_repair %>% 
  group_by(strain, sample, experiment) %>% 
  summarise(sum_ratio = sum(nucleotide_count_ratio, na.rm = TRUE)) %>%  
  summarise(mean_ratio = mean(sum_ratio, na.rm = TRUE), 
            sd_ratio = sd(sum_ratio, na.rm = TRUE))





write_tsv(nucleotide_processed_df_ratio_summary, file.path(root_dir, paste0(strain_name, "_repair_pathway_summary.tsv")))
write_tsv(nucleotide_processed_df_ratio_repair_summary, file.path(root_dir, paste0(strain_name, "_repair_pathway_only_repair_summary.tsv")))
write_tsv(nucleotide_processed_df_ratio_repair_fraction_summary, file.path(root_dir, paste0(strain_name, "_repair_pathway_only_repair_fraction_summary.tsv")))

log_step("Preparing plot data...")
samples_order <- c("T0", "TLG", "TLR")
repair_pathway_order <- c("NHEJ", "shHR", "lHR")

nucleotide_processed_df_ratio_summary_graphs <- nucleotide_processed_df_ratio_summary %>% 
  mutate(sample = factor(sample, levels = samples_order),
         repair_pathway = factor(repair_pathway, levels = repair_pathway_order)) 


nucleotide_processed_df_ratio_repair_summary_graphs <- nucleotide_processed_df_ratio_repair_summary %>% 
  mutate(repair_pathway = factor(repair_pathway, levels = repair_pathway_order)) 

log_step("Plotting...")
# Generate plots

bar_plot <- ggplot(nucleotide_processed_df_ratio_summary_graphs, aes(x = sample, y = mean_ratio, fill = repair_pathway)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio), 
                linewidth = 0.8, width = 0.2, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("gold","dodgerblue", "purple3")) +
  theme_classic(base_family = "Helvetica") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                   panel.background = element_blank(), 
                                                   plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  coord_cartesian(ylim = c(0, 1), expand=FALSE) +
  theme(aspect.ratio = 0.75) + 
  scale_x_discrete(name = expression("Sample")) +
  scale_y_continuous(name = expression("Fraction"),
                     #limits = c(0, 1),
                     breaks = seq(0,1,0.2)) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 20, angle = 0), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  labs(
    title = paste0("Repair_pathway - ", strain_name)
  )

ggsave(
  filename = paste0(strain,"/", "Repair_pathway_plot_", strain_name, ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)


log_step("Plotting...")
bar_plot <- ggplot(nucleotide_processed_df_ratio_repair_summary_graphs, aes(x = sample, y = mean_ratio, fill = repair_pathway)) +
  geom_col(position = "dodge2") +
  geom_errorbar(aes(ymin = mean_ratio - sd_ratio, ymax = mean_ratio + sd_ratio), 
                linewidth = 0.8, width = 0.2, colour = "gray10", position = position_dodge(width = 0.9)) +
  scale_fill_manual(values=c("gold","dodgerblue", "purple3")) +
  theme_classic(base_family = "Helvetica") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                   panel.background = element_blank(), 
                                                   plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="right") +  
  coord_cartesian(expand=FALSE) +
  coord_cartesian(ylim = c(0, 1), expand=FALSE) +
  theme(aspect.ratio = 0.75) + 
  scale_x_discrete(name = expression("Sample")) +
  scale_y_continuous(name = expression("Fraction"),
                     #limits = c(0, 1),
                     breaks = seq(0,1,0.2)) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
        axis.title.y = element_text(vjust = 1, size = 25)) + 
  theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 20, angle = 0), 
        axis.text.y = element_text(vjust = 0, size = 20)) +
  labs(
    title = paste0("Repair_pathway - repair - ", strain_name)
  )

ggsave(
  filename = paste0(strain,"/", "Repair_pathway_only_repair_plot_", strain_name, ".svg"),
  plot = bar_plot,
  #width = 8,
  #height = 3.6,
  device = svglite,
  bg = "transparent"
)



# Calculate % of each nucleotide
# Define a function to process wig files
process_wig_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  
  
  # Remove extension
  file_no_ext <- str_remove(file_base, "\\.wig$")
  
  # Extract sample and experiment
  
  parts <- str_split(file_no_ext, "_", simplify = TRUE)
  
  if (ncol(parts) < 2) {
    warning(paste("Filename does not match expected format:", file_base))
    return(NULL)
  }
  
  sample_name <- parts[1]
  experiment_name <- parts[2]  
  
  # Extract everything after experiment (E1_, E2_, E3_)
  nucleotide_name <- str_remove(file_no_ext, paste0("^", sample_name, "_", experiment_name, "_"))
  
  # ---- Read the file ----
  all_lines <- readLines(file_path)
  data_lines <- all_lines[!grepl("^(track|#|variableStep)", all_lines)]
  
  if (length(data_lines) == 0) {
    return(NULL)
  }
  
  wig_df <- read_tsv(
    I(data_lines),
    col_names = c("Pos", "A", "C", "G", "T", "N", "DEL", "INS"),
    show_col_types = FALSE
  )
  
  if (nrow(wig_df) == 0) {
    return(NULL)
  }
  
  # ---- Add metadata ----
  wig_df <- wig_df %>%
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name,
      nucleotide = nucleotide_name
    )
  
  return(wig_df)
}



log_step("Finding wig files...")
# Get all wig files recursively in root folder
wig_files <- list.files(
  path = root_dir,
  pattern = "\\.wig$",
  recursive = TRUE,
  full.names = TRUE
)


wig_processed_df <- purrr::map_dfr(wig_files, process_wig_files)

wig_processed_df <- wig_processed_df %>% 
  mutate(repair_pathway = ifelse(nucleotide == "200689_A_200753_C", "NHEJ", 
                                 ifelse(nucleotide == "200689_A_CT", "shHR+NHEJ", 
                                        ifelse(nucleotide == "200689_A_CT_G", "all", 
                                               ifelse(nucleotide == "200689_A_200753_T", "shHR", 
                                                      ifelse(nucleotide == "200689_G", "lHR", "-"))))))


wig_processed_df_percent <- wig_processed_df %>% 
  pivot_longer(cols = c(A, C, T, G, N, DEL, INS),
               names_to = "mutation",
               values_to = "mutation_count") %>% 
  group_by(Pos, strain, sample, experiment, nucleotide, repair_pathway) %>%
  mutate(
    total_count = sum(mutation_count),
    perc_mutation = (mutation_count / total_count) * 100
  ) %>%
  ungroup() 

write_tsv(wig_processed_df_percent, file.path(root_dir, paste0( "/",strain_name,  "_mutation_wig_percentage_df.tsv")))


wig_processed_df_percent_summary <- wig_processed_df_percent %>% 
  group_by(Pos, strain, sample, nucleotide, repair_pathway, mutation) %>% 
  summarise(mean_percentage = mean(perc_mutation, na.rm = TRUE), 
            sd_percentage = sd(perc_mutation, na.rm = TRUE)) %>% 
  ungroup() %>% 
  filter(mutation != "N")

write_tsv(wig_processed_df_percent_summary, file.path(root_dir, paste0("/",strain_name, "_mutation_wig_percentage_summary.tsv")))


# Calculate % of AF
# Define a function to process BCF files
process_bcf_files <- function(file_path) {
  # Extract filename and directory parts
  file_base <- basename(file_path)
  strain_name <- basename(dirname(file_path))  
  
  # Remove extension
  file_no_ext <- str_remove(file_base, "_BCF\\.tsv$")
  
  # Extract sample and experiment
  
  parts <- str_split(file_no_ext, "_", simplify = TRUE)
  
  if (ncol(parts) < 2) {
    warning(paste("Filename does not match expected format:", file_base))
    return(NULL)
  }
  
  sample_name <- parts[1]
  experiment_name <- parts[2]  
  
  # Extract everything after experiment (E1_, E2_, E3_)
  nucleotide_name <- str_remove(file_no_ext, paste0("^", sample_name, "_", experiment_name, "_"))
  
  # ---- Read the file ----
  all_lines <- readLines(file_path)
  data_lines <- all_lines[!grepl("^(track|#|variableStep)", all_lines)]
  
  if (length(data_lines) == 0) {
    return(NULL)
  }
  
  bcf_df <- read_tsv(
    I(data_lines),
    col_select = c(X2, X6),
    col_names = FALSE,
    show_col_types = FALSE
  )
  
  if (nrow(bcf_df) == 0) {
    return(NULL)
  }
  
  # ---- Add metadata ----
  bcf_df <- bcf_df %>%
    mutate(
      strain = strain_name,
      sample = sample_name,
      experiment = experiment_name,
      nucleotide = nucleotide_name
    ) %>% 
    rename("Pos" = !!names(.[1]), "AF_perc" = !!names(.[2]))
  
  return(bcf_df)
}



log_step("Finding bcf files...")
# Get all BCFv files recursively in root folder
bcf_files <- list.files(
  path = root_dir,
  pattern = "\\BCF.tsv",
  recursive = TRUE,
  full.names = TRUE
)


bcf_processed_df <- purrr::map_dfr(bcf_files, process_bcf_files) 


bcf_processed_df <- bcf_processed_df %>% 
  mutate(repair_pathway = ifelse(nucleotide == "A_200753_C", "NHEJ", 
                                 ifelse(nucleotide == "A_CT", "shHR+NHEJ", 
                                        ifelse(nucleotide == "A_CT_G", "all", 
                                               ifelse(nucleotide == "A_200753_T", "shHR", 
                                                      ifelse(nucleotide == "G", "lHR", "-")))))) %>% 
  mutate(mutation = "AF")


bcf_processed_df_percent_summary <- bcf_processed_df %>% 
  group_by(Pos, strain, sample, nucleotide, mutation, repair_pathway) %>%
  summarise(mean_percentage = mean(AF_perc, na.rm = TRUE), 
            sd_percentage = sd(AF_perc, na.rm = TRUE)) %>% 
  ungroup() 

write_tsv(bcf_processed_df, file.path(root_dir, paste0(strain_name, "_mutation_bcf_percentage_df.tsv")))
write_tsv(bcf_processed_df_percent_summary, file.path(root_dir, paste0(strain_name, "_mutation_bcf_percentage_summary.tsv")))

wig_processed_df_percent_summary_plots <- wig_processed_df_percent_summary %>% 
  select(c(Pos,strain, sample, repair_pathway, mutation, mean_percentage, sd_percentage))

bcf_processed_df_summary_plots <- bcf_processed_df_percent_summary %>% 
  select(c(Pos,strain, sample, repair_pathway, mutation, mean_percentage, sd_percentage))


wig_bcf_combined_summary <- bind_rows(wig_processed_df_percent_summary_plots, bcf_processed_df_summary_plots) %>% 
  select(Pos,strain, sample, repair_pathway, mutation, mean_percentage, sd_percentage) %>% 
  filter(mutation != "N") 

write_tsv(wig_bcf_combined_summary, file.path(root_dir, paste0("/",strain_name, "_mutation_wig_bcf_percentage_summary.tsv")))





mutation_order <- rev(c("A", "C", "G", "T", "AF", "DEL", "INS"))

wig_bcf_plot_df <- wig_bcf_combined_summary %>% 
  mutate(rel_pos = Pos - 200753) %>% 
  mutate(mutation = factor(mutation, levels = mutation_order))




log_step("Plotting heatmaps...")

wig_bcf_plot_df %>%
  group_split(repair_pathway, sample) %>%
  purrr::imap(function(df, i) {
    rp <- unique(df$repair_pathway)
    sm <- unique(df$sample)
    
    p <- ggplot(df, aes(x = rel_pos, y = mutation, fill = mean_percentage)) +
      geom_tile(color = "white", linewidth = 0.4) +
      geom_hline(yintercept = seq(0.5, length(unique(df$mutation)) + 0.5, by = 1),
                 color = "black", linewidth = 0.2) +
      geom_vline(xintercept = seq(-11.5, 9.5, by = 1),
                 color = "black", linewidth = 0.2) +
      scale_x_continuous(breaks = seq(-11, 9, by = 1), expand = c(0, 0)) +
      scale_fill_gradientn(
        colors = rev(c("#AF2418", "#E1AC40", "#EFD24D",
                                "#5D8B27", "#4EACE9", "#4573A1", "white")),
                                values = scales::rescale(c(0, 25, 50, 75, 100)),
        na.value = "gray90",
        name = "Normalized Coverage",
        limits = c(0, 100),
        oob = scales::squish
      ) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste0(rp, " - ", sm, " - ", strain_name),
        x = "Genomic position",
        y = "Base/Indel"
      )
    #ggtitle(paste("Repair pathway:", rp, "| Sample:", sm)) +
    
    
    ggsave(
      filename = paste0(strain,"/", "Repair_heatmap_plot_",rp, "_", sm, "_", strain_name, ".svg"),
      plot = p,
      width = 8,
      height = 3,
      device = svglite,
      bg = "transparent"
    )
    
  })

log_step("Plotting DNA logos...")

library(ggseqlogo)

wig_bcf_plot_df_logo <- wig_bcf_plot_df %>% 
  select(strain, sample, repair_pathway, mutation, mean_percentage, rel_pos) %>% 
  filter(mutation %in% c("A", "C", "T", "G")) %>%
  filter(sample == "TLG")

cs1 = make_col_scheme(chars=c('C', 'A', 'G', 'T'), 
                      cols=c('#355C95', '#459450', '#ECB74B', '#C53940'))



wig_bcf_plot_df_logo%>%
  group_split(repair_pathway, sample) %>%
  purrr::imap(function(df, i) {
    
    rp <- unique(df$repair_pathway)
    sm <- unique(df$sample)
    
    
    matrix <- df %>%
      mutate(freq = mean_percentage / 100) %>%        
      select(rel_pos, mutation, freq) %>%
      pivot_wider(names_from = rel_pos, values_from = freq) %>%
      column_to_rownames("mutation") %>%
      as.matrix()
    
    # ensure nucleotide order is standard A,C,G,T if needed
    matrix <- matrix[c("A", "C", "G", "T"), , drop = FALSE]
    
    # make the ggseqlogo plot
    p <- ggseqlogo(matrix, method = "custom", seq_type='dna',
                   stack_width = 0.9, font = "helvetica_regular", col_scheme=cs1) +
      theme_classic(base_family = "Helvetica") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                       panel.background = element_blank()) +
      ggtitle(paste0(rp, " - ", sm, " - ", strain_name)) +
      theme(aspect.ratio = 0.3) + 
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank()
      ) +
      xlab("Genomic position") +
      ylab("Freq") +
      scale_x_continuous(
        breaks = seq_along(colnames(matrix)),
        labels = colnames(matrix)
      )
    ggsave(
      filename = paste0(strain,"/", "Repair_logo_plot_",rp, "_", sm, "_", strain_name, ".svg"),
      plot = p,
      width = 8,
      height = 3,
      device = svglite,
      bg = "transparent"
    )
    
    
  })

#### Mutagenic rate

wig_processed_df_percent_mut_rate <- wig_processed_df_percent %>% 
  select(Pos, strain, sample, experiment, repair_pathway, mutation, perc_mutation) %>% 
  filter(mutation %in% c("DEL", "INS")) %>% 
  filter(Pos != 200753) %>% 
  group_by(strain, sample, experiment, repair_pathway, mutation) %>% 
  summarise(mean_mut_rate = mean(perc_mutation))


bcf_processed_df_mut_rate <- bcf_processed_df %>% 
  select(Pos, strain, sample, experiment, repair_pathway, mutation, AF_perc) %>% 
  mutate(perc_mutation = AF_perc) %>% 
  select(Pos, strain, sample, experiment, repair_pathway, mutation, perc_mutation) %>% 
  filter(mutation == "AF") %>% 
  filter(Pos != 200753) %>% 
  group_by(strain, sample, experiment, repair_pathway, mutation) %>% 
  summarise(mean_mut_rate = mean(perc_mutation))



wig_bcf_combined_summary_mut_rate <- bind_rows(wig_processed_df_percent_mut_rate, 
                                               bcf_processed_df_mut_rate) %>% 
  filter(sample == "TLG") %>% 
  group_by(strain, sample, repair_pathway, mutation) %>% 
  summarise(mean_mut_rate_global = mean(mean_mut_rate),
            sd_mut_rate_global = sd(mean_mut_rate))

write_tsv( wig_bcf_combined_summary_mut_rate, file.path(root_dir, paste0("/",strain_name, "_mutation_wig_bcf_mut_rate_summary.tsv")))


mutation_order_plot <- c("AF", "DEL", "INS")

wig_bcf_combined_summary_mut_rate <- wig_bcf_combined_summary_mut_rate %>% 
  mutate(mutation = factor(mutation, levels = mutation_order_plot))



wig_bcf_combined_summary_mut_rate %>% 
  group_split(repair_pathway, sample) %>%
  purrr::imap(function(df, i) {
    rp <- unique(df$repair_pathway)
    sm <- unique(df$sample)
    
    p <- ggplot(df, aes(x = strain, y = mean_mut_rate_global, fill = mutation)) +
      geom_col(position = "dodge2") +
      geom_errorbar(aes(ymin = mean_mut_rate_global - sd_mut_rate_global, ymax = mean_mut_rate_global + sd_mut_rate_global), 
                    linewidth = 0.8, width = 0.2, colour = "gray10", position = position_dodge(width = 0.9)) +
      scale_fill_manual(values=c("gold", "dodgerblue", "purple3")) +
      theme_classic(base_family = "Helvetica") + theme(panel.grid = element_line(color = "black", linewidth = 0.1),
                                                       panel.background = element_blank(), 
                                                       plot.background = element_rect(fill = "transparent", colour = NA)) +
      theme(legend.position="right") +  
      coord_cartesian(expand=FALSE) +
      coord_cartesian(ylim = c(0, 10), expand=FALSE) +
      theme(aspect.ratio = 0.7) + 
      scale_x_discrete(name = expression("Strain")) +
      scale_y_continuous(name = expression("Percentage"),
                         #limits = c(0, 1),
                         breaks = seq(0,10,2)) +
      theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 25), 
            axis.title.y = element_text(vjust = 1, size = 25)) + 
      theme(axis.text.x = element_text(hjust = 0.5, vjust = 0, size = 15, angle = 0), 
            axis.text.y = element_text(vjust = 0, size = 20)) +
      labs(
        title = paste0(rp, " - ", sm, " - ", strain_name)
      )
    
    ggsave(
      filename = paste0(strain,"/", "Mutagenic_rate_plot_",rp, "_", sm, "_", strain_name, ".svg"),
      plot = p,
      #width = 8,
      #height = 3,
      device = svglite,
      bg = "transparent"
    )
    
    p
  })
