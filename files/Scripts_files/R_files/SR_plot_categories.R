#### R script to order datasets from genomic categories and plot violin plots, 75nt, average all experiments
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

if (length(args) != 3) {
  log_step("ERROR: Incorrect number of arguments")
  log_step("Usage: script.R <root_dir> <strain> <analysis_suffix>")
  log_step(paste("Received:", paste(args, collapse=" ")))
  quit(status = 1)
}

root_dir <- args[1]
strain <- args[2]
analysis_suffix <- args[3]


strain <- sub("/$", "", strain)
strain_name <- basename(strain)



log_step(paste("STRAIN:", strain))
log_step(paste("ANALYSIS:", analysis_suffix))



# load libraries

library(ggplot2)
library(extrafont)
library(svglite)
library(purrr)
library(stringr)
library(readr)
library(plotrix)
library(tidyverse, warn.conflicts = FALSE)
library(tidyr, warn.conflicts = FALSE)
library(dplyr, warn.conflicts = FALSE)
options(dplyr.summarise.inform = FALSE)





# REORDER THE DATASET BY CATEGORIES
# Define the order of categories
categoria_order <- c("ORF", "intergenic", "long_terminal_repeat", "transposable_element_gene", 
                     "LTR_retrotransposon", "tRNA_gene", "rRNA_gene", "ncRNA_gene", 
                     "snRNA_gene", "snoRNA_gene", "ARS", "centromere", "telomere")

# Define the file pattern to search for
genomic_sorted_path <- paste0(root_dir, "/Genomic_sorted_", strain_name, "_", analysis_suffix, ".tsv")

# Check if the file exists
if (!file.exists(genomic_sorted_path)) {
  stop("File not found: ", genomic_sorted_path)
}

# Read the file
genomic_sorted_df <- read_tsv(genomic_sorted_path, show_col_types = FALSE)

# Factor and sort by category
category_sorted_df <- genomic_sorted_df %>%
  mutate(Category_A = factor(Category_A, levels = categoria_order, ordered = TRUE)) %>%
  arrange(Category_A)

write_tsv(category_sorted_df, paste0(root_dir, "/Category_sorted_", strain_name, "_", analysis_suffix, ".tsv"))



# PLOT VIOLIN PLOTS WITH THE DISTRIBUTION OF THE DATA FOR EACH GENOMIC CATEGORY



# Define the transformation function # for better visualization
transform_y <- function(y) {
  ifelse(y <= 2, 
         y * (0.75 / 2),                     # scale 0-2 to 0-0.75
         0.75 + ((y - 2) * (0.25 / (10 - 2))) # scale 2-10 to 0.75-1
  )
}

# Transformed ratio_medio column
category_sorted_df$avg_ratio_coverage_feature_trans <- transform_y(category_sorted_df$avg_ratio_coverage_feature)

# # Create output folder
out_dir <- file.path(root_dir, "/Plots_Analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Define the color palette with your specific color codes
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

p1 <- ggplot(category_sorted_df, aes(x = Category_A, y = avg_ratio_coverage_feature_trans, fill = Category_A)) +
  geom_jitter(aes(color = Category_A), width = 0.1, size = 0.1, alpha = 0.2) +
  geom_violin(trim = FALSE, size = 0.1, width = 0.8, alpha = 0.8, scale = "width") +
  geom_boxplot(width = 0.1, alpha = 1, size = 0.1, outlier.shape = NA) +
  theme_minimal() +
  labs(title = paste("Strain:", strain_name, "-", analysis_suffix),
       x = "Category", y = "Coverage", fill = "Category") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values = categoria_colors) +
  scale_color_manual(values = categoria_colors, guide = "none") + 
  scale_x_discrete(labels = c(
    "ORF" = "ORF",
    "intergenic" = "Intergenic",
    "long_terminal_repeat" = "LTR",
    "transposable_element_gene" = "TEG",
    "LTR_retrotransposon" = "Ty",
    "tRNA_gene" = "tRNA",
    "rRNA_gene" = "rRNA",
    "ncRNA_gene" = "ncRNA",
    "snRNA_gene" = "snRNA",
    "snoRNA_gene" = "snoRNA",
    "ARS" = "ARS",
    "centromere" = "Centromere",
    "telomere" = "Telomere"
  )) +
  coord_cartesian(ylim = c(0, 1.02), expand = FALSE) +
  scale_y_continuous(
    name = "Coverage",
    breaks = c(0, 0.1875, 0.375, 0.5625, 0.75, 0.8125, 0.875, 0.9375, 1),
    labels = c("0", "0.5", "1", "1.5", "2", "4", "6", "8", "10")
  )

# Save the plot
formats <- c("svg", "png")
for (fmt in formats) {
  ggsave(
    filename = file.path(out_dir, paste0("Violin_", strain_name, "_", analysis_suffix, ".", fmt)),
    plot = p1,
    device = fmt,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# PLOTTING GENOMIC COVERAGE

# Prepare the template dataframe: Generate column1, total number of rows in the dataframe
total_rows <- nrow(category_sorted_df)

# First segment: 0.5 increments
first_part <- seq(1, 6597.5, by = 0.5)

# Number of remaining positions
remaining <- total_rows - length(first_part)

# Second segment: 1 increments starting from 6599
second_part <- seq(6598.5, by = 1, length.out = remaining)

# Combine both parts
position_vector <- c(first_part, second_part)

# Create the final dataframe
Cov_df <- data.frame(
  Position = position_vector,
  Coverage = category_sorted_df$avg_ratio_coverage_feature,
  Coverage_eje = category_sorted_df$avg_ratio_coverage_feature - 1,
  SD = category_sorted_df$sd_ratio_coverage_feature,
  Series = category_sorted_df$Category_A
)

# Output filename
cov_out_name <- paste0(root_dir, "/Coverage_Template_Data_", strain_name, "_", analysis_suffix, ".tsv")

# Write the new TSV
write_tsv(Cov_df, file = cov_out_name)
message("Written: ", cov_out_name)

# Plot genomic coverage
p2 <- ggplot(Cov_df, aes(x = Position, y = Coverage_eje)) +
  geom_col(size = 0.1, linewidth = 0.1, color = "darkblue") +
  theme_classic(base_family = "Arial") +
  theme(panel.grid = element_line(color = "black", linewidth = 0.1),
        panel.background = element_blank(), 
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  theme(legend.position="none") +
  labs(title = paste("Strain:", strain_name, "-", analysis_suffix),
       x = "Position", y = "Coverage (avg_ratio)") +
  coord_cartesian(xlim = c(1,7921.5), ylim = c(-1, 1), expand=FALSE) +
  theme(aspect.ratio = 0.08) + 
  scale_x_continuous(breaks = c(1, 3285.5, 6598.5, 6981.5, 7072.5, 7122.5, 7397.5,
                                7422.5, 7439.5, 7445.5, 7522.5, 7874.5, 7890.5)) +
  scale_y_continuous(breaks = seq(-1, 1, 0.5),
                     labels = c("0", "0.5", "1", "1.5", "2")) +
  geom_vline(xintercept=c(1, 59.5, 287.5, 378, 796, 957.5, 1027, 1318.5, 1479,
                          1599.5, 1798.5, 1972.5, 2261.5, 2514, 2731.5, 3030),
             linetype="dashed", color = "black", linewidth=0.1) +
  geom_vline(xintercept=c(3285.5, 3342, 3567.5, 3669, 4081, 4247, 4319.5, 4623, 4781.5,
                          4899.5, 5101.5, 5284.5, 5564.5, 5824.5, 6040.5, 6343.5),
             linetype="dashed", color = "black", linewidth=0.1) +
  geom_vline(xintercept=c(3285.5, 6598.5, 6981.5, 7072.5, 7122.5, 7397.5,
                          7422.5, 7439.5, 7445.5, 7522.5, 7874.5, 7890.5),
             linetype="dashed", color = "black", linewidth=0.15) +
  theme(axis.title.x = element_text(hjust = 1, vjust = 0, size = 10),
        axis.title.y = element_text(vjust = 1, size = 10)) +
  theme(axis.text.x = element_text(vjust = 0, size = 0),
        axis.text.y = element_text(vjust = 0, size = 5))

# Save the plot
formats <- c("svg", "png")
for (fmt in formats) {
  ggsave(
    filename = file.path(out_dir, paste0("Genomic_Coverage_", strain_name, "_", analysis_suffix, ".", fmt)),
    plot = p2,
    device = fmt,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# QUANTIFICATION OF COVERAGE >1.2 AND <0.8

Cov_summary <- Cov_df %>%
  group_by(Series) %>%
  summarise(
    total_non_na = sum(!is.na(Coverage)),
    perc_above_1.2 = (sum(Coverage > 1.2, na.rm = TRUE) / total_non_na) * 100,
    perc_below_0.8 = (sum(Coverage < 0.8, na.rm = TRUE) / total_non_na) * 100,
    sd_above_1.2 = mean(SD[Coverage > 1.2], na.rm = TRUE),
    sd_below_0.8 = mean(SD[Coverage < 0.8], na.rm = TRUE)
  )

# Pivot to long format for plotting
Cov_summary_long <- Cov_summary %>%
  select(Series, perc_above_1.2, perc_below_0.8, sd_above_1.2, sd_below_0.8) %>%
  pivot_longer(cols = starts_with("perc_"), names_to = "Category", values_to = "Percentage") %>%
  mutate(
    Category = recode(Category,
                      "perc_above_1.2" = "> 1.2",
                      "perc_below_0.8" = "< 0.8"),
    Percentage = ifelse(Category == "< 0.8", -Percentage, Percentage)
  )

# Add SDs in same long format
Cov_sd_long <- Cov_summary %>%
  select(Series, sd_above_1.2, sd_below_0.8) %>%
  pivot_longer(cols = starts_with("sd_"), names_to = "Category", values_to = "SD") %>%
  mutate(
    Category = recode(Category,
                      "sd_above_1.2" = "> 1.2",
                      "sd_below_0.8" = "< 0.8")
  )

# Merge percentages and SDs
Cov_summary <- left_join(Cov_summary_long, Cov_sd_long, by = c("Series", "Category"))

# Plot the results
# Define the order for the y-axis categories
category_order <- c(
  "ORF", "intergenic", "long_terminal_repeat", "transposable_element_gene",
  "LTR_retrotransposon", "tRNA_gene", "rRNA_gene", "ncRNA_gene",
  "snRNA_gene", "snoRNA_gene", "ARS", "centromere", "telomere"
)
# Relevel Series factor in the summary dataframe
Cov_summary$Series <- factor(Cov_summary$Series, levels = rev(category_order))

p3 <- ggplot(Cov_summary, aes(x = Percentage, y = Series, fill = Category)) +
  geom_col(width = 0.6) +
  geom_errorbarh(aes(xmin = Percentage - SD, xmax = Percentage + SD, color = Category),
                 height = 0.3, linewidth = 0.25) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
  scale_color_manual(values = c("> 1.2" = "#A30000", "< 0.8" = "#2F2C7E")) +
  scale_fill_manual(values = c("> 1.2" = "#A30000", "< 0.8" = "#2F2C7E")) +
  scale_x_continuous(labels = abs,
                     name = "Percentage",
                     limits = c(-75, 75)) +
  scale_y_discrete(labels = c(
    "ORF" = "ORF",
    "intergenic" = "Intergenic",
    "long_terminal_repeat" = "LTR",
    "transposable_element_gene" = "TEG",
    "LTR_retrotransposon" = "Ty",
    "tRNA_gene" = "tRNA",
    "rRNA_gene" = "rRNA",
    "ncRNA_gene" = "ncRNA",
    "snRNA_gene" = "snRNA",
    "snoRNA_gene" = "snoRNA",
    "ARS" = "ARS",
    "centromere" = "Centromere",
    "telomere" = "Telomere"
  )) +
  labs(title = paste("Strain:", strain_name, "-", analysis_suffix),
       x = "Percentage",
       y = "Category") +
  theme_minimal() +
  theme(
    panel.grid.major.y = element_blank(),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10))
  )

# Save the plot
formats <- c("svg", "png")
for (fmt in formats) {
  ggsave(
    filename = file.path(out_dir, paste0("Percentages <0.8 and >1.2_", strain_name, "_", analysis_suffix, ".", fmt)),
    plot = p3,
    device = fmt,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# QUANTIFICATION OF COVERAGE >1.2 AND <0.8 BY CHROMOSOME

# Reorder the dataset by chromosome
# Define the order of chromosomes
chromosome_order <- c("CHRI", "CHRII", "CHRIII", "CHRIV", "CHRV", "CHRVI", 
                      "CHRVII", "CHRVIII", "CHRIX", "CHRX", "CHRXI", 
                      "CHRXII", "CHRXIII", "CHRXIV", "CHRXV", "CHRXVI")

# Sort Coverage file by chromosme
chromosome_sorted_df <- genomic_sorted_df %>%
  mutate(chromosome = factor(chromosome, levels = chromosome_order, ordered = TRUE)) %>%
  arrange(chromosome)

# Otuput file
chr_name <- paste0(root_dir, "/Chromosome_sorted_", strain_name, "_", analysis_suffix, ".tsv")

# Write Output
write_tsv(chromosome_sorted_df, chr_name)
message("Written: ", chr_name)

# Filtering the data (cov >1.2 and <0.8), determine their percentage against the total values per chromosome
Chr_summary <- chromosome_sorted_df %>%
  group_by(chromosome) %>%
  summarise(
    total_non_na = sum(!is.na(avg_ratio_coverage_feature)),
    above_1.2 = sum(avg_ratio_coverage_feature > 1.2, na.rm = TRUE),
    below_0.8 = sum(avg_ratio_coverage_feature < 0.8, na.rm = TRUE),
    perc_above = (above_1.2 / total_non_na) * 100,
    perc_below = (below_0.8 / total_non_na) * 100,
    sd_above = mean(sd_ratio_coverage_feature[avg_ratio_coverage_feature > 1.2], na.rm = TRUE),
    sd_below = mean(sd_ratio_coverage_feature[avg_ratio_coverage_feature < 0.8], na.rm = TRUE)
  ) %>%
  pivot_longer(
    cols = c(perc_above, perc_below),
    names_to = "Category",
    values_to = "Percentage"
  ) %>%
  mutate(
    SD = ifelse(Category == "perc_above", sd_above, sd_below),
    Category = ifelse(Category == "perc_above", "> 1.2", "< 0.8"),
  ) %>%
  select(chromosome, Category, Percentage, SD)

# Apply negation to the <0.8 values only for p4
Chr_summary_p4 <- Chr_summary %>%
  mutate(Percentage = ifelse(Category == "< 0.8", -Percentage, Percentage))

# Plot the results
# Relevel Chromosome factor in the summary dataframe
Chr_summary_p4$chromosome <- factor(Chr_summary$chromosome, levels = (chromosome_order))

p4 <- ggplot(Chr_summary_p4, aes(x = chromosome, y = Percentage, fill = Category)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = Percentage - SD, ymax = Percentage + SD, color = Category),
                width = 0.2, linewidth = 0.3) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  scale_fill_manual(values = c("> 1.2" = "#A30000", "< 0.8" = "#2F2C7E")) +
  scale_color_manual(values = c("> 1.2" = "#A30000", "< 0.8" = "#2F2C7E")) +
  scale_y_continuous(labels = abs,
                     name = "Percentage",
                     limits = c(-5, 15),
                     breaks = seq(-5, 15, 2.5)) +
  scale_x_discrete(
    name = "Chromosome",
    labels = c("I", "II", "III", "IV", "V", "VI", 
               "VII", "VIII", "IX", "X", "XI", 
               "XII", "XIII", "XIV", "XV", "XVI")
  ) +
  labs(title = paste("Strain:", strain_name, "-", analysis_suffix)) +
  theme_minimal() +
  theme(
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "right"
  )

# Save the plot
formats <- c("svg", "png")
for (fmt in formats) {
  ggsave(
    filename = file.path(out_dir, paste0("Chr_Percentages <0.8 and >1.2_", strain_name, "_", analysis_suffix, ".", fmt)),
    plot = p4,
    device = fmt,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# Plot the results in a circular format
# Relevel Chromosome factor in the summary dataframe
Chr_summary$chromosome <- factor(Chr_summary$chromosome, levels = chromosome_order)

p5 <- ggplot(Chr_summary, aes(x = chromosome, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(1), width = 1, alpha = 1) +
  coord_polar(start = -pi / 16) +
  scale_y_continuous(limits = c(-1.66, 10), breaks = c(0, 3.3, 6.6, 10), expand = c(0, 0)) +
  scale_x_discrete(
    name = "Chromosome",
    labels = c("I", "II", "III", "IV", "V", "VI", 
               "VII", "VIII", "IX", "X", "XI", 
               "XII", "XIII", "XIV", "XV", "XVI")
  ) +
  geom_vline(xintercept = seq(0.5, length(unique(Chr_summary$chromosome)) + 1.5, by = 1), color = "gray", linewidth = 0.1) +
  geom_errorbar(aes(ymin = Percentage - SD, ymax = Percentage + SD, color = Category),
                width = 0.2, linewidth = 0.3,
                position = position_dodge (1)
  ) +
  scale_color_manual(values = c("> 1.2" = "#BE1823", "< 0.8" = "#2F2C7E")) +  
  scale_fill_manual(values = c("> 1.2" = "#BE1823", "< 0.8" = "#2F2C7E")) +  
  theme_minimal() +
  theme(
    panel.grid.major.y = element_line(linewidth = 0.2),
    panel.grid.minor.y = element_line(linewidth = 0.2),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_line(linewidth = 0.2),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 12)
  ) +
  labs(title = paste("Strain:", strain_name, "-", analysis_suffix),
       y = "Percentage",
       x = "Chromosome"
  )

# Save the plot
formats <- c("svg", "png")
for (fmt in formats) {
  ggsave(
    filename = file.path(out_dir, paste0("Chr_Percentages_Circular_<0.8_and_>1.2_", strain_name, "_", analysis_suffix, ".", fmt)),
    plot = p5,
    device = fmt,
    width = 8,
    height = 5,
    dpi = 300
  )
}
