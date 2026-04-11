#### R script to process SR coverage profiles from chr III and V, 75nt and 18nt, average all experiments
### Loop for each subdir/timepoint file in MYWD
### 08/03/2026 - Lydia

args <- commandArgs(trailingOnly = TRUE)

sample <- args[1]
outfile <- args[2]
files <- args[-c(1,2,3)]

# sample <- "T0"
# outfile <- "/Users/lydia/GWSt/WD_test_v0.3.2/1_Wt/T0_E1_75nt_avg.tsv"
# files <- "/Users/lydia/GWSt/WD_test_v0.3.2/1_Wt/T0_E1_75nt.tsv"


# Read all experiment files
tables <- lapply(files, function(f) read.table(f, sep="\t", header=FALSE, stringsAsFactors=FALSE))

# Use the coordinates from the first file (assume column 3)
coords <- tables[[1]][,3]

# Chromosome column (optional, e.g., column 1)
chromosomes <- tables[[1]][,1]

# Extract values to average (column 4)
vals <- sapply(tables, function(x) as.numeric(x[,4]))

# Compute row-wise average
avg <- rowMeans(vals)

# Build output: coordinate, sample_time, average, chromosome
out <- data.frame(
  coords,
  sample,
  avg
)

write.table(
  out,
  outfile,
  sep="\t",
  quote=FALSE,
  row.names=FALSE,
  col.names=FALSE
)
