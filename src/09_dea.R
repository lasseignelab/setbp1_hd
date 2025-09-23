suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(ComplexHeatmap)
  library(statmod)
  library(UpSetR)
  # source(paste0(getwd(), "/functions.R"))
  # source(paste0(getwd(), "/src/functions.R"))
})

# create parser object
parser <- ArgumentParser()

parser$add_argument("-d", "--degs",
  type = "character",
  metavar = "PATH",
  help = "RSE file (rds format) for EDA from recount3 data (00_data_download or 01_data_preprocessing)" # nolint
)
parser$add_argument("-s", "--sra",
  type = "character",
  metavar = "PATH",
  help = "SRA or dataset for naming purposes" # nolint
)
parser$add_argument("-o", "--output_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to the directory to create and write EDA results"
)

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# source functions file if supplied
if (!is.null(args$functions)) source(args$functions)

# start timer
ptm <- proc.time()

#### Set variables for dataset ####
print(paste0("Reading in RSE from: ", args$rse))

rse <- readRDS(args$rse)
meta <- rse@metadata
sra <- args$sra
# create output dir if needed
filepath <- paste0(args$output_dir, sra, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

filepath <- "./results/dea/"

# degs_comb <- read_csv(, file = paste0(filepath, "combined_contrasts_labeled_degs.csv"))
sig_degs_comb <- read_csv(file = paste0(filepath, "shaw/combined_contrasts_sig_degs.csv"))


gene_contrast_matrix <- sig_degs_comb %>%
  select(gene, contrast) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = contrast, values_from = value, values_fill = 0)

# Convert to data frame for UpSetR
gene_contrast_matrix <- as.data.frame(gene_contrast_matrix)
rownames(gene_contrast_matrix) <- gene_contrast_matrix$gene
gene_contrast_matrix$gene <- NULL
str(gene_contrast_matrix)

# Plot using UpSetR
png(paste0(filepath, "upset_plot.png"), width = 30, height = 20, units = "cm", res = 300)
upset(gene_contrast_matrix, sets = colnames(gene_contrast_matrix), order.by = "freq")
dev.off()
