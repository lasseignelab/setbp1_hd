suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(ComplexHeatmap)
  library(statmod)
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

filepath <- "./results/dea/"
# create output dir if needed
filepath <- paste0(args$output_dir, sra, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
# save labeled and sig-filtered results
# degs_comb <- read_csv(, file = paste0(filepath, "combined_contrasts_labeled_degs.csv"))
sig_degs_comb <- read_csv(file = paste0(filepath, "shaw/combined_contrasts_sig_degs.csv"))

# volcano plots across contrasts
vol_plot <- degs_comb %>%
  # filter out highest p-adj genes for easier plotting
  #filter(adj.P.Val < 0.5) %>%
  ggplot(aes(y = -log10(`adj.P.Val`), x = logFC, color = significance)) +
  geom_point() +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
  facet_wrap(~ contrast)

vol_plot
ggsave(filename = paste0(filepath,  "degs_volcanoplot.png"))



