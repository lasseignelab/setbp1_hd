suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(ComplexHeatmap)
  library(statmod)
  library(gprofiler2)
  # source(paste0(getwd(), "/functions.R"))
  # source(paste0(getwd(), "/src/functions.R"))
})

# create parser object
parser <- ArgumentParser()

parser$add_argument("-r", "--rse",
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

rse <- readRDS("./data/processed/cardo_prepped_rse.rds")
sra <- "cardo"
filepath <- "./results/archive/"
# save labeled and sig-filtered results
# write_csv(degs_comb, file = paste0(filepath, "combined_contrasts_labeled_degs.csv"))
# write_csv(sig_degs_comb, file = paste0(filepath, "combined_contrasts_sig_degs.csv"))
library(gprofiler2)
day15_degs <- sig_degs_comb %>%
  filter(contrast == "KOvsWT_15" & significance == "Up") %>%
  # filter(contrast == "KOvsWT_15" & significance == "Down") %>%
  slice_min(n = 1000, order_by = adj.P.Val)

day15_kegg <- gost(query = day15_degs$gene, ordered_query = FALSE, sources = "KEGG")
day15_kegg$result$term_name

day15_kegg$result %>%
  slice_max(n = 30, order_by = precision) %>%
  ggplot(aes(x = precision, y = term_name, size = p_value, colour = intersection_size)) +
  geom_point()
ggsave(filename = paste0(filepath,  "day15_up_kegg_enrichment.png"))
# ggsave(filename = paste0(filepath,  "day15_down_kegg_enrichment.png"))

# KEGG for day 21
day21_degs <- sig_degs_comb %>%
  filter(contrast == "KOvsWT_21" & significance == "Up") %>%
  slice_min(n = 1000, order_by = adj.P.Val)

day21_kegg <- gost(query = day21_degs$gene, ordered_query = FALSE, sources = "KEGG")
day21_kegg$result$term_name

day21_kegg$result %>%
  slice_max(n = 30, order_by = precision) %>%
  ggplot(aes(x = precision, y = term_name, size = p_value, colour = intersection_size)) +
  geom_point()
ggsave(filename = paste0(filepath,  "day21_up_kegg_enrichment.png"))

# KEGG for day 34 
day34_degs <- sig_degs_comb %>%
  filter(contrast == "KOvsWT_34" & significance == "Up") %>%
  slice_min(n = 1000, order_by = adj.P.Val)

day34_kegg <- gost(query = day34_degs$gene, ordered_query = FALSE, sources = "KEGG")
day34_kegg$result$term_name

day34_kegg$result %>%
  slice_max(n = 30, order_by = precision) %>%
  ggplot(aes(x = precision, y = term_name, size = p_value, colour = intersection_size)) +
  geom_point()
ggsave(filename = paste0(filepath,  "day34_up_kegg_enrichment.png"))

