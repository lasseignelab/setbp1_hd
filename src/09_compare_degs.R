suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(ComplexHeatmap)
  library(statmod)
  library(gprofiler2)
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

degs <- read_csv(args$rse)
meta <- rse@metadata
sra <- args$sra

# create output dir if needed
filepath <- paste0(args$output_dir, sra, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
# save labeled and sig-filtered results
# degs_comb <- read_csv(, file = paste0(filepath, "combined_contrasts_labeled_degs.csv"))
filepath <- "./results/dea/"
cardo_degs <- read_csv(file = paste0(filepath, "cardo/combined_contrasts_sig_degs.csv"))
shaw_degs <- read_csv(file = paste0(filepath, "shaw/combined_contrasts_sig_degs.csv"))

unique(cardo_degs$contrast)
unique(shaw_degs$contrast)

# compare DEGs from NPCs across data sets
# long format
npc_degs <- cardo_degs %>%
  bind_rows(shaw_degs) %>%
  filter(contrast %in% c("KOvsWT_NPC", "NPC_HD")) %>%
  # find common DEGs -- they'll be duplicated
  add_count(gene, name = "n") %>%
  filter(n > 1)

str(npc_degs)
unique(npc_degs$gene_name)

npc_degs %>%
  ggplot(aes(x = logFC, y = gene_name, colour = contrast, size = adj.P.Val)) +
  geom_point()

ggsave(paste0(filepath, "combined/NPC_combined_degs.png"))


# wide format
npc_wide <- cardo_degs %>%
  filter(contrast == "KOvsWT_NPC") %>%
  inner_join(filter(shaw_degs, contrast == "NPC_HD"), by = "gene")

str(npc_wide)

npc_wide %>%
  ggplot(aes(x = logFC.x, y = logFC.y )) +
  geom_point()

ggsave(paste0(filepath, "combined/NPC_compare_degs.png"))


# npc_up_degs <- npc_degs %>%
  # filter(
