##### Permuation Analysis of SETBP1 Targets #####
## GOAL: To determine the presence and/or directional bias of SETBP1 targets toward pre-natal or post-natal window.
## RESULTS: SETBP1 targets are robustly up-regulated during pre-natal window (i.e., early brain development)
## NOTES: Both mean and median were tested for these permutation analyses to determine 

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(gprofiler2)
  library(ggplot2)
  library(tidyverse)
  library(vroom)
})
###### SET UP ######
# create parser object
parser <- ArgumentParser()

parser$add_argument("-i", "--indir",
  type = "character",
  help = "Path to the input data directory to access expression data"
)
parser$add_argument("-o", "--outdir",
  type = "character",
  help = "Path to the directory to write results"
)
# get command line options, if help option encountered print help and exit
args <- parser$parse_args()
indir <- args$indir
outdir <- args$outdir

# create output dir if needed
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
# SETBP1 targets
setbp1_targets_table <- read_csv(paste0(indir, "260303_setbp1_targets.csv"))
setbp1_targets <- unique(setbp1_targets_table$Target)
# setbp1_targets   <- read_excel("250609_setbp1_targets_names_only.xlsx", col_names = TRUE)
# setbp1_targets <- setbp1_targets[[1]] |> unique() |> as.character()

# PAN-GO for target prioritization
pango <- vroom(paste0(indir, "pango_annotations.zip"))
cat(str(pango))

setbp1_pango <- filter(pango, gene_symbol %in% setbp1_targets)
cat(str(setbp1_pango))
terms_of_interest <- c("translation", "chromatin regulation", "mitochondria")
 


