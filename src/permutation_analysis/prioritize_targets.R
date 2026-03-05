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

### gene functional pathway annotations ###
pathway_map <- list(
  "Chromatin_organization" = "GO:0006325",
  "Chromatin_binding" = "GO:0003682",
  "Prostaglandin_receptor_activity" = "GO:0004955",
  "Synapse_assembly" = "GO:0007416",
  "Gene_expression" = "GO:0010467",
  "Translation" = "GO:0006412",
  "Translation_regulator_activity" = "GO:0045182",
  "Mitochondrion" = "GO:0005739",
  "Mitochondrion_organization" = "GO:0007005"
)

target_pathway_ids <- unlist(pathway_map) 

pathway_df <- gconvert(
  query = target_pathway_ids,
  organism = "hsapiens",
  target = "HGNC",
  mthreshold = Inf,
  filter_na = TRUE
)
cat(str(pathway_df))

setbp1_pathways <- pathway_df %>%
  dplyr::select(input, target) %>%
  rename("STANDARD_NAME" = input, "GeneSymbol" = target) %>%
  filter(!is.na(GeneSymbol)) %>%
  filter(GeneSymbol %in% setbp1_targets)

write_csv(setbp1_pathways, paste0(indir, "setbp1_target_pathways.csv"))

print(str(setbp1_pathways))
print(unique(stringr::str_sort(setbp1_pathways$GeneSymbol)))

# terms_of_interest <- c("translation", "chromatin regulation", "mitochondria")

