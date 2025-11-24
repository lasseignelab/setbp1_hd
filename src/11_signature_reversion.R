suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(ComplexHeatmap)
  library(gprofiler2)
  library(signatureSearch)
  library(ExperimentHub)
  library(rhdf5)
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

# create output dir if needed
filepath <- paste0(args$output_dir, sra, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
# save labeled and sig-filtered results
# degs_comb <- read_csv(, file = paste0(filepath, "combined_contrasts_labeled_degs.csv"))
filepath <- "./results/dea/consensus_signature/"
# save majority/high freq signature
top_freq_degs <- read_csv(paste0(filepath, "majorityfreq_consensus_signature.csv"))
# save meta analysis signature
meta_p <- read_csv(paste0(filepath, "metaanalysis_consensus_signature.csv"))
# save scores across methods for DEGs + composite scores
composite_scores <- read_csv(paste0(filepath, "full_scored_degs.csv"))
# save composite score signature
composite_signature <- read_csv(paste0(filepath, "consensus_signature.csv"))

# ENTREZ gene conversion
gene_conv <- gconvert(composite_signature$gene, target = "ENTREZGENE_ACC")
composite_conv <- composite_signature %>% inner_join(gene_conv, by = c("gene" = "input"))
# up vs down DEGs in signature
composite_up <- composite_conv %>% filter(meta_logFC > 0)
composite_up
composite_down <- composite_conv %>% filter(meta_logFC < 0)
composite_down

eh <- ExperimentHub()
lincs2 <- eh[["EH7297"]]

# 83/114 up and 44/81 down genes shared in ref db!!
qsig_lincs <- qSig(query = list(upset = composite_up$target, downset = composite_down$target), gess_method = "LINCS", refdb = lincs2)

lincs <- gess_lincs(qsig_lincs, sortby="NCS", tau=FALSE, workers=1)
str(lincs)
unique(lincs@result$cell)
lincs_brain <- lincs@result %>%
  filter(cell %in% c("NEU", "NPC")) %>%
  filter(WTCS < 0) %>%
  arrange(WTCS)
head(lincs_brain)

# rm BRD compounds
brd_comps <- grep("BRD-", lincs_brain$pert, value = TRUE)
filtered_brain <- lincs_brain %>%
  filter(!pert %in% brd_comps)
filtered_brain
