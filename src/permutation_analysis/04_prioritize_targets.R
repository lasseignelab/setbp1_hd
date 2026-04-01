##### Prioritization of SETBP1 Targets #####
## GOAL: Prioritize SETBP1 targets by coexpression and pathway annotations

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
gtex_coexpr <- read_csv(paste0(outdir, "gtex_target_coexpression.csv"))
brainspan_coexpr <- read_csv(
  paste0(
    outdir, "SETBP1_target_coexpression_brainspan.csv"
  )
)
cat("\nTotal number of unique SETBP1 targets: ", length(setbp1_targets))

### TOP COEXPRESSED TARGETS ###
top_gtex_coexpr <- gtex_coexpr %>%
  filter(FDR < 0.05) %>%
  slice_max(abs(Spearman_Rho), prop = 0.10) %>%
  mutate(Method = "GTEx_Coexpr")
cat("\nTop quartile significantly coexpressed targets across tissues/ages:\n")
print(top_gtex_coexpr)

top_brainspan_coexpr <- brainspan_coexpr %>%
  filter(FDR < 0.05) %>%
  slice_max(abs(Spearman_Rho), prop = 0.10) %>%
  mutate(Method = "Brainspan_Coexpr")
cat("\nTop quartile significantly coexpressed targets across tissues/ages:\n")
print(top_brainspan_coexpr)

write_csv(
  rbind(top_gtex_coexpr, top_brainspan_coexpr),
  paste0(outdir, "top_target_coexpression.csv")
)

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
cat("\npathway_df structure:\n")
cat(str(pathway_df))

setbp1_pathways <- pathway_df %>%
  dplyr::select(input, target) %>%
  rename("Method" = input, "Target_Gene" = target) %>%
  filter(!is.na(Target_Gene)) %>%
  filter(Target_Gene %in% setbp1_targets)

write_csv(setbp1_pathways, paste0(outdir, "setbp1_target_pathways.csv"))

### COMBINED TARGET PRIORITIZATION ###
prioritized_targets <- top_gtex_coexpr %>%
  select(Target_Gene, Method) %>%
  bind_rows(setbp1_pathways)

write_csv(prioritized_targets, paste0(outdir, "setbp1_targets_prioritized.csv"))
cat("\nSaved final prioritized targets to: setbp1_targets_prioritized.csv\n\n")
cat("\nNumber of top prioritized targets: ", length(unique(prioritized_targets$Target_Gene)))
cat("\nFrequency of targets by prioritization method:\n")

print(sort(table(prioritized_targets$Target_Gene)))
