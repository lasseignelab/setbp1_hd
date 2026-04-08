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
setbp1_targets_table <- read_csv(paste0(indir, "260403_setbp1_targets.csv"))
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
  slice_max(abs(Spearman_Rho), prop = 0.25) %>%
  mutate(Method = "GTEx_Coexpr")
cat("\nTop quartile significantly coexpressed targets across tissues/ages:\n")
print(head(top_gtex_coexpr))

top_brainspan_coexpr <- brainspan_coexpr %>%
  filter(FDR < 0.05) %>%
  slice_max(abs(Spearman_Rho), prop = 0.25) %>%
  mutate(Method = "Brainspan_Coexpr")
cat("\nTop quartile significantly coexpressed targets across tissues/ages:\n")
print(head(top_brainspan_coexpr))

top_coexpr <- rbind(top_gtex_coexpr, top_brainspan_coexpr)

write_csv(
  top_coexpr, paste0(outdir, "prioritized_target_coexpression.csv")
)

### gene functional pathway annotations ###
pathway_map <- list(
  # Transcriptional and epigenetic regulation
  "Gene_expression" = "GO:0010467",
  "Positive_regulation_RNA_metabolic_process" = "GO:0051254",
  "Chromatin_organization" = "GO:0006325",
  "Chromatin_binding" = "GO:0003682",
  "Translation_regulator_activity" = "GO:0045182",
  # Cell proliferation and survival
  "Negative regulation of protein phosphatase activity" = "GO:0010923",
  "Regulation of protein kinase B signaling" = "GO:0051896",
  "Intracellular signal transduction" = "GO:0035556",
  "Regulation of cell cycle process" = "GO:0010564",
  "Transforming growth factor beta receptor signaling pathway" = "GO:0007179",
  # Neurodevelopment and circuitry
  "Regulation of neural precursor cell proliferation" = "GO:2000177",
  "Synapse_assembly" = "GO:0007416",
  "Prostaglandin_receptor_activity" = "GO:0004955"
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

write_csv(setbp1_pathways, paste0(outdir, "prioritized_target_pathways.csv"))

### COMBINED TARGET PRIORITIZATION ###
# Targets filtered by PHAROS development leve
pharos_druggable <- setbp1_targets_table %>%
  filter(`Target Development Level` %in% c("Tchem", "Tclin"))

# Convert pathway map to df
pathway_map_df <- data.frame(
  "Method" = target_pathway_ids,
  "GO_Term" = names(pathway_map)
)
str(pathway_map_df)

prioritized_targets <- top_coexpr %>%
  select(Target_Gene, Method, Spearman_Rho) %>%
  bind_rows(setbp1_pathways) %>%
  left_join(pathway_map_df, by = c("Method")) %>%
  group_by(Target_Gene) %>%
  summarize(
    across(
      everything(), ~ paste(na.omit(.), collapse = "; ")
    ),
    .groups = "drop"
  ) %>%
  filter(Target_Gene %in% pharos_druggable$Target) %>%
  mutate(GO_Term = gsub("_", " ", GO_Term))

write_csv(prioritized_targets, paste0(outdir, "prioritized_setbp1_targets.csv"))
cat("\nSaved final prioritized targets to: prioritized_setbp1_targets.csv\n\n")
cat(
  "\nNumber of top prioritized targets: ",
  length(unique(prioritized_targets$Target_Gene))
)
cat("\nStructure of saved table:\n")
print(str(prioritized_targets))

#### END ####

# end timer
fptm <- proc.time()
cat("\n\nTotal runtime:  ", (fptm[3] / 60), "\n")

#  session info
cat("\n\nSession Info:\n")
print(sessionInfo())
