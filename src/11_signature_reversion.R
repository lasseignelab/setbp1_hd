suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(signatureSearch)
  library(ExperimentHub)
  library(rhdf5)
})

# create parser object
parser <- ArgumentParser()

parser$add_argument("-i", "--input_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to data directory"
)
parser$add_argument("-o", "--output_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to the results directory to read consensus signature and write signatureSearch results" # nolint
)

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# source functions file if supplied
if (!is.null(args$functions)) source(args$functions)

# start timer
ptm <- proc.time()

#### Set variables for dataset ####

datapath <- args$input_dir
resdir <- args$output_dir
filepath <- paste0(resdir, "signature_reversion/")
querypath <- paste0(resdir, "consensus_signature/")
cachedir <- paste0(datapath, "ExperimentHubCache")
message("Location of cache directory: ", cachedir)

# create dirs if needed
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
if (!dir.exists(cachedir)) dir.create(cachedir, recursive = TRUE)

###### EXPERIMENT HUB SET UP ######

setExperimentHubOption("ASK", FALSE)
setExperimentHubOption("CACHE", cachedir)

# Force a clean load of the hub
Sys.setenv(EXPERIMENT_HUB_CACHE = cachedir)
eh <- ExperimentHub(localHub = FALSE)

# lincs 2020 filepath and db exploration
message("Loading in lincs2 EH db...")
lincs2 <- eh[["EH7297"]]

# load LINCS perturbation annotations
data(lincs_pert_info2)
message("Loading in lincs2 annotation EH db...")

###### READ IN DATA ######

# RRA signature
rra_degs <- read_csv(
  paste0(querypath, "rankaggregate_consensus_signature.csv")
)
# SETBP1 regulatory targets
# Curated targets by Sasha Taluri
setbp1_targets_all <- read_csv(paste0(
  datapath, "260403_setbp1_targets.csv"
))
# Prioritized targets from coexpression & pathway analysis
setbp1_targets_top <- read_csv(paste0(
  datapath, "setbp1_targets_prioritized.csv"
))

###### PREPARE QUERY SIGNATURE ######

# genes available in LINCS DB
db_genes <- h5read(lincs2, "rownames", drop = TRUE)

# up vs down DEGs in signature
degs_up <- rra_degs %>%
  mutate(entrez = as.character(entrez)) %>%
  filter(meta_logFC > 0.1) %>%
  filter(entrez %in% db_genes) %>%
  slice_min(Score, n = 100)
cat("Up genes summary... ", summary(degs_up), "\n", sep = " ")

degs_down <- rra_degs %>%
  mutate(entrez = as.character(entrez)) %>%
  filter(meta_logFC < -0.1) %>%
  filter(entrez %in% db_genes) %>%
  slice_min(Score, n = 100)
cat("Down genes summary... ", summary(degs_down), "\n", sep = " ")

###### SIGNATURE REVERSION ######

# query LINCS signatures
qsig_lincs <- qSig(
  query = list(upset = degs_up$entrez, downset = degs_down$entrez),
  gess_method = "LINCS",
  refdb = lincs2
)

# calculate signature reversion metrics of query vs LINCS DB
lincs <- gess_lincs(qsig_lincs, sortby = "WTCS", tau = FALSE, workers = 1)
str(lincs)
write_rds(lincs, paste0(filepath, "lincs_results.rds"))

# filter for brain/CNS-derived/relevant cells (HEK cells often used in neuro)
brain_cells <- c("NEU", "SHSY5Y", "NPC", "HEK293", "HEK293T")
lincs_brain <- lincs@result %>%
  filter(cell %in% brain_cells) %>%
  filter(WTCS < 0) %>%
  arrange(WTCS)
head(lincs_brain)
write_csv(lincs_brain, paste0(filepath, "lincs_braincells_results.csv"))

# annotate drugs
brain_drugs_anno <- lincs_brain %>%
  left_join(lincs_pert_info2, by = c("pert" = "pert_id"))
cat(
  "Number of unique reversion drugs in brain-relevant cells: ",
  length(unique(brain_drugs_anno$pert)), "\n",
  sep = " "
)

write_rds(
  brain_drugs_anno, paste0(filepath, "lincs_anno_braincell_results.rds")
)

###### TARGET ANALYSIS ######
# separate targets by rows to compare drug targets
brain_drug_targets <- brain_drugs_anno %>%
  dplyr::select(
    pert, pref_name, max_phase, cell, WTCS,
    WTCS_Pval, mergeTargets, Target_pathway
  ) %>%
  separate_longer_delim(mergeTargets, delim = "; ") %>%
  filter(!is.na(mergeTargets)) %>%
  add_count(mergeTargets) %>%
  # annotate targets that are prioritized SETBP1 regulatory targets
  mutate(
    SETBP1_target = ifelse(
      mergeTargets %in% setbp1_targets_all$Target, TRUE, FALSE
    ),
    SETBP1_top_target = ifelse(
      mergeTargets %in% setbp1_targets_top$Target_Gene, TRUE, FALSE
    )
  )

# filter for drugs that target any SETBP1 regulatory targets
drug_targets_setbp1 <- brain_drug_targets %>%
  filter(SETBP1_target == TRUE) %>%
  distinct()
str(drug_targets_setbp1)
print(n = 21, arrange(drug_targets_setbp1, pref_name))

cat(
  "Drugs that target SETBP1 targets: ",
  unique(drug_targets_setbp1$pref_name),
  "\n",
  sep = " "
)
cat(
  "Drug target SETBP1 targets: ",
  unique(drug_targets_setbp1$mergeTargets),
  "\n",
  sep = " "
)
write_csv(
  drug_targets_setbp1, paste0(filepath, "setbp1_drug_targets.csv")
)

# full drug annotations for setbp1-target-targeting drugs

setbp1_drugs <- brain_drugs_anno %>%
  filter(pert %in% drug_targets_setbp1$pert)


#### APPROVED/POST-MARKET DRUGS ####
approved_drugs <- brain_drugs_anno %>%
  # max WTCS to help with plot ordering
  group_by(pref_name) %>%
  mutate(drug_top_WTCS = min(WTCS)) %>%
  ungroup() %>%
  add_count(pert) %>%
  # filter for approved drugs
  filter(max_phase == 4)

write_csv(approved_drugs, paste0(filepath, "approved_drugs_res.csv"))


# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
