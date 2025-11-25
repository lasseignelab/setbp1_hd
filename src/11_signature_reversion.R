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
# filepath <- paste0(args$output_dir, sra, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
if (!dir.exists(cachedir)) dir.create(cachedir, recursive = TRUE)
# save labeled and sig-filtered results
# degs_comb <- read_csv(, file = paste0(filepath, "combined_contrasts_labeled_degs.csv"))
cachedir <- "/data/user/lizzyr/setbp1_hd/data/ExperimentHubCache"
filepath <- "./results/dea/signature_reversion/"

# # avoid lock file issues
# Sys.setenv(BIOCFILECACHE_NOLOCK = "TRUE")
# h5disableFileLocking()
# Sys.setenv(HDF5_USE_FILE_LOCKING = "FALSE")
#
# save majority/high freq signature
top_freq_degs <- read_csv(paste0(filepath, "majorityfreq_consensus_signature.csv"))
# save meta analysis signature
meta_p <- read_csv(paste0(filepath, "metaanalysis_consensus_signature.csv"))
# save scores across methods for DEGs + composite scores
composite_scores <- read_csv(paste0(filepath, "full_scored_degs.csv"))
setbp1_targets <- read_csv("./data/241203_setbp1_targets_names_only.csv", col_names = FALSE)
# save composite score signature
# composite_signature <- read_csv(paste0(filepath, "consensus_signature.csv"))

# up vs down DEGs in signature
composite_signature <- composite_scores %>%
  mutate(entrez = as.character(entrez))# %>%
  # slice_max(composite_score, n = 200)
composite_up <- composite_signature %>% filter(meta_logFC > 0)
composite_up
composite_down <- composite_signature %>% filter(meta_logFC < 0)
composite_down

# explore sigsearch data
setExperimentHubOption("ASK", FALSE)
setExperimentHubOption("CACHE", cachedir)

message("Using cache directory: ", cachedir)

# getExperimentHubOption("CACHE")
# removeCache(eh, ask = FALSE)
eh <- ExperimentHub(localHub = FALSE)
ssdata <- query(eh, c("signatureSearchData"))
ssdata
ssdata$title
as.list(ssdata[4])

# lincs 2020 filepath and db exploration
lincs2 <- eh[["EH7297"]]
rhdf5::h5ls(lincs2)
# lincs_db <- rhdf5::h5read(lincs2, "assay")
# rowdata <- rhdf5::h5read(lincs2, "rownames")
# coldata <- rhdf5::h5read(lincs2, "colData")
# str(coldata)
# dim(lincs_db)
#
data(lincs_pert_info2)

# 83/114 up and 44/81 down genes shared in ref db!!
# when pre-filtering for ENTREz genes... 70/105 up, 41/95 down???
# w/o top 200 selection... 111/187 up, 79/180 down
qsig_lincs <- qSig(query = list(upset = composite_up$entrez, downset = composite_down$entrez), gess_method = "LINCS", refdb = lincs2)

lincs <- gess_lincs(qsig_lincs, sortby="WTCS", tau=FALSE, workers=1)
str(lincs)
unique(lincs@result$cell)
write_rds(lincs, paste0(filepath, "gess_lincs_results.rds"))
# lincs <- read_rds(paste0(filepath, "gess_lincs_results.rds"))

# filter for brain/CNS-derived/relevant cells (HEK cells often used in neuroscience)
brain_cells <- c("NEU", "SHSY5Y", "NPC", "HEK293", "HEK293T")
lincs_brain <- lincs@result %>%
  filter(cell %in% brain_cells) %>%
  filter(WTCS < 0) %>%
  arrange(WTCS)
head(lincs_brain)

# annotate drugs
lincs_brain_drugs <- lincs_brain %>%
  left_join(lincs_pert_info2, by = c("pert" = "pert_id"))
head(lincs_brain_drugs)
unique(lincs_brain_drugs$mergeMOA)
table(lincs_brain_drugs$max_phase)

write_rds(lincs_brain_drugs, paste0(filepath, "gess_lincs_anno_braincell_results.rds"))

#### TARGET ANALYSIS ####
# separate targets by rows to compare drug targets
brain_drug_targets <- lincs_brain_drugs %>%
  dplyr::select(
    pert, pref_name, max_phase, cell, WTCS, WTCS_Pval, mergeTargets, Target_pathway) %>%
  separate_longer_delim(mergeTargets, delim = "; ")
brain_drug_targets

# find SETBP1 targets that are also drug targets
brain_drug_targets_setbp1 <- brain_drug_targets %>%
  filter(mergeTargets %in% setbp1_targets$X1)

brain_drug_targets_setbp1 
write_csv(brain_drug_targets_setbp1, paste0(filepath, "setbp1_drug_shared_targets.csv"))

brain_drug_targets_setbp1 %>%
  ggplot(aes(
    x = -WTCS, y = max_phase, colour = mergeTargets, shape = cell, size = -log(WTCS_Pval)
  )) +
  geom_point()

ggsave(paste0(filepath, "setbp1_drug_shared_targets_trialphase.png"), width = 19, height = 19, units = "cm")

# target set enrichment analysis
brain_drugs_filt <- filter(lincs_brain, !is.na(t_gn_sym))
tsea_go <- tsea_dup_hyperG(drugs = unique(brain_drugs_filt$pert), type = "GO", ont = "MF", pvalueCutoff = 0.5, minGSSize = 2, maxGSSize = 2000)
result(tsea_go)

tsea_kegg <- tsea_dup_hyperG(drugs = unique(brain_drugs_filt$pert), type = "KEGG", pvalueCutoff = 0.5, minGSSize = 2, maxGSSize = 2000)
result(tsea_kegg)

# drug set enrichment analysis
dsea_go <- dsea_hyperG(drugs = unique(lincs_brain$pert), type = "GO", ont = "MF")
result(dsea_go)

dsea_kegg <- dsea_hyperG(drugs = unique(lincs_brain$pert), type = "KEGG")
result(dsea_kegg)

#### APPROVED/POST-MARKET DRUGS ####
approved_drugs <- lincs_brain_drugs %>%
  # filter for approved drugs
  filter(max_phase == 4)
sort(table(approved_drugs$pref_name))
write_csv(approved_drugs, paste0(filepath, "approved_drugs_braincell_res.csv"))

approved_drugs %>%
  # max WTCS to help with plot ordering
  group_by(pref_name) %>%
  mutate(drug_max_WTCS = max(abs(WTCS))) %>%
  ggplot(aes(
    x = -WTCS, y = reorder(pref_name, drug_max_WTCS), alpha = 0.7, colour = cell, size = -log10(WTCS_Pval)
    )) +
  geom_point()
ggsave(paste0(filepath, "approved_drugs_braincells_full.png"), width = 16, height = 25, units = "cm")

approved_drugs %>%
  slice_min(WTCS, n = 30) %>%
  ggplot(aes(
    x = -WTCS, y = reorder(pref_name, -WTCS), alpha = 0.7, colour = cell, size = -log10(WTCS_Pval)
    )) +
  geom_point()
ggsave(paste0(filepath, "approved_drugs_braincells_topWTCS.png"), width = 16, height = 25, units = "cm")

### MOA ###
approved_moa <- data.frame(approved_drugs) %>%
  dplyr::select(pert, pref_name, WTCS, WTCS_Pval, cell, mergeMOA) %>%
  separate_longer_delim(cols = mergeMOA, delim = "; ") %>%
  mutate(MOA = stringr::str_to_upper(mergeMOA)) %>%
  add_count(MOA)
str(approved_moa)
sort(table(approved_moa$MOA))
length(unique(approved_moa$mergeMOA))

approved_moa %>%
  filter(n > 1) %>%
  ggplot(aes(
    y = MOA, fill = cell
  )) +
  geom_bar(position = "fill")
ggsave(paste0(filepath, "approved_drugs_top_moa.png"), width = 16, height = 25, units = "cm")

