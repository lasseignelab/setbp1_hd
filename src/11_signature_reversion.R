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
#
# # create parser object
# parser <- ArgumentParser()
#
# parser$add_argument("-d", "--degs",
#   type = "character",
#   metavar = "PATH",
#   help = "RSE file (rds format) for EDA from recount3 data (00_data_download or 01_data_preprocessing)" # nolint
# )
# parser$add_argument("-s", "--sra",
#   type = "character",
#   metavar = "PATH",
#   help = "SRA or dataset for naming purposes" # nolint
# )
# parser$add_argument("-o", "--output_dir",
#   type = "character",
#   metavar = "PATH",
#   help = "Path to the directory to create and write EDA results"
# )
#
# # get command line options, if help option encountered print help and exit
# args <- parser$parse_args()
#
# # source functions file if supplied
# if (!is.null(args$functions)) source(args$functions)
#
# # start timer
# ptm <- proc.time()
#
# #### Set variables for dataset ####
# print(paste0("Reading in RSE from: ", args$rse))
#
# degs <- read_csv(args$rse)

cachedir <- "/data/user/lizzyr/setbp1_hd/data/ExperimentHubCache"
filepath <- "./results/signature_reversion/composite_"

# create dirs if needed
# filepath <- paste0(args$output_dir, sra, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
if (!dir.exists(cachedir)) dir.create(cachedir, recursive = TRUE)

# majority/high freq signature
datapath <- "./results/consensus_signature/"
top_freq_degs <- read_csv(paste0(datapath, "majorityfreq_consensus_signature.csv"))
# RRA signature
rra_degs <- read_csv(paste0(datapath, "rankaggregate_consensus_signature.csv"))
# meta analysis signature
meta_p <- read_csv(paste0(datapath, "metaanalysis_consensus_signature.csv"))
# scores across methods for DEGs + composite scores
# full_scores <- read_csv(paste0(datapath, "rra/compositescore_full_table_degs.csv"))
setbp1_targets <- read_csv("./doc/241203_setbp1_targets_prioritized.csv")#, col_names = FALSE)
# save composite score signature
# TODO: FIX FILEPATH...
composite_signature <- read_csv(paste0(datapath, "compositescore_consensus_signature.csv"))

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

data(lincs_pert_info2)

db_genes <- h5read(lincs2, "rownames", drop = TRUE)

# up vs down DEGs in signature
# length(unique(composite_scores$gene))
# composite_signature <- composite_scores %>%
#   mutate(entrez = as.character(entrez))# %>%
  # slice_max(composite_score, n = 200)
# composite_signature <- meta_p
# composite_signature <- full_scores %>% filter(gene %in% rra_degs$gene)
composite_up <- rra_degs %>% #composite_signature %>%
  mutate(entrez = as.character(entrez)) %>%
  filter(entrez %in% db_genes) %>%
  slice_max(meta_logFC, n = 100)
    # meta_logFC > 0)
composite_up
composite_down <- rra_degs %>% #composite_signature %>% 
  mutate(entrez = as.character(entrez)) %>%
  filter(entrez %in% db_genes) %>%
  slice_min(meta_logFC, n = 100)
  # filter(meta_logFC < 0)
composite_down

# composite: 112/151 up, 103/149 down meta-p: 95/145 up, 88/155 down
# rra: 94/121 up, 77/106 down
qsig_lincs <- qSig(query = list(upset = composite_up$entrez, downset = composite_down$entrez), gess_method = "LINCS", refdb = lincs2)

lincs <- gess_lincs(qsig_lincs, sortby="WTCS", tau=FALSE, workers=1)
str(lincs)
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
# lincs_brain_drugs <- read_rds(paste0(filepath, "gess_lincs_anno_braincell_results.rds"))

#### TARGET ANALYSIS ####
# separate targets by rows to compare drug targets
brain_drug_targets <- lincs_brain_drugs %>%
  dplyr::select(
    pert, pref_name, max_phase, cell, WTCS, WTCS_Pval, mergeTargets, Target_pathway) %>%
  separate_longer_delim(mergeTargets, delim = "; ") %>%
  filter(!is.na(mergeTargets)) %>%
  add_count(mergeTargets) %>%
  mutate(
    SETBP1_target = ifelse(mergeTargets %in% setbp1_targets$target, TRUE, FALSE))

# drug target histogram
brain_drug_targets %>%
  distinct(mergeTargets, n, SETBP1_target) %>%
  # filter(n > 4) %>%
  # ggplot(aes(y = reorder(mergeTargets, n), x = n)) +
  ggplot(aes(y = n, fill = SETBP1_target)) +
  geom_bar()
  # geom_bar(stat = "identity", aes(fill = SETBP1_target))
  # geom_histogram(stat = "count")
ggsave(paste0(filepath, "drug_targets_histogram.png"), width = 19, height = 24, units = "cm")

# find SETBP1 targets that are also drug targets
drug_targets_setbp1 <- brain_drug_targets %>%
  inner_join(setbp1_targets, by = c("mergeTargets" = "target")) %>%
  # semi_join(brain_drug_targets %>% filter(SETBP1_target == TRUE), by = "pert") %>%
  # separate_longer_delim(mergeTargets, delim = "; ") %>%
  distinct()

unique(drug_targets_setbp1$pref_name) 
write_csv(drug_targets_setbp1, paste0(filepath, "setbp1_drug_shared_targets.csv"))

drug_targets_setbp1 %>%
  ggplot(aes(
    x = -WTCS, y = max_phase, colour = mergeTargets, shape = cell, size = -log(WTCS_Pval)
  )) +
  geom_point()

ggsave(paste0(filepath, "setbp1_drug_shared_targets_trialphase.png"), width = 19, height = 19, units = "cm")

drug_targets_setbp1 %>%
  ggplot(aes(y = mergeTargets, x = -WTCS, colour = pref_name, size = 4, shape = as.factor(max_phase))) +
  geom_point()

ggsave(paste0(filepath, "setbp1_drugs_shared_targets.png"), width = 19, height = 16, units = "cm")
# target set enrichment analysis
brain_drugs_filt <- filter(lincs_brain, !is.na(t_gn_sym))
tsea_go <- tsea_dup_hyperG(drugs = unique(brain_drugs_filt$pert), type = "GO", ont = "MF", pvalueCutoff = 0.5, minGSSize = 2, maxGSSize = 2000)
result(tsea_go)

tsea_kegg <- tsea_dup_hyperG(drugs = unique(brain_drugs_filt$pert), type = "KEGG", pvalueCutoff = 0.5, minGSSize = 2, maxGSSize = 2000)

# drug set enrichment analysis
dsea_go <- dsea_hyperG(drugs = unique(lincs_brain$pert), type = "GO", ont = "MF")
result(dsea_go)

# dsea_kegg <- dsea_hyperG(drugs = unique(lincs_brain$pert), type = "KEGG")
# result(dsea_kegg)

#### APPROVED/POST-MARKET DRUGS ####
approved_drugs <- lincs_brain_drugs %>%
  # filter for approved drugs
  filter(max_phase == 4)
unique(approved_drugs$pref_name)
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
    y = reorder(MOA, n), fill = cell
  )) +
  geom_bar(position = "stack")
ggsave(paste0(filepath, "approved_drugs_top_moa.png"), width = 20, height = 20, units = "cm")

# # compare sig results across methods
# cs <- read_rds(paste0(filepath, "gess_lincs_anno_braincell_results.rds"))
# rra <- read_rds(paste0(filepath, "rra_gess_lincs_anno_braincell_results.rds"))
# ma <- read_rds(paste0(filepath, "metap_gess_lincs_anno_braincell_results.rds"))
# length(unique(cs$pert)) # 870
# length(unique(rra$pert)) # 885
# length(unique(ma$pert)) # 1063
# length(intersect(cs$pert, rra$pert)) # 381 # now 542
# length(intersect(cs$pert, ma$pert)) # 288 # now 339
# length(intersect(rra$pert, ma$pert)) # 299
# # drugs in all methods
# consistent_drugs <- cs %>%
#   filter(pert %in% rra$pert) %>%
#   filter(pert %in% ma$pert)
# unique(consistent_drugs$pref_name)
# write_csv(consistent_drugs, paste0(filepath, "intersect_drugs_composite_rra_ma251201.csv"))
# prev_consis <- read_csv(paste0(filepath, "intersect_drugs_composite_rra_ma.csv"))
# unique(prev_consis$pref_name)
# setdiff(prev_consis$pref_name, consistent_drugs$pref_name)
