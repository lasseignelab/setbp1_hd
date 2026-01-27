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
cachedir <- paste0(datapath, "ExperimentHubCache")
resdir <- args$output_dir
filepath <- paste0(resdir, "signature_reversion/")
querypath <- paste0(resdir, "consensus_signature/")

### TEMP ###
brain_drugs_anno <- read_rds(paste0(filepath, "gess_lincs_anno_braincell_results.rds"))
brain_drugs_filt <- brain_drugs_anno %>%
                    dplyr::select(
                       pert, pref_name, max_phase, cell, WTCS, WTCS_Pval, mergeTargets, Target_pathway) %>%
                    filter(!is.na(pref_name))
write_csv(brain_drugs_filt, paste0(filepath, "brain_drugs_filtered.csv"))


# datapath <- "./data/"
# querypath <- "./results/consensus_signature/"
# cachedir <- "/data/user/lizzyr/setbp1_hd/data/ExperimentHubCache"
# filepath <- "./results/signature_reversion/"
# create dirs if needed
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
if (!dir.exists(cachedir)) dir.create(cachedir, recursive = TRUE)

###### EXPERIMENT HUB SET UP ######

setExperimentHubOption("ASK", FALSE)
setExperimentHubOption("CACHE", cachedir)

message("Using cache directory: ", cachedir)

# getExperimentHubOption("CACHE")
# removeCache(eh, ask = FALSE)
eh <- ExperimentHub()#localHub = FALSE)
# lincs 2020 filepath and db exploration
lincs2 <- eh[["EH7297"]]

# load LINCS perturbation annotations
# data(lincs_pert_info2)

###### READ IN DATA ######

# RRA signature
rra_degs <- read_csv(paste0(querypath, "rankaggregate_consensus_signature.csv"))
setbp1_targets_all <- read_csv(paste0(
  datapath, "241203_setbp1_targets_names_only.csv"), col_names = FALSE)
setbp1_targets_top <- read_csv(paste0(
  datapath, "241203_setbp1_targets_prioritized.csv"))

###### PREPARE QUERY SIGNATURE ######

# genes available in LINCS DB
db_genes <- h5read(lincs2, "rownames", drop = TRUE)

# up vs down DEGs in signature
composite_up <- rra_degs %>%
  mutate(entrez = as.character(entrez)) %>%
  filter(meta_logFC > 0.1) %>%
  filter(entrez %in% db_genes) %>%
  slice_min(Score, n = 100)
summary(composite_up)

composite_down <- rra_degs %>%
  mutate(entrez = as.character(entrez)) %>%
  filter(meta_logFC < -0.1) %>%
  filter(entrez %in% db_genes) %>%
  slice_min(Score, n = 100)
summary(composite_down)

###### SIGNATURE REVERSION ######

# query LINCS signatures
qsig_lincs <- qSig(
  query = list(upset = composite_up$entrez, downset = composite_down$entrez), 
  gess_method = "LINCS", 
  refdb = lincs2)

# calculate signature reversion metrics of query vs LINCS DB
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
write_csv(lincs_brain, paste0(filepath, "gess_lincs_braincells_results.csv"))

# annotate drugs
brain_drugs_anno <- lincs_brain %>%
  left_join(lincs_pert_info2, by = c("pert" = "pert_id"))
head(brain_drugs_anno)
unique(brain_drugs_anno$mergeMOA)
table(brain_drugs_anno$max_phase)
length(unique(brain_drugs_anno$pert))

write_rds(brain_drugs_anno, paste0(filepath, "gess_lincs_anno_braincell_results.rds"))
# brain_drugs_anno <- read_rds(paste0(filepath, "composite_gess_lincs_anno_braincell_results.rds"))
# brain_drugs_anno<- read_rds(paste0(filepath, "gess_lincs_anno_braincell_results.rds"))
# setdiff(brain_drugs_anno$pref_name, brain_drugs_anno$pref_name)
# setdiff(brain_drugs_anno$pref_name, brain_drugs_anno$pref_name)

###### TARGET ANALYSIS ######
# separate targets by rows to compare drug targets
brain_drug_targets <- brain_drugs_anno %>%
  dplyr::select(
    pert, pref_name, max_phase, cell, WTCS, WTCS_Pval, mergeTargets, Target_pathway) %>%
  separate_longer_delim(mergeTargets, delim = "; ") %>%
  filter(!is.na(mergeTargets)) %>%
  add_count(mergeTargets) %>%
  # annotate targets that are prioritized SETBP1 regulatory targets
  mutate(
    SETBP1_target = ifelse(mergeTargets %in% setbp1_targets_all$X1, TRUE, FALSE),
    SETBP1_top_target = ifelse(mergeTargets %in% setbp1_targets_top$target, TRUE, FALSE)
  )

# drug target histogram
brain_drug_targets %>%
  distinct(mergeTargets, n, SETBP1_target) %>%
  ggplot(aes(y = n, fill = SETBP1_target)) +
  geom_bar()
ggsave(paste0(filepath, "drug_targets_histogram.png"), width = 19, height = 24, units = "cm")

# filter for drugs that target any SETBP1 regulatory targets
drug_targets_setbp1 <- brain_drug_targets %>%
  filter(SETBP1_target == TRUE) %>%
  distinct()
str(drug_targets_setbp1)
print(n = 21, arrange(drug_targets_setbp1, pref_name))

unique(drug_targets_setbp1$pref_name) 
unique(drug_targets_setbp1$mergeTargets) 
write_csv(drug_targets_setbp1, paste0(filepath, "setbp1_drug_shared_targets.csv"))

drug_targets_setbp1 %>%
  ggplot(aes(
    x = -WTCS, y = max_phase, colour = mergeTargets, shape = cell, size = -log(WTCS_Pval)
  )) +
  geom_point()

ggsave(paste0(filepath, "setbp1_drug_shared_targets_trialphase.png"), width = 19, height = 19, units = "cm")

drug_targets_setbp1 %>%
  ggplot(aes(
    y = mergeTargets, 
    x = -WTCS
   )) +
  geom_point(aes(
    colour = SETBP1_top_target, 
    size = 4, 
    alpha = 0.7,
    shape = as.factor(max_phase)
    )
  ) +
  scale_colour_viridis_d(end = 0.7) +
  labs(
    x = "-(Weighted Connectivity Score)",
    y = "SETBP1 & Drug Target",
    colour = "Priority SETBP1 Target",
    shape = "Max Drug Approval Phase"
  ) +
  guides(size = "none", alpha = "none") +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "setbp1_drugs_shared_targets.png"), width = 19, height = 16, units = "cm")

drug_targets_setbp1 %>%
  filter(!is.na(pref_name)) %>%
  ggplot(aes(
    y = mergeTargets, 
    x = -WTCS
   )) +
  geom_point(aes(
    colour = pref_name,
    size = 4, 
    alpha = SETBP1_top_target,
    shape = SETBP1_top_target
    )
  ) +
  scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.45)) +
  scale_colour_viridis_d(end = 0.7) +
  labs(
    x = "-(Weighted Connectivity Score)",
    y = "SETBP1 & Drug Target",
    colour = "Priority SETBP1 Target",
    shape = "Max Drug Approval Phase"
  ) +
  guides(size = "none", alpha = "none") +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "setbp1_drugs_shared.png"), width = 19, height = 16, units = "cm")

# full drug annotations for setbp1-target-targeting drugs

setbp1_drugs <- brain_drugs_anno %>%
  filter(pert %in% drug_targets_setbp1$pert)
summary(setbp1_drugs)
unique(setbp1_drugs$indication_class)
unique(setbp1_drugs$mechanism_of_action)

setbp1_drugs %>%
  filter(!is.na(mechanism_of_action)) %>%
  add_count(mechanism_of_action) %>%
  ggplot(aes(
    y = reorder(mechanism_of_action, n) 
    # x = -WTCS
   )) +
  geom_bar(aes(
    fill = pref_name
    # alpha = SETBP1_top_target,
    # shape = SETBP1_top_target
    )
  ) +
  # scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.45)) +
  scale_fill_viridis_d(end = 0.7) +
  labs(
    y = "MOA",
    fill = "Drug with SETBP1 Target"
  ) +
  # guides(size = "none", alpha = "none") +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "moa_setbp1_drugs_shared.png"), width = 19, height = 16, units = "cm")

# target set enrichment analysis
brain_drugs_filt <- filter(lincs_brain, !is.na(t_gn_sym))
tsea_go <- tsea_dup_hyperG(drugs = unique(brain_drugs_filt$pert), type = "GO", ont = "MF", pvalueCutoff = 0.5, minGSSize = 2, maxGSSize = 2000)
result(tsea_go)

tsea_kegg <- tsea_dup_hyperG(drugs = unique(brain_drugs_filt$pert), type = "KEGG", pvalueCutoff = 0.5, minGSSize = 2, maxGSSize = 2000)
result(tsea_kegg)

# drug set enrichment analysis
dsea_go <- dsea_hyperG(drugs = unique(lincs_brain$pert), type = "GO", ont = "MF")
result(dsea_go)

# dsea_kegg <- dsea_hyperG(drugs = unique(lincs_brain$pert), type = "KEGG")
# result(dsea_kegg)

#### APPROVED/POST-MARKET DRUGS ####
approved_drugs <- brain_drugs_anno %>%
  # max WTCS to help with plot ordering
  group_by(pref_name) %>%
  mutate(drug_top_WTCS = min(WTCS)) %>%
  ungroup() %>%
  add_count(pert) %>%
  # filter for approved drugs
  filter(max_phase == 4)
head(approved_drugs[,61:63])
unique(approved_drugs$pref_name)
write_csv(approved_drugs, paste0(filepath, "approved_drugs_braincell_res.csv"))

approved_drugs %>%
  filter(n > 1) %>%
  ggplot(aes(
    x = -WTCS, y = reorder(pref_name, drug_top_WTCS), alpha = 0.7, colour = cell, size = -log10(WTCS_Pval)
    )) +
  geom_point()
ggsave(paste0(filepath, "approved_drugs_braincells_multiplecells.png"), width = 16, height = 25, units = "cm")

approved_drugs %>%
  ggplot(aes(
    x = -WTCS, y = reorder(pref_name, -drug_top_WTCS), alpha = 0.7, colour = cell, size = -log10(WTCS_Pval)
    )) +
  geom_point()
ggsave(paste0(filepath, "approved_drugs_braincells_full.png"), width = 16, height = 29, units = "cm")

approved_drugs %>%
  slice_min(drug_top_WTCS, n = 50) %>%
  mutate(SETBP1_Target = ifelse(pert %in% drug_targets_setbp1$pert, TRUE, FALSE)) %>%
  ggplot(aes(
    x = -WTCS, y = reorder(pref_name, -drug_top_WTCS))) +
  geom_point(aes(
    alpha = SETBP1_Target, 
    colour = cell, 
    size = -log10(WTCS_Pval), 
    shape = SETBP1_Target)
) +
  scale_size(range = c(4,7.5)) +
  scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.45)) +
  scale_colour_viridis_d(end = 0.9) +
  labs(
    x = "-(Weighted Connectivity Score)",
    y = "Drug Name",
    shape = "Priority SETBP1 Target",
    size = "-log10(P-value)",
    colour = "Cell Line"
  ) +
  guides(size = "none", alpha = "none") +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "approved_drugs_braincells_topWTCS.png"), width = 19, height = 25, units = "cm")

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
  filter(n > 2) %>%
  ggplot(aes(
    y = reorder(MOA, n), fill = cell
  )) +
  labs(y = "Mechanism of Action of Approved Drugs") +
  scale_fill_viridis_d(end = 0.9) +
  geom_bar(position = "stack") +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "approved_drugs_top_moa.png"), width = 20, height = 20, units = "cm")
