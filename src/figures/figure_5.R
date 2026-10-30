# script for plotting signature reversion drug figure
set.seed(42)

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  library(ggrepel)
  library(UpSetR)
  library(circlize)
  library(gprofiler2)
  library(cowplot)
})

# create parser object
parser <- ArgumentParser()

parser$add_argument("-i", "--input_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to the directory with input files"
)

parser$add_argument("-d", "--data_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to the data directory"
)

parser$add_argument("-o", "--output_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to the directory to write results"
)

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# start timer
ptm <- proc.time()

#### Set variables for dataset ####
indir <- args$input_dir
datadir <- args$data_dir
filepath <- args$output_dir

# create output dir if needed
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
print(paste0("Results will be saved to ", filepath))

# ggplot theme
mytheme <- theme(
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  panel.background = element_blank(),
  text = element_text(family = "Helvetica", face = "bold", size = 14)
)
#### READ IN DRUG DATA ####
# signature reversion results and annotations
approved_drugs <- read_csv(
  paste0(indir, "approved_drugs_res.csv")
)
drug_targets_setbp1 <- read_csv(
  paste0(indir, "setbp1_drug_targets.csv")
)
# DailyMed safety annotations for signature reversion drugs
drug_safety <- read_csv(
  paste0(datadir, "drug_safety_summary.csv")
)

### DAILYMED BBW/WITHDRAWN-FLAGGED DRUGS ###
bbw_drugs <- approved_drugs %>%
  left_join(drug_safety, by = c("pref_name" = "OT_Name")) %>%
  filter(N_Distinct_BBW_Texts > 0 | !is.na(Withdrawn_Description))

paste0("Withdrawn reasons include: ")
unique(bbw_drugs$Withdrawn_Description)

bbw_drugs <- unique(bbw_drugs$pref_name)
paste0("DailyMed-Flagged BBW/Withdrawn Drugs: ")
bbw_drugs

# Filter previously-approved drugs by safety (BBW/Withdrawn status)
approved_drugs <- approved_drugs %>%
  # filter for non-bbw_drugs drugs
  filter(!pref_name %in% bbw_drugs )

paste0("Total number of unique drugs passing safety filtering: ")
length(unique(approved_drugs$pref_name))

write_csv(approved_drugs, paste0(filepath, "approved_drugs_safetyfiltered.csv"))

### SETBP1 DRUG TARGETS ###
message("Checking if NA max_phase is causing 'removed 5 rows' warning...")
drug_targets_setbp1 %>%
  filter(is.na(max_phase) | is.na(WTCS) | is.na(mergeTargets)) %>%
  select(pert, pref_name, max_phase, WTCS, mergeTargets) %>%
  print()

message("Plotting shared drug/SETBP1 targets plots...")

approved_drugs_setbp1 <- drug_targets_setbp1 %>%
  filter(as.integer(max_phase) == 4) %>%
  filter(!pref_name %in% bbw_drugs)

write_csv(approved_drugs_setbp1, paste0(filepath, "setbp1_converge_targets_safetyfiltered.csv"))
 
p_targets_hist <- approved_drugs_setbp1 %>%
  distinct(mergeTargets, n, SETBP1_target) %>%
  ggplot(aes(y = n, fill = SETBP1_target)) +
  geom_bar()
 
ggsave(paste0(filepath, "drug_targets_histogram.png"),
       plot = p_targets_hist, width = 19, height = 24, units = "cm")

drug_targets_setbp1 %>%
  ggplot(aes(
    x = -WTCS,
    y = max_phase,
    colour = mergeTargets,
    shape = cell,
    size = -log(WTCS_Pval)
  )) +
  geom_point()

ggsave(paste0(
  filepath, "setbp1_drug_shared_targets_trialphase.png"
), width = 19, height = 19, units = "cm")


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
  )) +
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

ggsave(paste0(
  filepath, "setbp1_drugs_shared_phase.png"
), width = 19, height = 16, units = "cm")


setbp1_targets_plot <- approved_drugs_setbp1 %>%
  ggplot(aes(
    y = mergeTargets,
    x = -WTCS
  )) +
  geom_point(aes(
    colour = pref_name,
    size = 4,
    alpha = SETBP1_top_target,
    shape = SETBP1_top_target
  )) +
  scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.45)) +
  scale_colour_viridis_d(end = 0.7) +
  labs(
    x = "-(Weighted Connectivity Score)",
    y = "SETBP1 & Drug Target",
    colour = "Drug",
    shape = "Priority SETBP1 Target"
  ) +
  guides(size = "none", alpha = "none") +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

setbp1_targets_plot
ggsave(paste0(
  filepath, "setbp1_drugs_shared.png"
), width = 19, height = 16, units = "cm")

### APPROVED DRUGS ###
message("Plotting approved drug plots...")

approved_drugs %>%
  filter(n > 1) %>%
  ggplot(aes(
    x = -WTCS, y = reorder(pref_name, drug_top_WTCS),
    alpha = 0.7,
    colour = cell, size = -log10(WTCS_Pval)
  )) +
  geom_point()
ggsave(paste0(
  filepath, "approved_drugs_braincells_multiplecells.png"
), width = 16, height = 25, units = "cm")

approved_drugs %>%
  ggplot(aes(
    x = -WTCS, y = reorder(pref_name, -drug_top_WTCS),
    alpha = 0.7, colour = cell, size = -log10(WTCS_Pval)
  )) +
  geom_point()
ggsave(paste0(
  filepath, "approved_drugs_braincells_full.png"
), width = 16, height = 29, units = "cm")

### TOP WTCS APPROVED DRUGS, SETBP1 TARGET ANNOTATED ###
message("Plotting top WTCS-approved, SETBP1-annotated drug plots...")

top_wtcs_plot <- approved_drugs %>%
  slice_min(drug_top_WTCS, n = 50) %>%
  mutate(SETBP1_Target = ifelse(
    pert %in% drug_targets_setbp1$pert, TRUE, FALSE
  )) %>%
  ggplot(aes(
    x = -WTCS, y = reorder(pref_name, -drug_top_WTCS)
  )) +
  geom_point(aes(
    alpha = SETBP1_Target,
    colour = cell,
    size = -log10(WTCS_Pval),
    shape = SETBP1_Target
  )) +
  scale_size(range = c(4, 7.5)) +
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

top_wtcs_plot
ggsave(paste0(
  filepath, "approved_drugs_braincells_topWTCS.png"
), width = 19, height = 25, units = "cm")

### ANNOTATED DRUG TABLE ###
message("Writing drug tables with SETBP1 + safety annotations...")
converg_drug_anno <- approved_drugs %>%
  mutate(SETBP1_Target = ifelse(
    pert %in% drug_targets_setbp1$pert, TRUE, FALSE
  )) %>%
filter(SETBP1_Target == TRUE)

write_csv(
  converg_drug_anno, paste0(filepath, "setbp1_converge_drug_anno_table.csv")
)

full_drug_anno <- converg_drug_anno %>%
  left_join(drug_safety, by = c("pref_name" = "OT_Name"))

write_csv(full_drug_anno, paste0(filepath, "drug_full_anno_safetyfiltered_table.csv"))


### MOA ###
message("Plotting approved drug MOA plots...")

approved_moa <- data.frame(approved_drugs) %>%
  dplyr::select(
    pert, pref_name, WTCS, WTCS_Pval, cell, mergeMOA
  ) %>%
  separate_longer_delim(cols = mergeMOA, delim = "; ") %>%
  mutate(MOA = stringr::str_to_upper(mergeMOA)) %>%
  add_count(MOA)
str(approved_moa)
sort(table(approved_moa$MOA))
length(unique(approved_moa$mergeMOA))

approved_moa_plot <- approved_moa %>%
  filter(n > 2) %>%
  ggplot(aes(
    y = reorder(MOA, n), fill = cell
  )) +
  labs(y = "Mechanism of Action of Approved Drugs") +
  scale_fill_viridis_d(end = 0.9) +
  geom_bar(position = "stack") +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

approved_moa_plot
ggsave(paste0(
  filepath, "approved_drugs_top_moa.png"
), width = 20, height = 20, units = "cm")


#### FIGURE 5 ASSEMBLY ####
message("Plotting figure 5...")

png(paste0(filepath, "figure_5.png"),
  width = 350,
  height = 425,
  units = "mm",
  res = 300
)

cowplot::plot_grid(
  approved_moa_plot,
  top_wtcs_plot,
  setbp1_targets_plot,
  labels = "AUTO",
  rel_heights = c(1, 1.5, 1),
  ncol = 1
)

dev.off()

# Poster figure 4
message("Plotting psoter figure 4...")
png(paste0(filepath, "poster_figure_4.png"),
  width = 250,
  height = 600,
  units = "mm",
  res = 300
)

cowplot::plot_grid(
  top_wtcs_plot,
  setbp1_targets_plot,
  labels = "AUTO",
  rel_heights = c(2.5, 1),
  ncol = 1
)

dev.off()

# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
