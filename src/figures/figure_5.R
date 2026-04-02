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
  paste0(datadir, "drug_safety_data.csv")
)

### SETBP1 DRUG TARGETS ###
message("Plotting shared drug/SETBP1 targets plots...")
drug_targets_setbp1 %>%
  distinct(mergeTargets, n, SETBP1_target) %>%
  ggplot(aes(y = n, fill = SETBP1_target)) +
  geom_bar()

ggsave(paste0(
  filepath, "drug_targets_histogram.png"
), width = 19, height = 24, units = "cm")

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


setbp1_targets_plot <- drug_targets_setbp1 %>%
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
  )) +
  scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.45)) +
  scale_colour_viridis_d(end = 0.7) +
  labs(
    x = "-(Weighted Connectivity Score)",
    y = "SETBP1 & Drug Target",
    colour = "Drug",
    shape = "Approved Drug"
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
message("Writing full drug table with SETBP1 + safety annotations...")

full_drug_anno <- approved_drugs %>%
  mutate(SETBP1_Target = ifelse(
    pert %in% drug_targets_setbp1$pert, TRUE, FALSE
  )) %>%
  left_join(drug_safety, by = c("pref_name" = "Name"))

write_csv(full_drug_anno, paste0(filepath, "full_drug_safety_anno_table.csv"))

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
  ncol = 1
)

dev.off()

# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
