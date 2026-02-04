# lintr: disable object_usage_linter
# script for analyzing DEGs combined across contrasts
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
parser$add_argument("-d", "--degs",
  type = "character",
  help = "Parent directory containing each dataset's DEG CSV files"
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
print(paste0("Reading in DEGs from: ", args$degs))
indir <- args$input_dir
degs <- args$degs
outdir <- args$output_dir

# create output dir if needed
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

mytheme <- theme(
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  panel.background = element_blank(),
  text = element_text(family = "Helvetica", face = "bold", size = 14)
)
#### READ IN AND COMBINE DEGS ACROSS DATA SETS ####
# helper function to read in DEGs
read_deg <- function(dataset_name) {
  file <- file.path(degs, dataset_name, "unfiltered_degs.csv")
  if (!file.exists(file)) stop("Missing file: ", file)
  read_csv(file) %>% mutate(dataset = dataset_name)
}
all_degs <- bind_rows(
  read_deg("cardo"),
  read_deg("shaw"),
  read_deg("tanaka1"),
  read_deg("tanaka2")
)


#### INITIAL CODE TRANSFER FROM OTHER SCRIPT ####
# Read in data
brain_drugs_anno <- read_rds(
  paste0(filepath, "gess_lincs_anno_braincell_results.rds")
)
brain_drugs_filt <- brain_drugs_anno %>%
  dplyr::select(
    pert, pref_name, max_phase, cell, WTCS,
    WTCS_Pval, mergeTargets, Target_pathway
  ) %>%
  filter(!is.na(pref_name))
write_csv(brain_drugs_filt, paste0(filepath, "brain_drugs_filtered.csv"))


# drug target histogram
brain_drug_targets %>%
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
  filepath, "setbp1_drugs_shared_targets.png"
), width = 19, height = 16, units = "cm")


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
  )) +
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

ggsave(paste0(
  filepath, "setbp1_drugs_shared.png"
), width = 19, height = 16, units = "cm")


setbp1_drugs %>%
  filter(!is.na(mechanism_of_action)) %>%
  add_count(mechanism_of_action) %>%
  ggplot(aes(
    y = reorder(mechanism_of_action, n)
  )) +
  geom_bar(aes(
    fill = pref_name
  )) +
  scale_fill_viridis_d(end = 0.7) +
  labs(
    y = "MOA",
    fill = "Drug with SETBP1 Target"
  ) +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(
  filepath, "moa_setbp1_drugs_shared.png"
), width = 19, height = 16, units = "cm")


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

approved_drugs %>%
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

ggsave(paste0(
  filepath, "approved_drugs_braincells_topWTCS.png"
), width = 19, height = 25, units = "cm")

### MOA ###
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

ggsave(paste0(
  filepath, "approved_drugs_top_moa.png"
), width = 20, height = 20, units = "cm")


#### FIGURE 3 ASSEMBLY ####
bottom_row <- plot_grid(
  performance_plot, volcano_plot,
  labels = c("B", "C")
)
png(paste0(outdir, "figure_3.png"),
  width = 350,
  height = 425,
  units = "mm",
  res = 300
)

cowplot::plot_grid(
  ora_plot,
  bottom_row,
  labels = c("A", ""),
  ncol = 1
)

dev.off()

# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
