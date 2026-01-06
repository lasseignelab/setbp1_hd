# lintr: disable object_usage_linter
# script for analyzing DEGs combined across contrasts
set.seed(42)

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(ggplot2)
  library(UpSetR)
  library(circlize)
  library(gprofiler2)
})

# create parser object
# parser <- ArgumentParser()
#
# parser$add_argument("-d", "--degs",
#   type = "character",
#   help = "Parent directory containing each dataset's DEG CSV files"
# )
# parser$add_argument("-o", "--output_dir",
#   type = "character",
#   metavar = "PATH",
#   help = "Path to the directory to write DEA analysis results"
# )
#
# # get command line options, if help option encountered print help and exit
# args <- parser$parse_args()
#
# # start timer
# ptm <- proc.time()
#
# #### Set variables for dataset ####
# print(paste0("Reading in DEGs from: ", args$degs))
# sig_degs_comb <- read_csv(args$degs)
# degs <- args$degs
# create output dir if needed
# filepath <- paste0(args$output_dir, "/")
degs <- "./results/dea/"
filepath <- "././results/figures/"
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

mytheme <- theme(
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  panel.background = element_blank(),
  text = element_text(family = "Helvetica", face = "bold", size = 16)
  )
#### READ IN AND COMBINE DEGS ACROSS DATA SETS ####
# helper function to read in DEGs
read_deg <- function(dataset_name) {
  file <- file.path(degs, dataset_name, "unfiltered_degs.csv")
  if (!file.exists(file)) stop("Missing file: ", file)
  read_csv(file) %>% mutate(dataset = dataset_name)
}
all_degs <- bind_rows(
  read_deg("cardo"), read_deg("shaw"), read_deg("tanaka1"), read_deg("tanaka2")
)

# DEG counts per dataset
all_degs %>%
  filter(significance != "NonSig") %>%
  add_count(dataset, significance) %>%
  distinct(dataset, significance, n) %>%
  ggplot(aes(x = n, y = dataset, fill = significance)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_viridis_d(begin = .3) +
  geom_text(
    aes(label = n, x = n + (n * 1.01), fontface = "bold"), 
    position = position_dodge(0.9)
  ) +
  scale_x_log10(expand = expansion(mult = 0, 0.4)) +
  mytheme
ggsave(filename = paste0(filepath, "degs_summary_bar.png"), width = 24, height = 22, units = "cm")

# DEG counts per contrast by dataset
all_degs %>%
  filter(significance != "NonSig" & dataset %in% c("cardo", "shaw")) %>%
  filter(!contrast %in% c("PATH_vs_PATH", "NPC_HD", "KOvsWT_NPC", "KOvsWT_all")) %>%
  add_count(contrast) %>%
  distinct(contrast, dataset, n) %>%
  ggplot(aes(x = n, y = reorder(contrast, n), fill = dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis_d(end = .7) +
  geom_text(
    aes(label = n, x = n + 350, fontface = "bold"), 
    position = position_dodge(0.9)
  ) +
  labs(x = "Number of Significant DEGs", y = "Differential Expression Contrast (Comparison)") +
  # scale_x_log10(expand = expansion(mult = 0, 0.5)) +
  mytheme
ggsave(filename = paste0(filepath, "degs_contrast_summary_bar.png"), width = 24, height = 22, units = "cm")


# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
