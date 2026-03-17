##### Permuation Analysis Plotting #####

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
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

# Loading in permutation results
stats_results <- read_csv(paste0(outdir, "temporal_window_perm_stats.csv"))

# Plot for testing purposes
# Define a lookup table for descriptive labels and chronological order
window_metadata <- data.frame(
  window = c("Early_Prenatal", "Mid_Prenatal", "Late_Prenatal", "Infancy", "Early_Childhood", "Childhood", "Adult"),
  label = c("Early Prenatal\n(8-12 pcw)", 
            "Mid Prenatal\n(13-24 pcw)", 
            "Late Prenatal\n(25-40 pcw)", 
            "Infancy\n(0-1 yr)", 
            "Early Childhood\n(2-5 yrs)", 
            "Childhood\n(6-16 yrs)", 
            "Adult\n(18+ yrs)"),
  order = 1:7
)

# Merge and reorder
stats_results <- stats_results %>%
  left_join(window_metadata, by = "window") %>%
  mutate(label = reorder(label, order)) %>%
  filter(!is.na(label)) # Remove "Other" if it exists

# Define colors for clarity: Significant = Red, Not Significant = Grey
sig_colors <- c("TRUE" = "firebrick3", "FALSE" = "grey80")

ggplot(stats_results, aes(y = label, x = z_score)) +
  geom_segment(aes(y = label, yend = label, x = 0, xend = z_score), 
               color = "grey50", lwd = 1) +
  geom_point(aes(fill = p_val_adj < 0.05), 
             size = 6, shape = 21, color = "black") +
  geom_text(aes(label = ifelse(p_val_adj < 0.05, "*", ""),
                vjust = ifelse(z_score >= 0, -0.5, 1.5)), 
            size = 10, fontface = "bold") +
  scale_fill_manual(values = sig_colors, 
                    labels = c("FALSE" = "Not Significant", "TRUE" = "FDR < 0.05"),
                    name = "Significance") +
  scale_x_discrete(expand = expansion(add = 1.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "SETBP1 Target Enrichment Trajectory",
       subtitle = "Positive Z = Enrichment | Negative Z = Depletion",
       y = "Developmental Stage (BrainSpan)",
       x = "Enrichment Z-Score") +
  theme(
    axis.text.y = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
    legend.position = "right"
  )

ggsave(paste0(outdir, "SETBP1_temporal_trajectory.png"), width = 10, height = 7)

cat("Plot saved to SETBP1_temporal_trajectory.png")
