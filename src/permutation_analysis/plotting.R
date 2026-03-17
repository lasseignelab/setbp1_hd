##### Permuation Analysis Plotting #####

suppressPackageStartupMessages({
  library(argparse)
  library(readr)
  library(ggplot2)
  library(tidyverse)
  library(vroom)
  library(viridis)
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
gtex_results <- read_csv(paste0(outdir, "postnatal_gtex_permutation_results.csv"))

# Set up vs down colors
up_color <- viridis(20, option = "magma")[[17]]
dn_color <-  viridis(20, option = "magma")[[3]]

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
  order = 7:1
)

# Merge and reorder
stats_results <- stats_results %>%
  left_join(window_metadata, by = "window") %>%
  mutate(label = reorder(label, order)) %>%
  filter(!is.na(label)) # Remove "Other" if it exists

# Define colors for clarity: Significant = Red, Not Significant = Grey
sig_alpha <- c("TRUE" = 0.9, "FALSE" = 0.2)
enrich_colors <- c("TRUE" = up_color, "FALSE" = dn_color)

ggplot(stats_results, aes(y = label, x = z_score)) +
  geom_segment(aes(y = label, yend = label, x = 0, xend = z_score), 
               color = "grey50", lwd = 1) +
  geom_point(aes(fill = z_score > 0, alpha = p_val_adj < 0.05), 
             size = 6, shape = 21, color = "black") +
  geom_text(aes(label = ifelse(p_val_adj < 0.05, "*", ""),
                hjust = ifelse(z_score >= 0, -1.0, 2.0)), 
            vjust = 0.75, size = 10, fontface = "bold") +
  scale_alpha_manual(values = sig_alpha, 
                    labels = c("FALSE" = "Not Significant", "TRUE" = "FDR < 0.05"),
                    name = "Significance") +
  scale_fill_manual(values = enrich_colors, 
                    labels = c("FALSE" = "Depletion", "TRUE" = "Enrichment"),
                    name = "Z-Score Direction") +
  scale_x_continuous(expand = expansion(add = 1.5)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "SETBP1 Target Enrichment Trajectory",
       y = "Developmental Stage (BrainSpan)",
       x = "Enrichment Z-Score") +
  theme(
    axis.text = element_text(size = 10, face = "bold", hjust = 1),
    legend.position = "right"
  )

ggsave(paste0(outdir, "SETBP1_temporal_trajectory.png"), width = 10, height = 7)

cat("Plot saved to SETBP1_temporal_trajectory.png")

## HEATMAP ##
# Order the Age column chronologically so the x-axis flows naturally
gtex_results$Age <- factor(gtex_results$Age, 
                            levels = c("Pediatric", "20-29", "30-39", "40-49", "50-59", "60-69", "70-79", "Unknown"))

# Filter out 'Unknown' ages
plot_data <- subset(gtex_results, Age != "Unknown")

# Generate the Heatmap
ggplot(plot_data, aes(x = Age, y = Tissue, fill = Z_Score)) +
  geom_tile(color = "white", size = 0.5) +
  # Add significance stars for FDR < 0.05
  geom_text(aes(label = ifelse(FDR < 0.05, "*", "")), 
            color = "black", size = 6, vjust = 0.7) +
  # Diverging color palette: Blue (depleted) -> White (0) -> Red (enriched)
  scale_fill_gradient2(low = dn_color, mid = "white", high = up_color, midpoint = 0, 
                       name = "Z-Score\n(Enrichment)") +
  theme_minimal(base_size = 14) +
  labs(title = "SETBP1 Target Expression Across Tissues and Development",
       x = "Age Bracket",
       y = "Tissue") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank() # Removes grid lines for a cleaner tile look
  )

ggsave(paste0(outdir, "SETBP1_target_exp_gtex_heatmap.png"))

cat("GTEx postnatal heatmap saved. Finished!")
