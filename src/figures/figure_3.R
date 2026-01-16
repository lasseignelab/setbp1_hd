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
 outdir <- "././results/figures/"
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
  read_deg("cardo"), read_deg("shaw"), read_deg("tanaka1"), read_deg("tanaka2")
)
#### READ IN ORA, PERFORMANCE, AND RRA-META RESULTS ####
# consensus DEGs
consensus_long <- read_csv(paste0(indir, "consensus_degs_unfiltered_long.csv"))
# ORA results for signatures
results_df <- read_csv(paste0(indir, "ORA_signatures_results.csv"))
# performance metrics for signatures
performance_res <- read_csv(paste0(indir, "signature_performance_metrics.csv"))
# RRA signature with meta-analysis logFCs
rra_meta <- read_csv(paste0(indir, "RRA_meta_DEGs.csv"))

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
ggsave(
  filename = paste0(
    outdir, "degs_summary_bar.png"
  ), width = 24, height = 22, units = "cm"
)

# DEG counts per contrast by dataset
all_degs %>%
  filter(significance != "NonSig" & dataset %in% c("cardo", "shaw")) %>%
  filter(
    !contrast %in% c("PATH_vs_PATH", "NPC_HD", "KOvsWT_NPC", "KOvsWT_all")
  ) %>%
  add_count(contrast) %>%
  distinct(contrast, dataset, n) %>%
  ggplot(aes(x = n, y = reorder(contrast, n), fill = dataset)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_viridis_d(end = .7) +
  geom_text(
    aes(label = n, x = n + 350, fontface = "bold"),
    position = position_dodge(0.9)
  ) +
  labs(
    x = "Number of Significant DEGs",
    y = "Differential Expression Contrast (Comparison)"
  ) +
  mytheme
ggsave(
  filename = paste0(outdir, "degs_contrast_summary_bar.png"),
  width = 24, height = 22, units = "cm"
)


### COMPARE DEGs ACROSS CONTRASTS + DATA SETS ###
# prep names for plotting
consensus_degs <- consensus_long %>%
  mutate(
    Dataset = str_to_title(dataset),
    Contrast = gsub(contrast, replacement = " day", pattern = "_"),
    Contrast = gsub(Contrast, replacement = " vs. ", pattern = "vs")
  )
str(consensus_degs)
unique(consensus_degs$Contrast)

compare_contrasts <- consensus_degs %>%
  select(gene_name, Contrast) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = Contrast, values_from = value, values_fill = 0
  ) %>%
  distinct() %>%
  column_to_rownames(var = "gene_name")

# set colors for upset plot
contrast_dataset_map <- consensus_degs %>%
  distinct(Contrast, Dataset)

dataset_colors <- c(
  "Cardo" = viridis(10)[[1]],
  "Shaw" = viridis(10)[[7]]
)

degs_upset <- upset(compare_contrasts,
  order.by = "freq",
  set.metadata = list(data = contrast_dataset_map, plots = list(
    list(
      type = "text",
      column = "Dataset",
      assign = 10,
      colors = dataset_colors
    ),
    list(
      type = "matrix_rows",
      column = "Dataset",
      assign = 10,
      colors = dataset_colors
    )
  ))
)
png(
  paste0(
    outdir, "upset_plot_compare_datasets.png"
  ),
  width = 30, height = 20, units = "cm", res = 300
)
degs_upset
dev.off()

### ORA RESULTS ACROSS SIGNATURES ###
# prep data for plotting
plot_df <- results_df %>%
  filter(significant == TRUE) %>%
  filter(term_size < 300) %>%
  mutate(
    neg_log10_p = -log10(p_value),
    # Simplify term names for better readability
    term_name = gsub("Homo sapiens ", "", term_name)
  )

# counts of terms
plot_df %>%
  group_by(query, source) %>%
  summarise(
    n_significant_terms = n(),
    .groups = "drop"
  ) %>%
  ggplot(aes(x = query, y = n_significant_terms, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    x = "Query Condition",
    y = "Count of Significant Terms",
    fill = "Data Source"
  ) +
  mytheme

ggsave(paste0(outdir, "ORA_term_counts_barplot.png"),
  dpi = 300, height = 22, width = 18, units = "cm"
)
# Create the comparative dot plot
plot_df %>%
  ggplot(aes(x = query, y = reorder(term_name, neg_log10_p))) +
  geom_point(aes(size = intersection_size, color = neg_log10_p)) +
  scale_color_viridis_c(
    name = "-log10(p-value)", option = "magma", end = 0.85
  ) +
  labs(
    x = "Query Condition",
    y = "Enriched Term",
    size = "Gene Overlap Size"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_grid(source ~ ., scales = "free_y", space = "free_y")

ggsave(paste0(outdir, "ORA_comparison_dotplot.png"),
  dpi = 300, height = 25, width = 22, units = "cm"
)

# GO:BP only
ora_plot <- plot_df %>%
  filter(!source %in% c("WP", "GO:MF", "GO:CC")) %>%
  ggplot(aes(x = query, y = reorder(term_name, neg_log10_p))) +
  geom_point(aes(size = intersection_size, color = neg_log10_p)) +
  scale_color_viridis_c(
    name = "-log10(p-value)", option = "magma", end = 0.85
  ) +
  labs(
    x = "Signature",
    y = "Enriched Term",
    size = "Gene Overlap Size"
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", face = "bold"),
    axis.text = element_text(size = 12)
  ) +
  facet_grid(source ~ ., scales = "free_y", space = "free_y")

ora_plot
ggsave(paste0(outdir, "ORA_comparison_GOBP_KEGG_dotplot.png"),
  dpi = 300, height = 25, width = 22, units = "cm"
)

### PERFORMANCE OF SIGNATURES BY VALIDATION GENES ###

performance_plot <- ggplot(performance_res, aes(x = Recall, y = Precision, color = F1_Score)) +
  geom_point(aes(size = -log10(P_Value))) +
  geom_label_repel(
    aes(label = Signature),
    seed = 42,
    size = 6.5,
    min.segment.length = 0.05,
    fontface = "bold",
    box.padding = 1.55, # Distance from point
    point.padding = 0.65, # Distance from point edge
    segment.color = "grey50", # Line connecting text to point
    max.overlaps = Inf # Ensure all labels show
  ) +
  scale_color_viridis_c(option = "magma", end = 0.85) +
  # expand limits to ensure labels aren't cut off
  expand_limits(
    x = c(0.01, max(performance_res$Recall) + 0.01),
    y = c(0.08, max(performance_res$Precision) + 0.01)
  ) +
  labs(
    x = "Recall (Sensitivity)",
    y = "Precision (Positive Predictive Value)"
  ) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.98, .05),
    legend.justification = c("right", "bottom")
  ) +
  mytheme

performance_plot
ggsave(paste0(outdir, "signature_performance_validation.png"),
  dpi = 300, height = 22, width = 18, units = "cm"
)

### RRA META-ANALYSIS SIGNATURE DEG VOLCANO PLOT ###
lfc_th <- 0.1 # logFC threshold
pval_th <- -log10(0.05) # p-value threshold

core_genes_label <- rra_meta %>%
  slice_max(RRA_pval, n = 10)

volcano_plot <- ggplot(rra_meta, aes(x = meta_logFC, y = RRA_pval)) +
  geom_point(aes(color = Direction), alpha = 0.75, size = 2.5) +
  scale_color_manual(values = c(
    "Upregulated" = viridis(20, option = "magma")[[17]],
    "Downregulated" = viridis(20, option = "magma")[[3]],
    "Non-Significant/Inconsistent" = "gray80"
  )) +
  # horizontal threshold line for RRA significance
  geom_hline(
    yintercept = pval_th, linetype = "dashed", color = "black"
  ) +
  # vertical threshold lines for logFC magnitude
  geom_vline(
    xintercept = c(-lfc_th, lfc_th),
    linetype = "dashed",
    color = "black"
  ) +
  labs(
    x = "MA Aggregate logFC (Magnitude)",
    y = "-log10(RobustRankAggregate P-value)",
    color = "Consensus Status"
  ) +
  theme(
    legend.position = "inside",
    legend.background = element_rect(fill = (alpha("white", 0.6))),
    legend.position.inside = c(1, .05),
    legend.justification = c("right", "bottom")
  ) +
  mytheme +
  # add text labels for top genes
  ggrepel::geom_text_repel(
    data = core_genes_label,
    aes(label = gene_name),
    seed = 42,
    size = 5.5,
    fontface = "bold",
    min.segment.length = 0,
    box.padding = 1.15,
    point.padding = 0.70,
    segment.color = "grey50",
    max.overlaps = Inf
  )
volcano_plot

ggsave(paste0(outdir, "RRAp_metalogFC_volcano.png"),
  dpi = 300, height = 30, width = 26, units = "cm"
)

#### FIGURE 3 ASSEMBLY ####
bottom_row <- plot_grid(performance_plot, volcano_plot, labels = c("B", "C"))
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
