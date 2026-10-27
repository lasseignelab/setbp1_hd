# script used to find consensus human signature from cardo and shaw data sets
set.seed(42)

suppressPackageStartupMessages({
  library(argparse)
  library(recount3)
  library(ggplot2)
  library(UpSetR)
  library(ComplexHeatmap)
  library(statmod)
  library(gprofiler2)
  library(RobustRankAggreg)
  library(metafor)
  library(circlize)
  library(tidyverse)
  library(readxl)
})
###### SET UP ######
# create parser object
parser <- ArgumentParser()

parser$add_argument("-i", "--indir",
  type = "character",
  help = "Path to the input data directory to access validation data"
)
parser$add_argument("-d", "--degs",
  type = "character",
  help = "Parent directory containing each dataset's DEG CSV files"
)
parser$add_argument("-n", "--numsig",
  type = "integer",
  help = "Number of top DEGs by method used to retain for signature"
)
parser$add_argument("-o", "--outdir",
  type = "character",
  help = "Path to the directory to write consensus signature results"
)
# get command line options, if help option encountered print help and exit
args <- parser$parse_args()
indir <- args$indir
degs <- args$degs
numsig <- args$numsig
outdir <- args$outdir

# create output dir if needed
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# helper function to read in DEGs
read_deg <- function(dataset_name, type = "significant_degs.csv") {
  file <- file.path(degs, dataset_name, type)
  if (!file.exists(file)) stop("Missing file: ", file)
  read_csv(file) %>% mutate(dataset = dataset_name)
}

# helper function for safely plotting (leading to warnings instead of halting script)
safely_plot <- function(expr, label) {
  tryCatch(expr, error = function(e) {
    message("PLOT FAILED (", label, "): ", conditionMessage(e))
    try(dev.off(), silent = TRUE)
    NULL
  })
}

# function to calculate performance metrics of signature by validation genes
calculate_sig_performance <- function(sig_genes, true_genes, universe_genes) {
  # Standardize inputs
  sig_genes <- unique(sig_genes)
  true_genes <- unique(true_genes)

  # Contingency Table Components
  tp <- length(intersect(sig_genes, true_genes)) # True Positives
  fp <- length(setdiff(sig_genes, true_genes)) # False Positives
  fn <- length(setdiff(true_genes, sig_genes)) # False Negatives
  tn <- length(universe_genes) - (tp + fp + fn) # True Negatives

  # Metrics
  precision <- ifelse((tp + fp) == 0, 0, tp / (tp + fp))
  recall <- ifelse((tp + fn) == 0, 0, tp / (tp + fn))
  f1 <- ifelse((precision + recall) == 0, 0,
    2 * (precision * recall) / (precision + recall)
  )

  # Statistical enrichment (Hypergeometric test)
  p_val <- phyper(tp - 1, length(true_genes),
    length(universe_genes) - length(true_genes),
    length(sig_genes),
    lower.tail = FALSE
  )

  data.frame(
    Precision = precision,
    Recall = recall,
    F1_Score = f1,
    P_Value = p_val,
    TP = tp,
    FP = fp
  )
}

# start timer
ptm <- proc.time()

###### READ IN DATA ######
cardo_degs <- read_deg("cardo")
shaw_degs <- read_deg("shaw")

# full universe of genes (all passing low-counts filtering)
cardo_unfilt <- read_deg("cardo", "unfiltered_degs.csv")
shaw_unfilt <- read_deg("shaw", "unfiltered_degs.csv")
combined_unfilt <- bind_rows(cardo_unfilt, shaw_unfilt)

# validation data
tanaka_degs <- bind_rows(read_deg("tanaka1"), read_deg("tanaka2"))
wong_fibro <- read_excel(paste0(
  indir, "wong/SourceData/SourceData6_forFIg5c_DEG.xlsx"
), sheet = 3, skip = 1) %>%
  select(gene = ensemblID, logFC = log2FoldChange, gene_name = symbol) %>%
  mutate(contrast = "fibro")


wong_neu <- read_excel(paste0(
  indir, "wong/SourceData/SourceData10_forFig7b_DEG.xlsx"
), sheet = 3, skip = 1) %>%
  select(gene = `...1`, logFC = log2FoldChange, gene_name = symbol) %>%
  mutate(contrast = "neu")

###### FIND CONSENSUS DEGS ######
unique(cardo_degs$contrast)
unique(shaw_degs$contrast)

# exclude redundant/biased contrasts
contrasts_skip <- c(
  "KOvsWT_21", "KOvsWT_NPC", "KOvsWT_all", "NPC_HD", "PATH_vs_PATH"
)

# combine human data sets and filter
combined_degs <- bind_rows(cardo_degs, shaw_degs)

consensus_degs <- combined_degs %>%
  filter(!contrast %in% contrasts_skip) %>%
  group_by(gene) %>%
  mutate(
    n_contrasts = n_distinct(contrast),
    n_datasets  = n_distinct(dataset)
  ) %>%
  ungroup() %>%
  filter(n_contrasts > 1)

# How much of the current signature rests on one dataset:
consensus_degs %>%
  distinct(gene, n_datasets) %>%
  count(n_datasets) %>%
  print()

# Strict alternative — run this as a sensitivity analysis, not as the default.
# Report how many signature genes survive it.
consensus_degs_strict <- consensus_degs %>%
  filter(n_datasets >= 2)

message(
  "Genes in >=2 contrasts: ", n_distinct(consensus_degs$gene),
  " | in >=2 datasets: ", n_distinct(consensus_degs_strict$gene)
)

# long format for plotting
consensus_long <- consensus_degs %>%
  select(gene, gene_name, contrast, dataset, logFC, adj.P.Val, se_logFC) %>%
  pivot_longer(
    cols = c(logFC, adj.P.Val, se_logFC),
    names_to = "metric",
    values_to = "value"
  )

# save consensus DEGs
write_csv(consensus_degs, paste0(outdir, "consensus_degs_unfiltered.csv"))
write_csv(consensus_long, paste0(outdir, "consensus_degs_unfiltered_long.csv"))

consensus_wide <- consensus_long %>%
  distinct() %>%
  unite("contrast_dataset_metric", contrast, dataset, metric) %>%
  pivot_wider(names_from = contrast_dataset_metric, values_from = value)
print(paste0("Consensus genes dataframe: "))
str(consensus_wide)
print(paste0("Number of consensus genes: "))
length(unique(consensus_long$gene))

###### PRIORTIZE DEGS TO MAKE SIGNATURE ######
### opt 1) high frequency DEGs ###
print(paste0("Finding high-frequency ('majority vote') DEGs..."))
freq_degs <- consensus_long %>%
  select(!c(metric, value)) %>%
  distinct() %>%
  add_count(gene)

### opt 2) meta-analysis per gene ###
print(paste0("Finding effect size meta analysis DEGs..."))
meta_results <- consensus_degs %>%
  group_by(gene) %>%
  nest() %>%
  mutate(meta = map(data, ~ tryCatch(
    rma(yi = .x$logFC, sei = .x$se_logFC, method = "REML"),
    error = function(e) {
      message("Error for gene: ", unique(.x$gene_name), " ", e$message)
      NULL
    }
  ))) %>%
  mutate(
    meta_logFC = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$b)     else NA_real_),
    meta_pval  = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$pval)  else NA_real_),
    meta_ci_lb = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$ci.lb) else NA_real_),
    meta_ci_ub = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$ci.ub) else NA_real_),
    # --- new ---
    meta_se    = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$se)    else NA_real_),
    meta_k     = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$k)     else NA_real_),
    meta_i2    = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$I2)    else NA_real_),
    meta_tau2  = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$tau2)  else NA_real_),
    meta_qep   = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$QEp)   else NA_real_)
  ) %>%
  select(-data, -meta) %>%
  ungroup()

# Use an explicit select() rather than positional indexing (`meta_results[, -2]`),
# which silently drops the wrong column once you add fields.
meta_stats <- meta_results %>%
  left_join(distinct(consensus_degs[, c("gene", "gene_name")]), by = "gene") %>%
  distinct()

# BH-adjust the meta p-values. You adjust within each contrast's DE test but
# never across the ~thousands of per-gene meta-analyses.
meta_stats <- meta_stats %>%
  mutate(meta_padj = p.adjust(meta_pval, method = "BH"))

summary(meta_stats)
write_csv(meta_stats, paste0(outdir, "metaanalysis_degs.csv"))

#### sign concordance ####
sign_stats <- consensus_degs %>%
  group_by(gene) %>%
  summarise(
    n_contrasts     = dplyr::first(n_contrasts),
    n_datasets      = dplyr::first(n_datasets),
    n_up            = sum(logFC > 0),
    n_down          = sum(logFC < 0),
    prop_concordant = max(sum(logFC > 0), sum(logFC < 0)) / dplyr::n(),
    sign_consistent = (sum(logFC > 0) == 0) | (sum(logFC < 0) == 0),
    max_abs_logFC   = max(abs(logFC)),
    .groups = "drop"
  )

# Sanity check on the two headline genes:
sign_stats %>%
  left_join(distinct(consensus_degs[, c("gene", "gene_name")]), by = "gene") %>%
  filter(gene_name %in% c("LIN28A", "SPON1", "ROBO1", "PSMG3", "ISL2")) %>%
  print()


### opt 3) robust rang aggregate -- create ranked lists per contrast ###
print(paste0("Finding robust rank aggregate (RRA) DEGs..."))
universe_n <- length(unique(combined_unfilt$gene))

ranked_lists <- lapply(unique(consensus_degs$contrast), function(c) {
  consensus_degs %>%
    filter(contrast == c) %>%
    arrange(adj.P.Val) %>%   # abs() was a no-op; adj.P.Val is already >= 0
    pull(gene)
})

rra <- aggregateRanks(ranked_lists, N = universe_n, exact = TRUE)
write_csv(rra, paste0(outdir, "rankaggregate_degs.csv"))

## Composite Scores from RRA, Meta Analysis, and Freq-Weight ##
print(
  paste0(
    "Finding composite scores by RRA, meta analysis, and frequency weight..."
  )
)

# convert DEGs to ENTREZ -- required for signaturesearch downstream
entrez_degs <- gconvert(
  query = consensus_degs$gene, target = "ENTREZGENE_ACC"
) %>%
  select(gene = input, entrez = target)

composite_scores <- freq_degs %>%
  select(gene, n) %>%
  left_join(entrez_degs, by = "gene") %>%
  distinct() %>%
  left_join(meta_stats, by = "gene") %>%
  left_join(rra, by = c("gene" = "Name")) %>%
  mutate(
    freq_weight = log2(n),
    rank_pval = Score,
    composite_score = -log10(meta_pval) + -log10(rank_pval) + freq_weight
  )

# save scores across methods for DEGs + composite scores
write_csv(
  composite_scores, paste0(outdir, "composite_aggregate_all_scores_degs.csv")
)

###### FILTER AND SAVE SIGNATURES ######

# save majority/high freq signature
top_freq <- composite_scores %>%
  select(gene, gene_name, entrez, n) %>%
  mutate(method = "majority") %>%
  slice_max(n = numsig, order_by = n)
print(
  paste0(
    length(
      unique(top_freq$gene)
    ), " DEGs with frequencies: ", unique(top_freq$n)
  )
)

write_csv(top_freq, paste0(outdir, "majorityfreq_consensus_signature.csv"))

# save meta analysis signature
meta_a <- composite_scores %>%
  select(gene, gene_name, entrez, meta_pval, meta_logFC) %>%
  distinct() %>%
  slice_min(meta_pval, n = numsig) %>%
  mutate(method = "meta_a")

write_csv(meta_a, paste0(outdir, "metaanalysis_consensus_signature.csv"))

# save RRA signature
rra_results <- composite_scores %>%
  distinct() %>%
  slice_min(order_by = Score, n = numsig) %>%
  mutate(method = "RRA") %>%
  select(gene, gene_name, entrez, Score, meta_logFC, method)

write_csv(rra_results, paste0(outdir, "rankaggregate_consensus_signature.csv"))

# save composite score signature
composite_signature <- composite_scores %>%
  mutate(method = "CompositeScore") %>%
  slice_max(composite_score, n = numsig)

write_csv(
  composite_signature, paste0(
    outdir, "compositescore_consensus_signature.csv"
  )
)

# RRA-META
lfc_th    <- 0.1
rra_th    <- 0.05    # threshold on the RRA Score (a p-value), not on -log10
meta_p_th <- 0.05
i2_th     <- 75      # I^2 (%) above this = substantial between-contrast heterogeneity

rra_meta <- rra %>%
  inner_join(meta_stats, by = c("Name" = "gene")) %>%
  left_join(sign_stats,  by = c("Name" = "gene")) %>%
  mutate(
    RRA_pval = -log10(Score),
    ci_width = meta_ci_ub - meta_ci_lb,
    ci_excludes_zero = sign(meta_ci_lb) == sign(meta_ci_ub),

    passes_rra  = Score < rra_th,
    passes_lfc  = abs(meta_logFC) > lfc_th,
    passes_meta = meta_pval < meta_p_th,
    passes_sign = sign_consistent,
    passes_het  = is.na(meta_i2) | meta_k < 3 | meta_i2 < i2_th,
    passes_data = n_datasets >= 2,

    # Your original call, kept so the two are directly comparable
    Direction = case_when(
      meta_logFC >  lfc_th & Score < rra_th ~ "Upregulated",
      meta_logFC < -lfc_th & Score < rra_th ~ "Downregulated",
      TRUE ~ "Non-Significant/Inconsistent"
    ),

    Tier = case_when(
      passes_rra & passes_lfc & passes_meta &
        passes_sign & passes_het & passes_data      ~ "Tier1_high_confidence",
      passes_rra & passes_lfc & passes_meta         ~ "Tier2_supported",
      Direction != "Non-Significant/Inconsistent"   ~ "Tier3_exploratory",
      TRUE                                          ~ "Excluded"
    )
  )

table(rra_meta$Tier)
table(rra_meta$Direction, rra_meta$Tier)

# Which criterion is doing the work:
rra_meta %>%
  filter(Direction != "Non-Significant/Inconsistent") %>%
  summarise(
    n               = dplyr::n(),
    fail_meta_p     = sum(!passes_meta,  na.rm = TRUE),
    fail_sign       = sum(!passes_sign,  na.rm = TRUE),
    fail_het        = sum(!passes_het,   na.rm = TRUE),
    fail_one_dataset = sum(!passes_data, na.rm = TRUE)
  ) %>%
  print()

rra_meta_filt <- rra_meta %>% filter(Direction != "Non-Significant/Inconsistent")
tier1 <- rra_meta %>% filter(Tier == "Tier1_high_confidence")

write_csv(rra_meta,      paste0(outdir, "RRA_meta_DEGs.csv"))
write_csv(rra_meta_filt, paste0(outdir, "RRA_meta_consensus_signature.csv"))
write_csv(tier1,         paste0(outdir, "RRA_meta_tier1_signature.csv"))


###### VALIDATION ######
## compare overlap between options for finding signature ##
compare_sigs <- top_freq %>%
  bind_rows(
    rra_results, meta_a, rra_results, composite_signature
  ) %>%
  select(gene, method) %>%
  mutate(method = ifelse(is.na(method), "Composite", method)) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = method, values_from = value, values_fill = 0) %>%
  column_to_rownames(var = "gene")

# plot DEGs across data sets for signature genes
dataset_sig_degs <- consensus_wide %>%
  select(gene, gene_name, contains("_logFC")) %>%
  inner_join(
    composite_scores %>%
      select(gene),
    by = "gene"
  ) %>%
  select(!contains("se_logFC")) %>%
  distinct()

safely_plot({
png(
  paste0(
    outdir, "upset_plot_compare_sigs.png"
  ),
  width = 30, height = 20, units = "cm", res = 300
)
upset(compare_sigs, sets = colnames(compare_sigs), order.by = "freq")
dev.off()
}, "upset_signatures")

# compare to Tanaka et al DEG
tanaka_genes <- gorth(
  query = tanaka_degs$gene,
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)

valid_degs <- tanaka_degs %>%
  inner_join(tanaka_genes, by = c("gene" = "input_ensg")) %>%
  select(
    gene = ortholog_ensg,
    logFC,
    padj = adj.P.Val, contrast, dataset
  ) %>%
  bind_rows(wong_fibro, wong_neu) %>%
  mutate(
    contrast = coalesce(contrast, "unknown"),
    dataset = coalesce(dataset, "wong")
  ) %>%
  filter(!is.na(gene))
write_csv(valid_degs, paste0(outdir, "validation_datasets_DEGs.csv"))

valid_degs_wide <- valid_degs %>%
  unite("contrast_dataset", contrast, dataset) %>%
  distinct(gene, contrast_dataset, .keep_all = TRUE) %>%
  select(gene, gene_name, contrast_dataset, logFC) %>%
  pivot_wider(
    names_from = contrast_dataset, values_from = logFC, values_fill = 0
  ) %>%
  # Use full_join to keep everything from both sides
  full_join(dataset_sig_degs, by = "gene", suffix = c(".valid", ".sig")) %>%
  mutate(
    # Merge the two name columns: if valid name is missing, use sig name
    gene_name = coalesce(gene_name.valid, gene_name.sig),
    # Clean up numeric NAs (for genes missing in one of the sets)
    across(where(is.numeric), ~ replace_na(.x, 0))
  ) %>%
  # Remove the temporary suffix columns
  select(-gene_name.valid, -gene_name.sig) %>%
  relocate(gene_name, .after = gene) %>%
  group_by(gene) %>%
  slice_max(
    order_by = rowSums(across(where(is.numeric))),
    n = 1,
    with_ties = FALSE
  ) %>%
  ungroup()

write_csv(
  valid_degs_wide, paste0(outdir, "validation_datasets_consensus_DEGs.csv")
)

# --- create column annotations --- #
mat_cols <- colnames(valid_degs_wide)[3:ncol(valid_degs_wide)]

# define groups
col_groups <- case_when(
  grepl("_cardo_logFC$", mat_cols) ~ "cardo",
  grepl("_shaw_logFC$", mat_cols) ~ "shaw",
  mat_cols == "meta_logFC" ~ "meta",
  TRUE ~ "validation"
)

# set colors by group
col_anno_colors <- c(
  "cardo"      = "#984EA3", # Purple
  "shaw"       = "#377EB8", # Blue
  "meta"       = "#4DAF4A", # Green
  "validation" = "black" # Black
)

# create column annotation object
column_ha <- HeatmapAnnotation(
  Source = col_groups,
  col = list(Source = col_anno_colors),
  show_annotation_name = FALSE
)

# color scale for logFC
max_lfc <- max(abs(meta_stats$meta_logFC))
col_fun <- circlize::colorRamp2(
  c(-max_lfc, 0, max_lfc), c("blue", "white", "red")
)

# prepare Matrix
plot_mat <- as.matrix(valid_degs_wide[, 3:ncol(valid_degs_wide)])
rownames(plot_mat) <- valid_degs_wide$gene

# plot
safely_plot({
pdf(paste0(outdir, "heatmap_all_degs_validation.pdf"), width = 10, height = 20)
Heatmap(
  plot_mat,
  name = "Log2FC",
  top_annotation = column_ha, # Add the new column annotation here
  row_labels = valid_degs_wide$gene_name,
  col = col_fun,
  cluster_columns = TRUE,
  show_row_dend = FALSE,
  row_names_gp = gpar(fontsize = 6),
  border = TRUE
)
dev.off()
}, "all_degs_validation")

# --- create row (gene) annotations --- #
### gene functional pathway annotations ###
pathway_map <- list(
  "Forebrain_Dev" = "GO:0030900", # forebrain dev
  "PI3K_AKT" = "REAC:R-HSA-5653656", # Reactome: PI3K/AKT Signaling
  "NF-kappa_B" = "KEGG:04064", # NF-kappa B signaling
  "PPAR" = "KEGG:03320", # PPAR signaling
  "MAPK" = "KEGG:04010", # MAPK signaling
  "Wnt" = "KEGG:04310" # KEGG: Wnt signaling pathway
)

target_pathway_ids <- unlist(pathway_map)

pathway_df <- gconvert(
  query = target_pathway_ids,
  organism = "hsapiens",
  target = "HGNC",
  mthreshold = Inf,
  filter_na = TRUE
)

pathways <- pathway_df %>%
  dplyr::select(input, target) %>%
  rename("STANDARD_NAME" = input, "GeneSymbol" = target) %>%
  filter(!is.na(GeneSymbol))
str(pathways)

# filter for DEGs that are also in at least one validation set
valid_degs_min <- valid_degs_wide %>%
  select(!any_of(c(
    "PATH3vsWT_shaw_logFC", "IPSC_HD_shaw_logFC",
    "HDvsWT_tanaka1", "HDvsWT_tanaka2"
  ))) %>%
  filter(rowSums(across(any_of(c("fibro_wong", "neu_wong")))) != 0) %>%
  filter(gene_name %in% pathways$GeneSymbol)
 
message("valid_degs_min rows: ", nrow(valid_degs_min))

# define row annotation based on established pathways
pathways_pivoted <- pathways %>%
  distinct() %>%
  filter(STANDARD_NAME %in% target_pathway_ids) %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = STANDARD_NAME, values_from = value, values_fill = 0
  )

# identify missing target columns and add them filled with 0
missing_cols <- setdiff(target_pathway_ids, names(pathways_pivoted))
if (length(missing_cols) > 0) {
  for (col in missing_cols) {
    pathways_pivoted[[col]] <- 0
  }
}

# join with DEGs and add pathway names
pathways_wide <- valid_degs_wide %>%
  rename("GeneSymbol" = gene_name) %>%
  left_join(pathways_pivoted, by = "GeneSymbol") %>%
  mutate(across(where(is.numeric), ~ replace_na(., 0))) %>%
  # assign a single, prioritized pathway using case_when
  mutate(
    Pathway = case_when(
      # Highest Priority
      `KEGG:04310` == 1 ~ "Wnt",
      `GO:0030900` == 1 ~ "Forebrain_Dev",
      `REAC:R-HSA-5653656` == 1 ~ "PI3K_AKT",
      `KEGG:04064` == 1 ~ "NF-kappa_B",
      `KEGG:04010` == 1 ~ "MAPK",
      `KEGG:03320` == 1 ~ "PPAR",
      TRUE ~ "Other"
    )
  ) %>%
  # prepare for ComplexHeatmap annotation
  column_to_rownames(var = "gene") %>%
  select(GeneSymbol, Pathway, all_of(target_pathway_ids))
str(pathways_wide)

# filter pathways for genes in plot matrix
plot_annotation <- pathways_wide %>%
  filter(row.names(.) %in% valid_degs_min$gene) %>%
  select(Pathway)

# Define a distinct color for each pathway
unique_pathways <- unique(plot_annotation$Pathway)
pathway_colors <- setNames(
  c(RColorBrewer::brewer.pal(
    min(length(unique_pathways), 8), "Set1"
  ), "#D3D3D3"),
  c(unique_pathways, "Other")
)

# create row annotation
gene_ha <- rowAnnotation(
  Pathway = plot_annotation$Pathway,
  col = list(Pathway = pathway_colors),
  show_annotation_name = TRUE
)

### plot heatmap ###
heatmap_obj <- Heatmap(
  as.matrix(valid_degs_min[, 3:length(valid_degs_min)]),
  name = "Log2(Fold Change)",
  col = col_fun,
  row_labels = valid_degs_min$gene_name,
  row_names_gp = gpar(fontsize = 6),
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  right_annotation = gene_ha,
  border = TRUE,
  column_title = "Significant SETBP1-HD DEGs Across Studies"
)

# save
safely_plot({
pdf(paste0(outdir, "validation_degs_min.pdf"), width = 10, height = 22)
draw(heatmap_obj, column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()
}, "validation_degs_min")

heatmap_by_method <- function(full_deg_df, genes, gene_data, title) {
  # filter and select sig degs
  plot_data <- full_deg_df %>%
    filter(gene %in% genes) %>%
    column_to_rownames("gene")
  print(head(plot_data))

  annotation_df <- gene_data %>%
    filter(row.names(.) %in% genes) %>%
    select(Pathway)

  gene_ha <- rowAnnotation(
    Pathway = annotation_df$Pathway,
    col = list(Pathway = pathway_colors)
  )

  # create plot
  hm <- Heatmap(
    as.matrix(plot_data[, 2:length(plot_data)]),
    name = "Log2(Fold Change)",
    col = col_fun,
    row_labels = plot_data$gene_name,
    row_names_gp = gpar(fontsize = 6),
    cluster_rows = TRUE,
    cluster_columns = TRUE,
    right_annotation = gene_ha,
    top_annotation = column_ha,
    border = TRUE,
    column_title = paste0(title, "-Based SETBP1-HD Consensus Signature")
  )
  list(plot = hm)
}

# plot heatmap for each method
sig_list <- list(
  Composite = c(composite_signature)$gene,
  RRA = c(rra_results)$gene,
  Meta = c(meta_a)$gene,
  Freq = c(top_freq)$gene
)
hm_list <- lapply(seq_along(sig_list), function(sig) {
  heatmap_by_method(
    full_deg_df = valid_degs_wide,
    genes = sig_list[[sig]],
    gene_data = pathways_wide,
    title = names(sig_list)[[sig]]
  )
})

# save all plots to one PDF
safely_plot({
pdf(
  paste0(outdir, "heatmap_validation_annotation_by_sig.pdf"),
  width = 10, height = 18
)
for (plot in hm_list) {
  print(plot)
}
dev.off()
}, "validation_anno_by_sig")

### RRA and Meta Analysis Volcano Plot ###
lfc_th <- 0.1 # logFC threshold
pval_th <- -log10(0.05) # p-value threshold

# compare RRA alone to RRA-meta
excluded <- setdiff(rra_results$gene_name, rra_meta_filt$gene_name)
excluded_genes <- rra_meta %>%
  filter(gene_name %in% excluded)

rra_meta_exc <- heatmap_by_method(
  valid_degs_wide,
  excluded_genes$Name,
  pathways_wide,
  title = "RRA Genes Excluded by Meta logFC Threshold"
)
# save all plots to one PDF
safely_plot({
pdf(
  paste0(outdir, "heatmap_rra_meta_excluded.pdf"),
  width = 10, height = 18
)
print(rra_meta_exc)
dev.off()
}, "rra_meta_excluded_genes")


###### ORA OF CONSENSUS SIGNATURES ######
# run ORA for each signature
ora_all <- gost(sig_list, sources = c("GO", "KEGG", "WP"), evcodes = TRUE)
results_df <- ora_all$result

write_csv(results_df, paste0(outdir, "ORA_signatures_results.csv"))

##### SIGNATURE PERFORMANCE VS VALIDATION SET #####
# background genes (all possible DEGs for signature)
background_g <- unique(c(
  cardo_degs$gene_name, shaw_degs$gene_name
))
print(paste0("Number of background genes: ", length(background_g)))

true_positives <- unique(c(
  pathways$GeneSymbol,
  tanaka_genes$ortholog_name,
  wong_neu$gene_name, wong_fibro$gene_name
))
print(paste0("Number of validation 'TP' genes: ", length(true_positives)))

sig_list <- list(
  Composite = c(composite_signature)$gene_name,
  RRA = c(rra_results)$gene_name,
  Meta = c(meta_a)$gene_name,
  Freq = c(top_freq)$gene_name
)

# calculate scores for each signature
performance_results <- map_df(names(sig_list), function(name) {
  calculate_sig_performance(sig_list[[name]], true_positives, background_g) %>%
    mutate(Signature = name)
}) %>%
  select(Signature, everything()) %>%
  arrange(desc(F1_Score))

print(paste0("Performance metrics for signatures vs validation genes: "))
print(performance_results)

write_csv(
  performance_results, paste0(outdir, "signature_performance_metrics.csv")
)

### CLAUDE RECOMMENDED DIAGNOSTICS AND VALIDATION ###
rank_disc <- rra_meta %>%
  filter(!is.na(meta_pval)) %>%
  mutate(
    rra_rank  = rank(Score,     ties.method = "average"),
    meta_rank = rank(meta_pval, ties.method = "average")
  )

rho <- cor.test(rank_disc$rra_rank, rank_disc$meta_rank, method = "spearman")
print(rho)

# Does CI width grow as RRA rank improves? (It should not.)
print(cor.test(rank_disc$rra_rank, rank_disc$ci_width, method = "spearman"))

p_disc <- ggplot(rank_disc, aes(x = rra_rank, y = -log10(meta_pval))) +
  geom_point(aes(colour = Tier), alpha = 0.6, size = 1.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  scale_x_continuous(trans = "log10") +
  labs(
    x = "RRA rank (log scale, 1 = best)",
    y = expression(-log[10]~"(meta-analysis p)"),
    title = "RRA rank vs meta-analysis evidence",
    subtitle = paste0("Spearman rho = ", signif(rho$estimate, 3))
  ) +
  theme_bw()

safely_plot({
ggsave(paste0(outdir, "rra_vs_meta_discordance.png"), p_disc,
       width = 18, height = 13, units = "cm", dpi = 300)
}, "rra_meta_discordance")

# Forest plots for the genes the claims rest on. This is the figure that
# shows a reader immediately why a gene did or did not survive.
forest_genes <- c("LIN28A", "SPON1", "ROBO1", "PSMG3", "ISL2", "DRGX")

safely_plot({
pdf(paste0(outdir, "forest_headline_genes.pdf"), width = 7, height = 5)
for (g in forest_genes) {
  d <- consensus_degs %>% filter(gene_name == g)
  if (nrow(d) < 2) next
  fit <- tryCatch(rma(yi = d$logFC, sei = d$se_logFC, method = "REML"),
                  error = function(e) NULL)
  if (is.null(fit)) next
  forest(fit, slab = paste(d$dataset, d$contrast, sep = ": "),
         xlab = "log2 fold change", main = g)
}
dev.off()
}, "forest_headline_genes")


# ==============================================================================
# G. LEAVE-ONE-DATASET-OUT SENSITIVITY  --  new block
# ==============================================================================
# Directly tests dataset dominance: rerun RRA with each dataset removed and see
# which signature genes survive. A gene that disappears when one dataset is
# dropped was never a consensus gene.

loo_rra <- function(drop_dataset) {
  sub <- consensus_degs %>% filter(dataset != drop_dataset)
  if (nrow(sub) == 0) return(NULL)
  lists <- lapply(unique(sub$contrast), function(c) {
    sub %>% filter(contrast == c) %>% arrange(adj.P.Val) %>% pull(gene)
  })
  lists <- lists[lengths(lists) > 0]
  if (length(lists) < 2) return(NULL)
  aggregateRanks(lists, N = universe_n, exact = TRUE) %>%
    rename(!!paste0("Score_drop_", drop_dataset) := Score)
}

loo_results <- reduce(
  lapply(unique(consensus_degs$dataset), loo_rra),
  full_join, by = "Name"
)

loo_compare <- rra_meta_filt %>%
  select(Name, gene_name, Score, meta_logFC, meta_pval, Tier, n_datasets) %>%
  left_join(loo_results, by = "Name") %>%
  mutate(across(starts_with("Score_drop_"), ~ .x < rra_th,
                .names = "survives_{.col}"))

write_csv(loo_compare, paste0(outdir, "leave_one_dataset_out_comparison.csv"))

loo_compare %>%
  summarise(across(starts_with("survives_"), ~ sum(.x, na.rm = TRUE))) %>%
  print()


# ==============================================================================
# H. PERFORMANCE BY TIER  --  extends your existing validation block
# ==============================================================================
# Reuses calculate_sig_performance() to show empirically whether the extra
# criteria improve precision against the validation genes. This is the
# justification for the filter, not an assertion that it is better.

tiered_sig_list <- list(
  Current_RRA_meta = rra_meta_filt$gene_name,
  Tier1            = tier1$gene_name,
  Tier1_plus_2     = rra_meta %>%
    filter(Tier %in% c("Tier1_high_confidence", "Tier2_supported")) %>%
    pull(gene_name),
  RRA_only         = rra_results$gene_name,
  Meta_only        = meta_a$gene_name
)

tier_performance <- map_df(names(tiered_sig_list), function(nm) {
  calculate_sig_performance(
    tiered_sig_list[[nm]], true_positives, background_g
  ) %>%
    mutate(Signature = nm, N_genes = length(unique(tiered_sig_list[[nm]])))
}) %>%
  select(Signature, N_genes, everything()) %>%
  arrange(desc(F1_Score))

print(tier_performance)
write_csv(tier_performance, paste0(outdir, "tier_performance_metrics.csv"))


# ==============================================================================
# I. SEPARATE ABSENCE  --  fix for valid_degs_wide
# ==============================================================================
# `values_fill = 0` and `replace_na(.x, 0)` make "not significant in this
# contrast" indistinguishable from "logFC is zero". In the heatmaps both render
# white, which reads as "measured, no change". Keep NA and build a companion
# presence matrix so the two can be told apart.
#
# In pivot_wider(): drop `values_fill = 0`
# In the across(): replace `~ replace_na(.x, 0)` with nothing (keep NA)
# Then, for Heatmap(), pass na_col = "grey90" so absence is visibly distinct.

presence_mat <- consensus_degs %>%
  mutate(present = 1L) %>%
  select(gene, contrast, dataset, present) %>%
  unite("contrast_dataset", contrast, dataset) %>%
  distinct() %>%
  pivot_wider(names_from = contrast_dataset, values_from = present,
              values_fill = 0L)

write_csv(presence_mat, paste0(outdir, "deg_presence_matrix.csv"))

#### END ####

# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
