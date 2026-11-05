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

# helper function for safely plotting (warns instead of halting the script).
# NOTE: plot objects are not auto-printed inside `expr`, so every plot call
# within it must be wrapped in print() or draw() or the saved file is blank.
safely_plot <- function(expr, label) {
  tryCatch(expr, error = function(e) {
    message("PLOT FAILED (", label, "): ", conditionMessage(e))
    try(dev.off(), silent = TRUE)
    NULL
  })
}

# helper function to cluster a logFC matrix that uses NA for "not a DEG in this
# contrast". Distances are computed with NA treated as 0 (as before), while the
# heatmap itself shows NA in grey so absence is distinguishable from logFC = 0.
cluster_na_matrix <- function(m) {
  m0 <- m
  m0[is.na(m0)] <- 0
  list(
    rows = if (nrow(m0) > 1) hclust(dist(m0)) else FALSE,
    cols = if (ncol(m0) > 1) hclust(dist(t(m0))) else FALSE
  )
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
# KOvsWT_21 (Cardo day 21): logFC is anti-correlated with Cardo day 15 and
#   disagrees in direction with Shaw for most shared DEGs (see sig_diagnostic.R)
contrasts_skip <- c(
  "KOvsWT_21", "KOvsWT_NPC", "KOvsWT_all", "NPC_HD", "PATH_vs_PATH"
)

# combine human data sets and filter
combined_degs <- bind_rows(cardo_degs, shaw_degs)

# keep genes that are DEGs in more than one contrast, and record how many
# contrasts and datasets support each gene
consensus_degs <- combined_degs %>%
  filter(!contrast %in% contrasts_skip) %>%
  group_by(gene) %>%
  mutate(
    n_contrasts = n_distinct(contrast),
    n_datasets  = n_distinct(dataset)
  ) %>%
  ungroup() %>%
  filter(n_contrasts > 1)

# number of consensus genes supported by 1 vs 2 datasets
consensus_degs %>%
  distinct(gene, n_datasets) %>%
  count(n_datasets) %>%
  print()

# stricter sensitivity set: genes that are DEGs in both datasets
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
# random-effects meta-analysis; REML can fail to converge with k = 2-4
# contrasts (e.g. GRM3), so fall back to DerSimonian-Laird for those genes
fit_meta <- function(d) {
  fit <- tryCatch(
    rma(yi = d$logFC, sei = d$se_logFC, method = "REML"),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    fit <- tryCatch(
      rma(yi = d$logFC, sei = d$se_logFC, method = "DL"),
      error = function(e) {
        message("Meta-analysis failed for gene: ", unique(d$gene_name), " ",
                e$message)
        NULL
      }
    )
  }
  fit
}

meta_results <- consensus_degs %>%
  group_by(gene) %>%
  nest() %>%
  mutate(meta = map(data, fit_meta)) %>%
  mutate(
    meta_method = map_chr(meta, ~ if (!is.null(.x)) .x$method else NA_character_),
    meta_logFC = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$b)     else NA_real_),
    meta_pval  = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$pval)  else NA_real_),
    meta_ci_lb = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$ci.lb) else NA_real_),
    meta_ci_ub = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$ci.ub) else NA_real_),
    # standard error, number of contrasts, and heterogeneity statistics
    meta_se   = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$se)    else NA_real_),
    meta_k     = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$k)     else NA_real_),
    meta_i2    = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$I2)    else NA_real_),
    meta_tau2  = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$tau2)  else NA_real_),
    meta_qep   = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$QEp)   else NA_real_)
  ) %>%
  select(-data, -meta) %>%
  ungroup()

meta_stats <- meta_results %>%
  left_join(distinct(consensus_degs[, c("gene", "gene_name")]), by = "gene") %>%
  distinct()

# BH-adjust meta-analysis p-values across all genes tested
meta_stats <- meta_stats %>%
  mutate(meta_padj = p.adjust(meta_pval, method = "BH"))

summary(meta_stats)
write_csv(meta_stats, paste0(outdir, "metaanalysis_degs.csv"))

#### sign concordance of logFC across contrasts ####
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

# print sign concordance for key genes
sign_stats %>%
  left_join(distinct(consensus_degs[, c("gene", "gene_name")]), by = "gene") %>%
  filter(gene_name %in% c("LIN28A", "SPON1", "ROBO1", "PSMG3", "ISL2", "CACNB2")) %>%
  print()


### opt 3) robust rank aggregate -- create ranked lists per contrast ###
print(paste0("Finding robust rank aggregate (RRA) DEGs..."))
# Each list is the COMPLETE set of genes tested in that contrast (all genes
# passing low-count filtering), not only the consensus DEGs. Ranking only the
# filtered DEGs against N = universe makes the Score mostly reflect how many
# contrasts a gene is a DEG in (presence/absence), not its rank position.
# P.Value gives the same order as adj.P.Val but without BH ties; |t| breaks
# the remaining ties among underflowed (P = 0) genes. Ranking is
# direction-agnostic; sign agreement is handled by the meta-analysis/tiers.
tested_unfilt <- combined_unfilt %>%
  filter(!contrast %in% contrasts_skip)

universe_n <- n_distinct(tested_unfilt$gene)

rank_contrast_lists <- function(df) {
  lists <- df %>%
    distinct(contrast, gene, .keep_all = TRUE) %>%
    arrange(P.Value, desc(abs(t))) %>%
    split(.$contrast) %>%
    map(~ .x$gene)
  lists[lengths(lists) > 0]
}

ranked_lists <- rank_contrast_lists(tested_unfilt)
message(
  "RRA: ", length(ranked_lists), " contrasts, ", universe_n, " genes ranked"
)

# aggregate over all tested genes, then keep the consensus genes
rra_all <- aggregateRanks(ranked_lists, N = universe_n, exact = TRUE)
rra <- rra_all %>%
  filter(Name %in% consensus_degs$gene)
write_csv(rra, paste0(outdir, "rankaggregate_degs.csv"))

## Composite Scores from RRA, Meta Analysis, and Freq-Weight ##
print(
  paste0(
    "Finding composite scores by RRA, meta analysis, and frequency weight..."
  )
)

# convert DEGs to ENTREZ -- required for signaturesearch downstream
# gconvert can return several Entrez IDs per Ensembl gene; keep the first so
# each gene occupies one row (and one slot in each top-N signature)
entrez_raw <- gconvert(
  query = unique(consensus_degs$gene), target = "ENTREZGENE_ACC"
)
message(
  "Ensembl genes with >1 Entrez ID: ",
  sum(table(entrez_raw$input) > 1)
)
entrez_degs <- entrez_raw %>%
  select(gene = input, entrez = target) %>%
  distinct(gene, .keep_all = TRUE)

# composite: average percentile rank of the three components. Summing them
# raw let -log10(meta_pval) (range ~0-240) swamp -log10(RRA Score) (~0-6)
# and log2(n) (~1-2). Higher composite_score = stronger evidence.
composite_scores <- freq_degs %>%
  distinct(gene, n) %>%
  left_join(entrez_degs, by = "gene") %>%
  left_join(meta_stats, by = "gene") %>%
  left_join(rra, by = c("gene" = "Name")) %>%
  mutate(
    pct_meta = percent_rank(-meta_pval),
    pct_rra  = percent_rank(-Score),
    pct_freq = percent_rank(n),
    composite_score = (pct_meta + pct_rra + pct_freq) / 3
  )
stopifnot(!anyDuplicated(composite_scores$gene))

# save scores across methods for DEGs + composite scores
write_csv(
  composite_scores, paste0(outdir, "composite_aggregate_all_scores_degs.csv")
)

###### FILTER AND SAVE SIGNATURES ######

# save majority/high freq signature
top_freq <- composite_scores %>%
  select(gene, gene_name, entrez, n) %>%
  left_join(select(meta_stats, gene, meta_pval), by = "gene") %>%
  mutate(method = "majority") %>%
  # n has heavy ties; break them on meta-analysis p so the signature is
  # exactly numsig genes (size-matched with the other methods)
  arrange(desc(n), meta_pval) %>%
  slice_head(n = numsig) %>%
  select(-meta_pval)
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
  slice_min(meta_pval, n = numsig, with_ties = FALSE) %>%
  mutate(method = "meta_a")

write_csv(meta_a, paste0(outdir, "metaanalysis_consensus_signature.csv"))

# save RRA signature
rra_results <- composite_scores %>%
  distinct() %>%
  slice_min(order_by = Score, n = numsig, with_ties = FALSE) %>%
  mutate(method = "RRA") %>%
  select(gene, gene_name, entrez, Score, meta_logFC, method)

write_csv(rra_results, paste0(outdir, "rankaggregate_consensus_signature.csv"))

# save composite score signature
composite_signature <- composite_scores %>%
  mutate(method = "CompositeScore") %>%
  slice_max(composite_score, n = numsig, with_ties = FALSE)

write_csv(
  composite_signature, paste0(
    outdir, "compositescore_consensus_signature.csv"
  )
)

###### RRA-META SIGNATURE + CONFIDENCE TIERS ######
# signature_call: RRA Score + meta-analysis |logFC| (the signature membership
#   call, with the sign of meta_logFC attached). It does NOT test the
#   meta-analysis p-value, so "Not in signature" does not mean a gene is
#   non-significant or inconsistent.
# meta_direction: direction from the meta-analysis alone (meta_logFC sign,
#   gated on meta_padj and sign consistency across contrasts).
# Tier: adds meta-analysis p, sign consistency, heterogeneity, and dataset
#   support to grade confidence within the signature. Only Tier1 requires
#   sign consistency, low heterogeneity and support from both datasets;
#   Tier2 genes may be single-dataset and/or heterogeneous (see the
#   n_datasets, passes_sign and het_exempt_k2 columns).
lfc_th    <- 0.1     # meta-analysis |logFC| threshold
rra_th    <- 0.05    # threshold on the RRA Score (a p-value)
meta_p_th <- 0.05    # threshold on the (unadjusted) meta-analysis p-value
i2_th     <- 75      # I^2 (%) above this = substantial between-contrast heterogeneity

rra_meta <- rra %>%
  inner_join(meta_stats, by = c("Name" = "gene")) %>%
  left_join(sign_stats,  by = c("Name" = "gene")) %>%
  mutate(
    # -log10 of the RRA Score for plotting; threshold Score itself
    neg_log10_RRA = -log10(Score),
    ci_width = meta_ci_ub - meta_ci_lb,
    ci_excludes_zero = sign(meta_ci_lb) == sign(meta_ci_ub),

    passes_rra  = Score < rra_th,
    passes_lfc  = abs(meta_logFC) > lfc_th,
    passes_meta = meta_pval < meta_p_th,
    passes_sign = sign_consistent,
    # I^2 is unstable with k = 2 contrasts, so those genes are exempt from the
    # heterogeneity filter; het_exempt_k2 flags the ones that would fail it
    het_exempt_k2 = meta_k < 3 & !is.na(meta_i2) & meta_i2 >= i2_th,
    passes_het  = is.na(meta_i2) | meta_k < 3 | meta_i2 < i2_th,
    passes_data = n_datasets >= 2,

    signature_call = case_when(
      meta_logFC >  lfc_th & Score < rra_th ~ "Upregulated",
      meta_logFC < -lfc_th & Score < rra_th ~ "Downregulated",
      TRUE ~ "Not in signature"
    ),

    meta_direction = case_when(
      is.na(meta_padj) | meta_padj >= meta_p_th ~ "Not significant",
      !sign_consistent                          ~ "Inconsistent sign",
      meta_logFC > 0                            ~ "Up",
      TRUE                                      ~ "Down"
    ),

    Tier = case_when(
      passes_rra & passes_lfc & passes_meta &
        passes_sign & passes_het & passes_data      ~ "Tier1_high_confidence",
      passes_rra & passes_lfc & passes_meta         ~ "Tier2_supported",
      signature_call != "Not in signature"          ~ "Tier3_exploratory",
      TRUE                                          ~ "Excluded"
    )
  )

table(rra_meta$Tier)
table(rra_meta$signature_call, rra_meta$Tier)
table(rra_meta$meta_direction, rra_meta$signature_call)

# number of signature genes failing each additional tier criterion
rra_meta %>%
  filter(signature_call != "Not in signature") %>%
  summarise(
    n                = dplyr::n(),
    fail_meta_p      = sum(!passes_meta,  na.rm = TRUE),
    fail_sign        = sum(!passes_sign,  na.rm = TRUE),
    fail_het         = sum(!passes_het,   na.rm = TRUE),
    het_exempt_k2    = sum(het_exempt_k2, na.rm = TRUE),
    fail_one_dataset = sum(!passes_data,  na.rm = TRUE)
  ) %>%
  print()

# same, by tier
rra_meta %>%
  group_by(Tier) %>%
  summarise(
    n                = dplyr::n(),
    fail_sign        = sum(!passes_sign,  na.rm = TRUE),
    fail_het         = sum(!passes_het,   na.rm = TRUE),
    het_exempt_k2    = sum(het_exempt_k2, na.rm = TRUE),
    single_dataset   = sum(!passes_data,  na.rm = TRUE)
  ) %>%
  print()

rra_meta_filt <- rra_meta %>% filter(signature_call != "Not in signature")
tier1 <- rra_meta %>% filter(Tier == "Tier1_high_confidence")

write_csv(rra_meta,      paste0(outdir, "RRA_meta_DEGs.csv"))
write_csv(rra_meta_filt, paste0(outdir, "RRA_meta_consensus_signature.csv"))
write_csv(tier1,         paste0(outdir, "RRA_meta_tier1_signature.csv"))


###### VALIDATION ######
## compare overlap between options for finding signature ##
compare_sigs <- top_freq %>%
  bind_rows(
    rra_results, meta_a, composite_signature
  ) %>%
  select(gene, method) %>%
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
print(upset(compare_sigs, sets = colnames(compare_sigs), order.by = "freq"))
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

# NA = not a DEG in that contrast/dataset (kept distinct from logFC = 0)
# Tanaka rows carry no gene_name (Wong rows do), so pivot on gene alone and
# attach one name per gene afterwards; pivoting on gene + gene_name split
# genes found in both validation studies across two rows.
valid_names <- bind_rows(
  distinct(dataset_sig_degs, gene, gene_name),
  valid_degs %>% filter(!is.na(gene_name)) %>% distinct(gene, gene_name),
  tanaka_genes %>% select(gene = ortholog_ensg, gene_name = ortholog_name)
) %>%
  # priority: consensus name, then Wong symbol, then Tanaka ortholog name
  filter(!is.na(gene), !is.na(gene_name)) %>%
  distinct(gene, .keep_all = TRUE)

valid_degs_wide <- valid_degs %>%
  unite("contrast_dataset", contrast, dataset) %>%
  distinct(gene, contrast_dataset, .keep_all = TRUE) %>%
  select(gene, contrast_dataset, logFC) %>%
  pivot_wider(names_from = contrast_dataset, values_from = logFC) %>%
  # keep genes from either the validation sets or the consensus set
  full_join(select(dataset_sig_degs, -gene_name), by = "gene") %>%
  left_join(valid_names, by = "gene") %>%
  relocate(gene_name, .after = gene)
stopifnot(!anyDuplicated(valid_degs_wide$gene))

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
# na.rm: genes whose meta-analysis did not converge (e.g. GRM3) have NA logFC
max_lfc <- max(abs(meta_stats$meta_logFC), na.rm = TRUE)
col_fun <- circlize::colorRamp2(
  c(-max_lfc, 0, max_lfc), c("blue", "white", "red")
)

# prepare Matrix
plot_mat <- as.matrix(valid_degs_wide[, 3:ncol(valid_degs_wide)])
rownames(plot_mat) <- valid_degs_wide$gene
plot_clust <- cluster_na_matrix(plot_mat)

# plot
safely_plot({
pdf(paste0(outdir, "heatmap_all_degs_validation.pdf"), width = 10, height = 20)
draw(Heatmap(
  plot_mat,
  name = "Log2FC",
  top_annotation = column_ha,
  row_labels = valid_degs_wide$gene_name,
  col = col_fun,
  na_col = "grey90",
  cluster_rows = plot_clust$rows,
  cluster_columns = plot_clust$cols,
  show_row_dend = FALSE,
  row_names_gp = gpar(fontsize = 6),
  border = TRUE
))
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
  filter(if_any(any_of(c("fibro_wong", "neu_wong")), ~ !is.na(.x))) %>%
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
min_mat <- as.matrix(valid_degs_min[, 3:length(valid_degs_min)])
min_clust <- cluster_na_matrix(min_mat)
heatmap_obj <- Heatmap(
  min_mat,
  name = "Log2(Fold Change)",
  col = col_fun,
  na_col = "grey90",
  row_labels = valid_degs_min$gene_name,
  row_names_gp = gpar(fontsize = 6),
  cluster_rows = min_clust$rows,
  cluster_columns = min_clust$cols,
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
  mat <- as.matrix(plot_data[, 2:length(plot_data)])
  clust <- cluster_na_matrix(mat)
  hm <- Heatmap(
    mat,
    name = "Log2(Fold Change)",
    col = col_fun,
    na_col = "grey90",
    row_labels = plot_data$gene_name,
    row_names_gp = gpar(fontsize = 6),
    cluster_rows = clust$rows,
    cluster_columns = clust$cols,
    right_annotation = gene_ha,
    top_annotation = column_ha,
    border = TRUE,
    column_title = paste0(title, "-Based SETBP1-HD Consensus Signature")
  )
  hm
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
  draw(plot)
}
dev.off()
}, "validation_anno_by_sig")

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
draw(rra_meta_exc)
dev.off()
}, "rra_meta_excluded_genes")


### heatmaps annotated by RRA-meta confidence tier ###
tier_colors <- c(
  "Tier1_high_confidence" = "#1B7837",
  "Tier2_supported"       = "#7FBF7B",
  "Tier3_exploratory"     = "#C7E9C0",
  "Excluded"              = "#969696",
  "Not in RRA-meta"       = "white"
)

# builds a logFC heatmap (per-contrast + validation + meta-analysis logFC) for
# the given Ensembl IDs, with a Tier row annotation
heatmap_by_tier <- function(genes, title, show_names = TRUE) {
  plot_data <- valid_degs_wide %>%
    filter(gene %in% genes) %>%
    left_join(
      rra_meta %>% select(Name, meta_logFC, Tier),
      by = c("gene" = "Name")
    ) %>%
    mutate(
      Tier = factor(
        replace_na(Tier, "Not in RRA-meta"), levels = names(tier_colors)
      )
    )

  mat <- plot_data %>%
    select(-gene, -gene_name, -Tier) %>%
    as.matrix()
  rownames(mat) <- plot_data$gene
  clust <- cluster_na_matrix(mat)

  mat_groups <- case_when(
    grepl("_cardo_logFC$", colnames(mat)) ~ "cardo",
    grepl("_shaw_logFC$", colnames(mat)) ~ "shaw",
    colnames(mat) == "meta_logFC" ~ "meta",
    TRUE ~ "validation"
  )

  Heatmap(
    mat,
    name = "Log2(Fold Change)",
    col = col_fun,
    na_col = "grey90",
    row_labels = plot_data$gene_name,
    show_row_names = show_names,
    row_names_gp = gpar(fontsize = 7),
    cluster_rows = clust$rows,
    cluster_columns = clust$cols,
    show_row_dend = FALSE,
    left_annotation = rowAnnotation(
      Tier = plot_data$Tier,
      col = list(Tier = tier_colors[levels(droplevels(plot_data$Tier))])
    ),
    top_annotation = HeatmapAnnotation(
      Source = mat_groups,
      col = list(Source = col_anno_colors),
      show_annotation_name = FALSE
    ),
    border = TRUE,
    column_title = title
  )
}

# 1) all signature genes (Tier 1-3; i.e. not Excluded)
tiered_genes <- rra_meta %>% filter(Tier != "Excluded") %>% pull(Name)

safely_plot({
pdf(paste0(outdir, "heatmap_tiered_signature.pdf"), width = 10, height = 20)
draw(heatmap_by_tier(
  tiered_genes,
  title = paste0("RRA-Meta Signature Genes by Tier (n = ", length(tiered_genes), ")"),
  show_names = length(tiered_genes) <= 150
))
dev.off()
}, "tiered_signature")

# 2) genes of interest (all tiers, including Excluded)
goi <- c(
  "TRIM71", "HMGA2", "IGF2BP1", "IGF2BP2", "IGF2BP3", "MYCN",
  "CCND1", "ARID3B", "DNMT3A", "DNMT3B"
)
goi_genes <- valid_degs_wide %>%
  filter(gene_name %in% goi | str_starts(gene_name, "CACN")) %>%
  pull(gene)

# genes absent here were not a DEG in any contrast or validation set
message(
  "Genes of interest not found among DEGs: ",
  paste(setdiff(goi, valid_degs_wide$gene_name), collapse = ", ")
)

safely_plot({
pdf(paste0(outdir, "heatmap_genes_of_interest_tiers.pdf"), width = 9, height = 8)
draw(heatmap_by_tier(
  goi_genes,
  title = "Genes of Interest (incl. CACN*) by Tier"
))
dev.off()
}, "goi_tiers")


###### ORA OF CONSENSUS SIGNATURES ######
# run ORA for each signature
ora_all <- gost(sig_list, sources = c("GO", "KEGG", "WP"), evcodes = TRUE)
results_df <- ora_all$result

write_csv(results_df, paste0(outdir, "ORA_signatures_results.csv"))

##### SIGNATURE PERFORMANCE VS VALIDATION SET #####
# background = every gene tested in the included contrasts (the same universe
# the RRA ranks), not only significant DEGs, and not the skipped contrasts
background_g <- unique(na.omit(tested_unfilt$gene_name))
print(paste0("Number of background genes: ", length(background_g)))

# two reference sets, both restricted to the background:
#  - deg_only: DEGs from the independent Tanaka (mouse) and Wong studies
#  - pathway_plus_deg: deg_only plus every gene in the six broad pathways
#    used for heatmap annotation (much larger; pathway membership is not
#    independent validation, so method rankings should hold under deg_only)
validation_sets <- list(
  deg_only = intersect(unique(c(
    tanaka_genes$ortholog_name,
    wong_neu$gene_name, wong_fibro$gene_name
  )), background_g),
  pathway_plus_deg = intersect(unique(c(
    pathways$GeneSymbol,
    tanaka_genes$ortholog_name,
    wong_neu$gene_name, wong_fibro$gene_name
  )), background_g)
)
print(lengths(validation_sets))

# kept for the downstream code/figures that use a single reference set
true_positives <- validation_sets$pathway_plus_deg

sig_list <- list(
  Composite = c(composite_signature)$gene_name,
  RRA = c(rra_results)$gene_name,
  Meta = c(meta_a)$gene_name,
  Freq = c(top_freq)$gene_name
)

# calculate scores for each signature against each reference set
performance_by_set <- map_df(names(validation_sets), function(vs) {
  map_df(names(sig_list), function(name) {
    calculate_sig_performance(
      sig_list[[name]], validation_sets[[vs]], background_g
    ) %>%
      mutate(Signature = name, Validation_set = vs)
  })
}) %>%
  select(Validation_set, Signature, everything()) %>%
  arrange(Validation_set, desc(F1_Score))

print(paste0("Performance metrics for signatures vs validation genes: "))
print(performance_by_set)

write_csv(
  performance_by_set,
  paste0(outdir, "signature_performance_by_validation_set.csv")
)
# pathway_plus_deg set only, in the format figure_3.R reads
performance_results <- performance_by_set %>%
  filter(Validation_set == "pathway_plus_deg") %>%
  select(-Validation_set)
write_csv(
  performance_results, paste0(outdir, "signature_performance_metrics.csv")
)

###### SIGNATURE DIAGNOSTICS ######
### agreement between RRA rank and meta-analysis evidence ###
rank_disc <- rra_meta %>%
  filter(!is.na(meta_pval)) %>%
  mutate(
    rra_rank  = rank(Score,     ties.method = "average"),
    meta_rank = rank(meta_pval, ties.method = "average")
  )

rho <- cor.test(rank_disc$rra_rank, rank_disc$meta_rank, method = "spearman")
print(rho)

# check whether meta-analysis CI width tracks RRA rank
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

### forest plots of per-contrast logFC for key genes ###
forest_genes <- c("LIN28A", "SPON1", "ROBO1", "PSMG3", "ISL2", "DRGX", "CACNB2")

safely_plot({
pdf(paste0(outdir, "forest_headline_genes.pdf"), width = 7, height = 5)
for (g in forest_genes) {
  d <- consensus_degs %>% filter(gene_name == g)
  if (nrow(d) < 2) next
  fit <- fit_meta(d)
  if (is.null(fit)) next
  forest(fit, slab = paste(d$dataset, d$contrast, sep = ": "),
         xlab = "log2 fold change", main = g)
}
dev.off()
}, "forest_headline_genes")


### leave-one-dataset-out sensitivity ###
# rerun RRA with each dataset removed to test whether a signature gene depends
# on a single dataset

# (full per-contrast ranked lists, as in the main RRA)
loo_rra <- function(drop_dataset) {
  sub <- tested_unfilt %>% filter(dataset != drop_dataset)
  if (nrow(sub) == 0) return(NULL)
  lists <- rank_contrast_lists(sub)
  if (length(lists) < 2) return(NULL)
  aggregateRanks(lists, N = n_distinct(sub$gene), exact = TRUE) %>%
    filter(Name %in% consensus_degs$gene) %>%
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

# genes that survive dropping either dataset
loo_robust <- loo_compare %>%
  filter(if_all(starts_with("survives_"), ~ .x %in% TRUE))
message(
  "Signature genes robust to dropping either dataset: ",
  nrow(loo_robust), " of ", nrow(loo_compare)
)
loo_robust %>% count(Tier) %>% print()


### performance of signature tiers vs validation genes ###

tiered_sig_list <- list(
  Current_RRA_meta = rra_meta_filt$gene_name,
  Tier1            = tier1$gene_name,
  Tier1_plus_2     = rra_meta %>%
    filter(Tier %in% c("Tier1_high_confidence", "Tier2_supported")) %>%
    pull(gene_name),
  RRA_only         = rra_results$gene_name,
  Meta_only        = meta_a$gene_name
)

tier_performance <- map_df(names(validation_sets), function(vs) {
  map_df(names(tiered_sig_list), function(nm) {
    calculate_sig_performance(
      tiered_sig_list[[nm]], validation_sets[[vs]], background_g
    ) %>%
      mutate(
        Signature = nm, Validation_set = vs,
        N_genes = length(unique(tiered_sig_list[[nm]]))
      )
  })
}) %>%
  select(Validation_set, Signature, N_genes, everything()) %>%
  arrange(Validation_set, desc(F1_Score))

print(tier_performance)
write_csv(tier_performance, paste0(outdir, "tier_performance_metrics.csv"))

# recall and F1 scale with signature size, so also compare precision (with a
# binomial CI) and enrichment over the background rate, overall and at matched N
precision_ci <- function(sig_genes, true_genes) {
  sig_genes <- unique(sig_genes)
  tp <- length(intersect(sig_genes, true_genes))
  bt <- binom.test(tp, length(sig_genes))
  tibble(
    n_genes = length(sig_genes), tp = tp,
    precision = unname(bt$estimate),
    ci_lb = bt$conf.int[1], ci_ub = bt$conf.int[2]
  )
}

# validation sets are already restricted to the background
base_rates <- lengths(validation_sets) / length(background_g)
message("Base rate of validation genes in background:")
print(round(base_rates, 4))

top_n_by <- function(df, order_col, n = 100) {
  df %>%
    arrange(.data[[order_col]]) %>%
    slice_head(n = n) %>%
    pull(gene_name)
}

precision_sig_list <- c(
  tiered_sig_list,
  list(
    RRA_top100   = top_n_by(rra_meta, "Score", 100),
    MetaP_top100 = top_n_by(filter(rra_meta, !is.na(meta_pval)), "meta_pval", 100)
  )
)

precision_results <- map_df(names(validation_sets), function(vs) {
  map_df(names(precision_sig_list), function(nm) {
    precision_ci(precision_sig_list[[nm]], validation_sets[[vs]]) %>%
      mutate(Signature = nm, Validation_set = vs)
  })
}) %>%
  mutate(enrichment_over_base = precision / base_rates[Validation_set]) %>%
  select(Validation_set, Signature, everything()) %>%
  arrange(Validation_set, desc(precision))

print(precision_results)
write_csv(precision_results, paste0(outdir, "signature_precision_metrics.csv"))


### DEG presence (1) / absence (0) per contrast for consensus genes ###
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
