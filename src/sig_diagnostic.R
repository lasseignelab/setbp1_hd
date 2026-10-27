suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(metafor)
  library(ComplexHeatmap)
  library(circlize)
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
resdir <- args$output_dir
filepath <- paste0(resdir, "signature_reversion/")
querypath <- paste0(resdir, "consensus_signature/")
outdir <- paste0(resdir, "diagnostics/")

# create dirs if needed
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

consensus_degs <- read_csv(paste0(querypath, "consensus_degs_unfiltered.csv"))

# ==============================================================================
# A. ARE THE CARDO CONTRASTS ANTI-CORRELATED?
# ==============================================================================
# In your output, essentially every printed gene flips sign between day 15 and
# day 21: LIN28A (-2.20 / +5.19), ABAT (+2.67 / -4.23), DNER (+2.82 / -4.50),
# EFNA5 (-1.86 / +2.79), PDK4 (+5.23 / -5.43), CELF4 (+2.00 / -5.72),
# CNTN1 (+1.10 / -4.91), THSD7A (+4.09 / -4.28), MYL9 (-0.78 / +5.50).
#
# A genome-wide negative correlation between two timepoints of the SAME
# KO-vs-WT comparison is not a biological pattern. Test it directly.

lfc_mat <- consensus_degs %>%
  select(gene, contrast, dataset, logFC) %>%
  unite("contrast_dataset", contrast, dataset) %>%
  distinct(gene, contrast_dataset, .keep_all = TRUE) %>%
  pivot_wider(names_from = contrast_dataset, values_from = logFC) %>%
  column_to_rownames("gene")

# Pairwise correlation on genes measured in both contrasts (NA, not 0)
lfc_cor <- cor(lfc_mat, use = "pairwise.complete.obs", method = "pearson")
print(round(lfc_cor, 2))

# Number of genes behind each pairwise correlation
n_pairs <- crossprod(!is.na(as.matrix(lfc_mat)))
print(n_pairs)

pdf(paste0(outdir, "contrast_logFC_correlation.pdf"), width = 8, height = 7)
Heatmap(
  lfc_cor,
  name = "Pearson r",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(sprintf("%.2f", lfc_cor[i, j]), x, y,
                    gp = grid::gpar(fontsize = 9))
  },
  column_title = "logFC correlation between contrasts (shared genes)"
)
dev.off()

# Decision rule: any pair of contrasts from the SAME dataset with r < -0.5 over
# hundreds of genes is a technical problem, not biology.


# ==============================================================================
# B. WHICH CARDO CONTRAST IS THE OUTLIER?
# ==============================================================================
# Use Shaw as the external reference. For genes DE in both datasets, ask how
# often each Cardo contrast agrees in sign with each Shaw contrast.

cardo_cols <- grep("_cardo$", colnames(lfc_mat), value = TRUE)
shaw_cols  <- grep("_shaw$",  colnames(lfc_mat), value = TRUE)

sign_agreement <- expand_grid(cardo = cardo_cols, shaw = shaw_cols) %>%
  rowwise() %>%
  mutate(
    n_shared = sum(!is.na(lfc_mat[[cardo]]) & !is.na(lfc_mat[[shaw]])),
    prop_same_sign = mean(
      sign(lfc_mat[[cardo]]) == sign(lfc_mat[[shaw]]),
      na.rm = TRUE
    )
  ) %>%
  ungroup()

print(sign_agreement)

# Expect ~0.5 under no relationship. Consistently < 0.4 for one Cardo contrast
# while the others are > 0.6 identifies the inverted or mislabelled one.

# Spot-check the genes the claims rest on:
check_genes <- c("LIN28A", "SPON1", "ROBO1", "PSMG3", "ISL2", "DRGX", "FOXG1")
consensus_degs %>%
  filter(gene_name %in% check_genes) %>%
  select(gene_name, dataset, contrast, logFC, adj.P.Val, se_logFC) %>%
  arrange(gene_name, dataset, contrast) %>%
  print(n = 50)


# ==============================================================================
# C. WHY THE META-ANALYSIS COLLAPSED
# ==============================================================================
# Your summary shows median I^2 = 98.49 and median QEp = 0. That is near-total
# between-contrast disagreement across the whole gene set, which is why
# meta_pval has a median of 0.699 and meta_padj a median of 1.0.
#
# Part of this is the SE disparity. Cardo SEs in your printout are 0.04-0.06;
# Shaw SEs are 0.25-0.36. Inverse-variance weighting gives Cardo roughly 35x
# the per-contrast weight, so the pooled estimate is driven by contrasts that
# contradict each other.

se_summary <- consensus_degs %>%
  group_by(dataset, contrast) %>%
  summarise(
    n_genes    = dplyr::n(),
    median_se  = median(se_logFC, na.rm = TRUE),
    median_lfc = median(abs(logFC), na.rm = TRUE),
    median_padj = median(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(relative_weight = (1 / median_se^2) / min(1 / median_se^2))

print(se_summary)

# Cardo adj.P.Val printing as exactly 0 alongside SEs of 0.04 is worth a look
# on its own — it suggests very low residual variance, which can happen when
# technical replicates are modelled as biological. Check the Cardo design:
#   table(cardo_metadata$genotype, cardo_metadata$timepoint, cardo_metadata$batch)


# ==============================================================================
# D. TWO-STAGE META-ANALYSIS  --  the fix for your actual concern
# ==============================================================================
# You switched to contrast-counting so that one bad timepoint would not nullify
# a whole dataset. That is the right concern, but contrast-counting makes the
# problem worse: a gene at -2.2 (d15) and +5.19 (d21) is counted as two
# independent supports when it is really one dataset contradicting itself.
#
# The standard answer is a hierarchical / two-stage pool:
#   Stage 1: pool contrasts WITHIN each dataset -> one estimate per dataset
#   Stage 2: pool the dataset-level estimates
#
# A noisy timepoint is absorbed at stage 1 (widening that dataset's CI) instead
# of either nullifying the dataset or masquerading as replication.

stage1 <- consensus_degs %>%
  group_by(gene, dataset) %>%
  group_modify(~ {
    if (nrow(.x) == 1) {
      tibble(ds_logFC = .x$logFC, ds_se = .x$se_logFC, ds_k = 1, ds_i2 = NA_real_)
    } else {
      fit <- tryCatch(rma(yi = .x$logFC, sei = .x$se_logFC, method = "REML"),
                      error = function(e) NULL)
      if (is.null(fit)) {
        tibble(ds_logFC = NA_real_, ds_se = NA_real_, ds_k = nrow(.x), ds_i2 = NA_real_)
      } else {
        tibble(ds_logFC = as.numeric(fit$b), ds_se = as.numeric(fit$se),
               ds_k = fit$k, ds_i2 = as.numeric(fit$I2))
      }
    }
  }) %>%
  ungroup()

write_csv(stage1, paste0(outdir, "stage1_within_dataset_meta.csv"))

stage2 <- stage1 %>%
  filter(!is.na(ds_logFC), !is.na(ds_se), ds_se > 0) %>%
  group_by(gene) %>%
  filter(dplyr::n() >= 2) %>%          # requires BOTH datasets
  group_modify(~ {
    fit <- tryCatch(rma(yi = .x$ds_logFC, sei = .x$ds_se, method = "REML"),
                    error = function(e) NULL)
    if (is.null(fit)) return(tibble())
    tibble(
      h_logFC = as.numeric(fit$b),  h_pval  = as.numeric(fit$pval),
      h_ci_lb = as.numeric(fit$ci.lb), h_ci_ub = as.numeric(fit$ci.ub),
      h_i2    = as.numeric(fit$I2), h_qep   = as.numeric(fit$QEp)
    )
  }) %>%
  ungroup() %>%
  mutate(h_padj = p.adjust(h_pval, method = "BH"))

summary(stage2)
write_csv(stage2, paste0(outdir, "stage2_hierarchical_meta.csv"))

# Compare: does pooling within-dataset first rescue the heterogeneity problem?
message("Median I^2, flat meta:         ", round(median(meta_stats$meta_i2, na.rm = TRUE), 1))
message("Median I^2, within-dataset:    ", round(median(stage1$ds_i2, na.rm = TRUE), 1))
message("Median I^2, across datasets:   ", round(median(stage2$h_i2, na.rm = TRUE), 1))
message("Genes with h_padj < 0.05:      ", sum(stage2$h_padj < 0.05, na.rm = TRUE))

hier_signature <- rra %>%
  inner_join(stage2, by = c("Name" = "gene")) %>%
  left_join(distinct(consensus_degs[, c("gene", "gene_name")]), by = c("Name" = "gene")) %>%
  filter(Score < 0.05, h_padj < 0.05, sign(h_ci_lb) == sign(h_ci_ub)) %>%
  mutate(Direction = if_else(h_logFC > 0, "Upregulated", "Downregulated"))

write_csv(hier_signature, paste0(outdir, "hierarchical_consensus_signature.csv"))
message("Hierarchical signature size: ", nrow(hier_signature))


# ==============================================================================
# E. LEAVE-ONE-OUT INTERSECTION
# ==============================================================================
# You printed 318 surviving without Cardo and 687 without Shaw, but not the
# intersection. Genes robust to dropping EITHER dataset are the real consensus.

loo_robust <- loo_compare %>%
  filter(if_all(starts_with("survives_"), ~ .x %in% TRUE))

message("Robust to dropping either dataset: ", nrow(loo_robust), " of ", nrow(loo_compare))
write_csv(loo_robust, paste0(outdir, "loo_robust_genes.csv"))

loo_robust %>% count(Tier) %>% print()


# ==============================================================================
# F. A FAIR SIGNATURE COMPARISON
# ==============================================================================
# F1 and recall are confounded by signature size here: recall is TP / 2081, so a
# 19-gene signature cannot exceed 0.009 by construction. Comparing F1 across
# sets of 944 and 19 genes measures size, not quality.
#
# Two fixes: report precision with a confidence interval, and compare precision
# at MATCHED N by taking the top N from each ranking.

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

base_rate <- length(intersect(true_positives, background_g)) / length(background_g)
message("Base rate of validation genes in background: ", round(base_rate, 4))

map_df(names(tiered_sig_list), function(nm) {
  precision_ci(tiered_sig_list[[nm]], true_positives) %>% mutate(Signature = nm)
}) %>%
  mutate(enrichment_over_base = precision / base_rate) %>%
  select(Signature, everything()) %>%
  arrange(desc(precision)) %>%
  print()

# Matched-N comparison: top 100 by each ranking criterion
top_n_by <- function(df, order_col, n = 100) {
  df %>% arrange(.data[[order_col]]) %>% slice_head(n = n) %>% pull(gene_name)
}

matched <- list(
  RRA_top100        = top_n_by(rra_meta, "Score", 100),
  MetaP_top100      = top_n_by(filter(rra_meta, !is.na(meta_pval)), "meta_pval", 100),
  Hierarchical_top100 = hier_signature %>% arrange(h_padj) %>% slice_head(n = 100) %>% pull(gene_name)
)

map_df(names(matched), function(nm) {
  precision_ci(matched[[nm]], true_positives) %>% mutate(Signature = nm)
}) %>%
  mutate(enrichment_over_base = precision / base_rate) %>%
  print()

# NOTE on the validation set itself: true_positives is dominated by membership
# in six broad pathways (1,777 gene-pathway pairs). That rewards signatures
# containing many generically annotated signalling genes. It is a reasonable
# sanity check but it cannot adjudicate between a permissive and a stringent
# filter, because stringency and pathway-annotation density are unrelated.
# The SGS reciprocity test is a better discriminator.

#### END ####
# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
