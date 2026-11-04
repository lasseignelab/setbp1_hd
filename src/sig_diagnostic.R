# script to run contrast-level diagnostics on the DEGs used for the consensus
# signature. Sections A-C read the DEA results directly and include ALL
# candidate contrasts (including Cardo day 21, KOvsWT_21), which is the
# evidence for excluding KOvsWT_21 in 10_consensus_signature.R.
# Section D runs an alternative two-stage (within- then across-dataset)
# meta-analysis on the final consensus DEGs written by 10_consensus_signature.R.
set.seed(42)

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(metafor)
  library(ComplexHeatmap)
  library(circlize)
})

###### SET UP ######
parser <- ArgumentParser()

parser$add_argument("-d", "--degs",
  type = "character",
  metavar = "PATH",
  help = "Parent directory containing each dataset's DEG CSV files"
)
parser$add_argument("-c", "--consensus_dir",
  type = "character",
  metavar = "PATH",
  help = "Directory with 10_consensus_signature.R outputs"
)
parser$add_argument("-o", "--outdir",
  type = "character",
  metavar = "PATH",
  help = "Path to the directory to write diagnostic results"
)

args <- parser$parse_args()
degs <- args$degs
consensus_dir <- args$consensus_dir
outdir <- args$outdir

if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# start timer
ptm <- proc.time()

read_deg <- function(dataset_name, type = "significant_degs.csv") {
  file <- file.path(degs, dataset_name, type)
  if (!file.exists(file)) stop("Missing file: ", file)
  read_csv(file, show_col_types = FALSE) %>% mutate(dataset = dataset_name)
}

###### READ IN DATA ######
# same redundant/biased contrasts skipped as in 10_consensus_signature.R, but
# KOvsWT_21 is deliberately KEPT here so it can be evaluated
contrasts_skip <- c("KOvsWT_NPC", "KOvsWT_all", "NPC_HD", "PATH_vs_PATH")

all_degs <- bind_rows(read_deg("cardo"), read_deg("shaw")) %>%
  filter(!contrast %in% contrasts_skip)

# genes that are DEGs in more than one contrast (as in the consensus script)
candidate_degs <- all_degs %>%
  group_by(gene) %>%
  filter(n_distinct(contrast) > 1) %>%
  ungroup()

table(candidate_degs$dataset, candidate_degs$contrast)

###### A. logFC CORRELATION BETWEEN CONTRASTS ######
# contrasts of the same KO-vs-WT comparison at different timepoints should be
# positively correlated; a strong negative correlation indicates a problem
lfc_mat <- candidate_degs %>%
  select(gene, contrast, dataset, logFC) %>%
  unite("contrast_dataset", contrast, dataset) %>%
  distinct(gene, contrast_dataset, .keep_all = TRUE) %>%
  pivot_wider(names_from = contrast_dataset, values_from = logFC) %>%
  column_to_rownames("gene")

# pairwise correlation on genes that are DEGs in both contrasts
lfc_cor <- cor(lfc_mat, use = "pairwise.complete.obs", method = "pearson")
print(round(lfc_cor, 2))

# number of shared genes behind each pairwise correlation
n_pairs <- crossprod(!is.na(as.matrix(lfc_mat)))
print(n_pairs)

write_csv(
  as_tibble(lfc_cor, rownames = "contrast"),
  paste0(outdir, "contrast_logFC_correlation.csv")
)
write_csv(
  as_tibble(n_pairs, rownames = "contrast"),
  paste0(outdir, "contrast_shared_gene_counts.csv")
)

pdf(paste0(outdir, "contrast_logFC_correlation.pdf"), width = 8, height = 7)
draw(Heatmap(
  lfc_cor,
  name = "Pearson r",
  col = colorRamp2(c(-1, 0, 1), c("blue", "white", "red")),
  cell_fun = function(j, i, x, y, w, h, fill) {
    grid::grid.text(
      sprintf("%.2f (n=%d)", lfc_cor[i, j], n_pairs[i, j]), x, y,
      gp = grid::gpar(fontsize = 7)
    )
  },
  column_title = "logFC correlation between contrasts (shared DEGs)"
))
dev.off()

###### B. SIGN AGREEMENT OF EACH CARDO CONTRAST WITH SHAW ######
# Shaw as an external reference: for genes that are DEGs in both, how often
# does each Cardo contrast agree in direction with each Shaw contrast?
# (~0.5 expected under no relationship)
cardo_cols <- grep("_cardo$", colnames(lfc_mat), value = TRUE)
shaw_cols <- grep("_shaw$", colnames(lfc_mat), value = TRUE)

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
write_csv(sign_agreement, paste0(outdir, "cardo_vs_shaw_sign_agreement.csv"))

# per-contrast logFC for key genes
check_genes <- c("LIN28A", "SPON1", "ROBO1", "PSMG3", "ISL2", "DRGX", "FOXG1")
candidate_degs %>%
  filter(gene_name %in% check_genes) %>%
  select(gene_name, dataset, contrast, logFC, adj.P.Val, se_logFC) %>%
  arrange(gene_name, dataset, contrast) %>%
  print(n = 50)

###### C. STANDARD ERRORS AND META-ANALYSIS WEIGHTS BY CONTRAST ######
# inverse-variance weighting gives contrasts with small SEs more weight in the
# per-gene meta-analysis; summarise how unequal those weights are
se_summary <- candidate_degs %>%
  group_by(dataset, contrast) %>%
  summarise(
    n_genes     = dplyr::n(),
    median_se   = median(se_logFC, na.rm = TRUE),
    median_lfc  = median(abs(logFC), na.rm = TRUE),
    median_padj = median(adj.P.Val, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(relative_weight = (1 / median_se^2) / min(1 / median_se^2))

print(se_summary)
write_csv(se_summary, paste0(outdir, "contrast_se_summary.csv"))

###### D. TWO-STAGE META-ANALYSIS (FINAL CONSENSUS DEGS) ######
# Stage 1: pool contrasts within each dataset -> one estimate per dataset
# Stage 2: pool dataset-level estimates (genes that are DEGs in both datasets)
# Compared against the single-stage meta-analysis in 10_consensus_signature.R
consensus_degs <- read_csv(
  paste0(consensus_dir, "consensus_degs_unfiltered.csv"),
  show_col_types = FALSE
)
meta_stats <- read_csv(
  paste0(consensus_dir, "metaanalysis_degs.csv"),
  show_col_types = FALSE
)
rra <- read_csv(
  paste0(consensus_dir, "rankaggregate_degs.csv"),
  show_col_types = FALSE
)

stage1 <- consensus_degs %>%
  group_by(gene, dataset) %>%
  group_modify(~ {
    if (nrow(.x) == 1) {
      tibble(
        ds_logFC = .x$logFC, ds_se = .x$se_logFC, ds_k = 1, ds_i2 = NA_real_
      )
    } else {
      fit <- tryCatch(
        rma(yi = .x$logFC, sei = .x$se_logFC, method = "REML"),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        tibble(
          ds_logFC = NA_real_, ds_se = NA_real_,
          ds_k = nrow(.x), ds_i2 = NA_real_
        )
      } else {
        tibble(
          ds_logFC = as.numeric(fit$b), ds_se = as.numeric(fit$se),
          ds_k = fit$k, ds_i2 = as.numeric(fit$I2)
        )
      }
    }
  }) %>%
  ungroup()

write_csv(stage1, paste0(outdir, "stage1_within_dataset_meta.csv"))

stage2 <- stage1 %>%
  filter(!is.na(ds_logFC), !is.na(ds_se), ds_se > 0) %>%
  group_by(gene) %>%
  filter(dplyr::n() >= 2) %>%
  group_modify(~ {
    fit <- tryCatch(
      rma(yi = .x$ds_logFC, sei = .x$ds_se, method = "REML"),
      error = function(e) NULL
    )
    if (is.null(fit)) return(tibble())
    tibble(
      h_logFC = as.numeric(fit$b), h_pval = as.numeric(fit$pval),
      h_ci_lb = as.numeric(fit$ci.lb), h_ci_ub = as.numeric(fit$ci.ub),
      h_i2 = as.numeric(fit$I2), h_qep = as.numeric(fit$QEp)
    )
  }) %>%
  ungroup() %>%
  mutate(h_padj = p.adjust(h_pval, method = "BH"))

summary(stage2)
write_csv(stage2, paste0(outdir, "stage2_hierarchical_meta.csv"))

message(
  "Median I^2, single-stage meta: ",
  round(median(meta_stats$meta_i2, na.rm = TRUE), 1)
)
message(
  "Median I^2, within-dataset:    ",
  round(median(stage1$ds_i2, na.rm = TRUE), 1)
)
message(
  "Median I^2, across datasets:   ",
  round(median(stage2$h_i2, na.rm = TRUE), 1)
)
message(
  "Genes with h_padj < 0.05:      ",
  sum(stage2$h_padj < 0.05, na.rm = TRUE)
)

hier_signature <- rra %>%
  inner_join(stage2, by = c("Name" = "gene")) %>%
  left_join(
    distinct(consensus_degs[, c("gene", "gene_name")]),
    by = c("Name" = "gene")
  ) %>%
  filter(Score < 0.05, h_padj < 0.05, sign(h_ci_lb) == sign(h_ci_ub)) %>%
  mutate(Direction = if_else(h_logFC > 0, "Upregulated", "Downregulated"))

write_csv(
  hier_signature, paste0(outdir, "hierarchical_consensus_signature.csv")
)
message("Hierarchical signature size: ", nrow(hier_signature))

#### END ####
# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
