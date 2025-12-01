suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(UpSetR)
  library(ComplexHeatmap)
  library(statmod)
  library(gprofiler2)
  library(RobustRankAggreg)
  library(metafor)
})

# create parser object
parser <- ArgumentParser()

parser$add_argument("-d", "--degs",
  type = "character",
  metavar = "PATH",
  help = "Parent directory containing each dataset's DEG CSV files" # nolint
)
parser$add_argument("-o", "--outdir",
  type = "character",
  metavar = "PATH",
  help = "Path to the directory to write consensus signature results"
)
# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# create output dir if needed
if (!dir.exists(args$outdir)) dir.create(args$outdir, recursive = TRUE)

# helper function to read in DEGs
read_deg <- function(dataset_name) {
  file <- file.path(args$degs, dataset_name, "combined_contrasts_sig_degs.csv")
  if (!file.exists(file)) stop("Missing file: ", file)
  read_csv(file) %>% mutate(dataset = dataset_name)
}

# start timer
ptm <- proc.time()

#### Set variables for dataset ####
cardo_degs <- read_deg("cardo")
shaw_degs <- read_deg("shaw")
tanaka_degs <- bind_rows(read_deg("tanaka1"), read_deg("tanaka2"))

unique(cardo_degs$contrast)
unique(shaw_degs$contrast)

# contrasts to keep
contrasts_keep <- c(
  "KOvsWT_15", "KOvsWT_21", "KOvsWT_34", "VUSvsWT", "PATH2vsWT", "PATH3vsWT"
)

# combine human data sets and filter
combined_degs <- bind_rows(cardo_degs, shaw_degs)

consensus_genes <- intersect(cardo_degs$gene, shaw_degs$gene)
consensus_degs <- combined_degs %>%
  filter(contrast %in% contrasts_keep) %>%
  filter(gene %in% consensus_genes)

# save consensus DEGs
write_csv(consensus_degs, paste0(args$outdir, "consensus_degs_unfiltered.csv"))

# convert DEGs to ENTREZ -- required for signaturesearch downstream
entrez_degs <- gconvert(
  query = consensus_degs$gene, target = "ENTREZGENE_ACC"
) %>%
  select(gene = input, entrez = target)
# long format for plotting
consensus_long <- consensus_degs %>%
  select(gene, gene_name, contrast, dataset, logFC, adj.P.Val, se_logFC) %>%
  pivot_longer(
    cols = c(logFC, adj.P.Val, se_logFC),
    names_to = "metric",
    values_to = "value"
  )

consensus_wide <- consensus_long %>%
  distinct() %>%
  unite("contrast_dataset_metric", contrast, dataset, metric) %>%
  pivot_wider(names_from = contrast_dataset_metric, values_from = value)
print(paste0("Consensus genes dataframe: "))
str(consensus_wide)
print(paste0("Number of consensus genes: "))
length(unique(consensus_long$gene))

### COMPARE DEGs ACROSS CONTRASTS + DATA SETS ###
# compare to Tanaka et al DEGk
tanaka_genes <- gorth(
  query = tanaka_degs$gene,
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)

tanaka_anno <- tanaka_genes %>%
  select(gene_name = ortholog_name) %>%
  mutate(contrast = "Tanaka_DEG")

compare_contrasts <- consensus_long %>%
  select(gene_name, contrast) %>%
  bind_rows(tanaka_anno) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = contrast, values_from = value, values_fill = 0
  ) %>%
  column_to_rownames(var = "gene_name")

str(compare_contrasts)

# set colors for upset plot
contrast_dataset_map <- consensus_long %>%
  distinct(contrast, dataset) %>%
  bind_rows(
    data.frame(
      dataset = c("Tanaka"),
      contrast = c("Tanaka_DEG")
    )
  )

dataset_colors <- c(
  "cardo" = "orange",
  "shaw" = "steelblue",
  "Tanaka" = "gray90"
)

png(
  paste0(
    args$outdir, "upset_plot_compare_shaw_cardo.png"
  ),
  width = 30, height = 20, units = "cm", res = 300
)
upset(compare_contrasts,
  set.metadata = list(data = contrast_dataset_map, plots = list(
    list(
      type = "matrix_rows",
      column = "dataset",
      assign = 10,
      colors = dataset_colors
    )
  ))
)
dev.off()

png(
  paste0(
    args$outdir, "upset_plot_compare_contrasts.png"
  ),
  width = 30, height = 20, units = "cm", res = 300
)
upset(compare_contrasts, sets = colnames(compare_contrasts), order.by = "freq")
dev.off()

### PRIORTIZE/NARROW DEGS TO MAKE SIGNATURE ###
# opt 1) high frequency DEGs
print(paste0("Finding high-frequency ('majority vote') DEGs..."))
freq_degs <- consensus_long %>%
  select(!c(metric, value)) %>%
  distinct() %>%
  add_count(gene) %>%
  left_join(entrez_degs, by = "gene")

top_freq_degs <- freq_degs %>%
  mutate(method = "majority") %>%
  filter(!is.na(entrez)) %>%
  filter(n > 3)
str(top_freq_degs)
length(unique(top_freq_degs$gene))

# save majority/high freq signature
write_csv(
  top_freq_degs, paste0(args$outdir, "majorityfreq_consensus_signature.csv")
)

# opt 2) meta-analysis per gene
print(paste0("Finding effect size meta analysis DEGs..."))
meta_results <- consensus_degs %>%
  add_count(gene) %>%
  filter(n > 1) %>%
  group_by(gene) %>%
  nest() %>%
  mutate(meta = map(data, ~ tryCatch(
    rma(yi = .x$logFC, sei = .x$se_logFC, method = "REML"),
    error = function(e) NULL
  ))) %>%
  mutate(
    meta_logFC = map_dbl(
      meta, ~ if (!is.null(.x)) as.numeric(.x$b) else NA_real_
    ),
    meta_pval = map_dbl(
      meta, ~ if (!is.null(.x)) as.numeric(.x$pval) else NA_real_
    ),
    meta_ci_lb = map_dbl(
      meta, ~ if (!is.null(.x)) as.numeric(.x$ci.lb) else NA_real_
    ),
    meta_ci_ub = map_dbl(
      meta, ~ if (!is.null(.x)) as.numeric(.x$ci.ub) else NA_real_
    )
  ) %>%
  select(-data)

meta_stats <- as.data.frame(meta_results[, -2]) %>%
  left_join(entrez_degs, by = "gene") %>%
  distinct()
str(meta_stats)
summary(meta_stats$meta_pval)

meta_p <- meta_stats %>%
  filter(!is.na(entrez)) %>%
  slice_min(meta_pval, n = 300) %>%
  mutate(method = "meta_p")

# save meta analysis & signature
write_csv(meta_stats, paste0(args$outdir, "metaanalysis_degs.csv"))
write_csv(meta_p, paste0(args$outdir, "metaanalysis_consensus_signature.csv"))

# opt 3) robust rang aggregate -- create ranked lists per contrast
print(paste0("Finding robust rank aggregate (RRA) DEGs..."))
ranked_lists <- lapply(unique(consensus_degs$contrast), function(c) {
  consensus_degs %>%
    filter(contrast == c) %>%
    arrange(abs(adj.P.Val)) %>%
    pull(gene)
})
head(ranked_lists)
# run RRA
rra <- aggregateRanks(ranked_lists, exact = TRUE)
str(rra)
summary(rra)

rra_results <- rra %>%
  left_join(entrez_degs, by = c("Name" = "gene")) %>%
  distinct() %>%
  filter(!is.na(entrez)) %>%
  slice_min(order_by = Score, n = 300) %>%
  mutate(method = "RRA") %>%
  select(gene = Name, entrez, method)

# save RRA & signature
write_csv(rra, paste0(args$outdir, "rankaggregate_degs.csv"))
write_csv(
  rra_results, paste0(args$outdir, "rankaggregate_consensus_signature.csv")
)

## Composite Scores from RRA, Meta Analysis, and Freq-Weight ##
print(
  paste0(
    "Finding composite scores by RRA, meta analysis, and frequency weight..."
  )
)
composite_scores <- top_freq_degs %>%
  select(gene, gene_name, n) %>%
  distinct() %>%
  inner_join(meta_stats, by = "gene") %>%
  inner_join(rra, by = c("gene" = "Name")) %>%
  mutate(
    freq_weight = log2(n),
    rank_pval = Score,
    composite_score = -log10(meta_pval) + -log10(rank_pval) + freq_weight
  )
str(composite_scores)
length(unique(composite_scores$gene))

# save scores across methods for DEGs + composite scores
write_csv(
  composite_scores, paste0(args$outdir, "compositescore_full_table_degs.csv")
)

# top 200 by composite score
composite_signature <- composite_scores %>%
  filter(!is.na(entrez)) %>%
  slice_max(composite_score, n = 300)
composite_signature

# save composite score signature
write_csv(
  composite_signature, paste0(
    args$outdir, "compositescore_consensus_signature.csv"
  )
)

## compare overlap between options for finding signature ##
compare_sigs_up <- top_freq_degs %>%
  bind_rows(
    rra_results, meta_p, rra_results, composite_signature
  ) %>%
  select(gene, method) %>%
  mutate(method = ifelse(is.na(method), "Composite", method)) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = method, values_from = value, values_fill = 0) %>%
  column_to_rownames(var = "gene")

str(compare_sigs_up)

png(
  paste0(
    args$outdir, "upset_plot_compare_sigs.png"
  ),
  width = 30, height = 20, units = "cm", res = 300
)
upset(compare_sigs_up, sets = colnames(compare_sigs_up), order.by = "freq")
dev.off()

# up vs down DEGs in signature
composite_up <- composite_signature %>% filter(meta_logFC > 0)
composite_down <- composite_signature %>% filter(meta_logFC < 0)
print(paste0(length(unique(composite_up$gene)), " UP composite DEGs"))
print(paste0(length(unique(composite_down$gene)), " DOWN composite DEGs"))

verification_degs <- tanaka_degs %>%
  inner_join(tanaka_genes, by = c("gene" = "input_ensg")) %>%
  select(name = ortholog_name, target = ortholog_ensg, logFC, adj.P.Val) %>%
  distinct() %>%
  mutate(
    Composite = ifelse(target %in% composite_signature$gene, TRUE, FALSE),
    RRA = ifelse(target %in% rra_results$gene, TRUE, FALSE),
    Meta = ifelse(target %in% meta_p$gene, TRUE, FALSE),
    Freq = ifelse(target %in% top_freq_degs$gene, TRUE, FALSE)
  ) %>%
  filter(if_any(Composite:Freq, ~ . == TRUE))

print(paste0("Verification DEGs: "))
verification_degs

#### ORA of Consensus DEGs ####
# function for ORA plotting
run_ora_plot <- function(genes, title) {
  ora_res <- tryCatch(
    {
      gost(genes, sources = c("GO", "KEGG", "WP"))$result
    },
    error = function(e) {
      message(paste("ORA failed for", title, ":", e$message))
      NULL
    }
  )
  # skip plotting for no ORA results
  if (is.null(ora_res) || nrow(ora_res) == 0) {
    message(paste("No ORA results for", title, "-- skipping plot"))
    NULL
  } else {

    # filter and select top terms
    plot_data <- ora_res %>%
      filter(term_size < 750) %>%
      slice_max(recall, n = 40)

    # create plot
    plot <- plot_data %>%
      ggplot(
        aes(x = recall, y = term_name,
        size = intersection_size, colour = precision)
      ) +
      geom_point() +
      labs(title = title, x = "Recall", y = NULL) +
      theme_minimal()

    list(result = ora_res, plot = plot)
  }
}
# up vs down meta ORA
meta_up <- meta_p %>% filter(meta_logFC > 0)
meta_down <- meta_p %>% filter(meta_logFC < 0)

# run ORA for each signature
sig_list <- list(composite_signature, rra_results, meta_p, top_freq_degs)
ora_list <- lapply(seq_along(sig_list), function(sig) {
  run_ora_plot(sig_list[[sig]]$gene, sig_list[[sig]]$method)
})

names(ora_list) <- c("composite", "RRA", "meta", "freq")
res_list <- lapply(seq_along(ora_list), function(ora) {
  ora_list[[ora]]$result %>%
    mutate(method = paste0(names(ora_list)[ora]))
})

# simplify results for term names in > 1 sig ORA results
res_df <- bind_rows(res_list)
str(res_df)

terms_filtered <- res_df %>%
  group_by(term_name) %>%
  mutate(term_count = n_distinct(method)) %>%
  filter(term_count > 1) %>%
  ungroup()

# plot
terms_plot <- terms_filtered %>%
  filter(term_size < 1750) %>%
  ggplot(
    aes(
      y = reorder(term_name, term_size), x = -log10(p_value),
      size = recall, colour = method
    )
  ) +
  geom_point(alpha = 0.6)
ggsave(
  paste0(
    args$outdir, "ORA_combined_terms.png"
  ), width = 20, height = 28, units = "cm")

ora_plots <- lapply(seq_along(ora_list), function(ora) {
  ora_list[[ora]]$plot
})

# save all plots to one PDF
pdf(paste0(args$outdir, "ORA_all_plots.pdf"), width = 10, height = 8)
for (plot in ora_plots) {
  print(plot)
}
dev.off()

#### END ####

# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
