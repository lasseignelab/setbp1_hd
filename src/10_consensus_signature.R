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
  help = "RSE file (rds format) for EDA from recount3 data (00_data_download or 01_data_preprocessing)" # nolint
)
parser$add_argument("-s", "--sra",
  type = "character",
  metavar = "PATH",
  help = "SRA or dataset for naming purposes" # nolint
)
parser$add_argument("-o", "--output_dir",
  type = "character",
  metavar = "PATH",
  help = "Path to the directory to create and write EDA results"
)

# get command line options, if help option encountered print help and exit
args <- parser$parse_args()

# source functions file if supplied
if (!is.null(args$functions)) source(args$functions)

# start timer
ptm <- proc.time()

#### Set variables for dataset ####
print(paste0("Reading in RSE from: ", args$rse))

degs <- read_csv(args$rse)

# create output dir if needed
filepath <- paste0(args$output_dir, sra, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
# save labeled and sig-filtered results
# degs_comb <- read_csv(, file = paste0(filepath, "combined_contrasts_labeled_degs.csv"))
degs <- "./results/dea/"
filepath <- "./results/dea/consensus_signature/"
cardo_degs <- read_csv(paste0(degs, "cardo/combined_contrasts_sig_degs.csv")) %>%
  mutate(dataset = "cardo")
shaw_degs <- read_csv(paste0(degs, "shaw/combined_contrasts_sig_degs.csv")) %>%
  mutate(dataset = "shaw")
tanaka1_degs <- read_csv(paste0(degs, "tanaka1/combined_contrasts_sig_degs.csv"))
tanaka2_degs <- read_csv(paste0(degs, "tanaka2/combined_contrasts_sig_degs.csv"))
tanaka_degs <- bind_rows(tanaka1_degs, tanaka2_degs)
setbp1_targets <- read_csv("./data/241203_setbp1_targets_names_only.csv", col_names = FALSE)

unique(cardo_degs$contrast)
unique(shaw_degs$contrast)

# contrasts to keep
contrasts_keep <- c("KOvsWT_15", "KOvsWT_21", "KOvsWT_34", "VUSvsWT", "PATH2vsWT", "PATH3vsWT")

# combine data sets
combined_degs <- bind_rows(cardo_degs, shaw_degs)

# consensus DEGs
consensus_genes <- intersect(cardo_degs$gene, shaw_degs$gene)
consensus_degs <- combined_degs %>%
  filter(contrast %in% contrasts_keep) %>%
  filter(gene %in% consensus_genes)
str(consensus_degs)

# compare to SETBP1 targets
intersect(consensus_degs$gene_name, setbp1_targets$X1)

# save consensus DEGs
write_csv(consensus_degs, paste0(filepath, "consensus_degs_unfiltered.csv"))

# long format for plotting
consensus_long <- consensus_degs %>%
  select(gene, gene_name, contrast, dataset, logFC, adj.P.Val, se_logFC) %>%
  pivot_longer(cols = c(logFC, adj.P.Val, se_logFC),
               names_to = "metric",
               values_to = "value") 
str(consensus_long)
table(consensus_long$contrast)

consensus_wide <- consensus_long %>%
  unite("contrast_dataset_metric", contrast, dataset, metric) %>%
  pivot_wider(names_from = contrast_dataset_metric, values_from = value)
str(consensus_wide)
length(unique(consensus_long$gene))

### COMPARE DEGs ACROSS CONTRASTS + DATA SETS ###
# compare to Tanaka et al DEGk
tanaka_genes <- gorth(query = tanaka_degs$gene, source_organism = "mmusculus", target_organism = "hsapiens")

tanaka_targets <- setbp1_targets %>%
  select(gene_name = X1) %>%
  mutate(contrast = "SETBP1_target") %>%
  bind_rows(data.frame(gene_name = tanaka_genes$ortholog_name, contrast = "Tanaka_DEG"))

compare_contrasts <- consensus_long %>%
  select(gene_name, contrast) %>%
  bind_rows(tanaka_targets) %>%
    distinct() %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = contrast, values_from = value, values_fill = 0) %>%
  column_to_rownames(var = "gene_name")

str(compare_contrasts)

# set colors for upset plot
contrast_dataset_map <- consensus_long %>%
  distinct(contrast, dataset) %>%
  bind_rows(
    data.frame(
      dataset = c("SETBP1", "Tanaka"), contrast = c("SETBP1_target", "Tanaka_DEG")))

dataset_colors <- c("cardo" = "orange",
"shaw" = "steelblue",
"SETBP1" = "gray60",
"Tanaka" = "gray90")

# contrast_colors <- ifelse(
#   contrasts$dataset == "cardo", "#1f77b4", ifelse(
#     contrasts$dataset == "shaw", "#ff7f0e", "black"))

png(paste0(filepath, "upset_plot_compare_shaw_cardo.png"), width = 30, height = 20, units = "cm", res = 300)
upset(compare_contrasts, 
  set.metadata = list(data = contrast_dataset_map, plots = list(
    list(type = "matrix_rows", column = "dataset", assign = 10, colors = dataset_colors)
  )))
dev.off()

png(paste0(filepath, "upset_plot_compare_contrasts.png"), width = 30, height = 20, units = "cm", res = 300)
upset(compare_contrasts, sets = colnames(compare_contrasts), order.by = "freq")
dev.off()


### PRIORTIZE/NARROW DEGS TO MAKE SIGNATURE ###
# opt 1) high frequency DEGs
freq_degs <- consensus_long %>%
  select(!c(metric, value)) %>%
  distinct() %>%
  add_count(gene)

top_freq_degs <- freq_degs %>%
  mutate(method = "majority") %>%
  filter(n > 4)
str(top_freq_degs)
length(unique(top_freq_degs$gene))

# opt 3) meta-analysis per gene
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
    meta_logFC = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$b) else NA_real_),
    meta_pval   = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$pval) else NA_real_),
    meta_ci_lb  = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$ci.lb) else NA_real_),
    meta_ci_ub  = map_dbl(meta, ~ if (!is.null(.x)) as.numeric(.x$ci.ub) else NA_real_)
  ) %>%
  select(-data)

meta_stats <- as.data.frame(meta_results[,-2])
str(meta_stats)
summary(meta_stats$meta_pval)

meta_lfc <- meta_stats %>%
  slice_max(abs(meta_logFC), n = 200) %>%
  mutate(method = "meta_lfc")

meta_p <- meta_stats %>%
  slice_min(meta_pval, n = 200) %>%
  mutate(method = "meta_p")

# opt 4) robust rang aggregate -- create ranked lists per contrast
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
  slice_min(order_by = Score, n = 200) %>%
  mutate(method = "RRA") %>%
  select(gene = Name, method)

## Composite Scores from RRA, Meta Analysis, and Freq-Weight ##
composite_scores <- top_freq_degs %>%
  select(gene, gene_name, n) %>%
  distinct() %>%
  inner_join(meta_stats, by = "gene") %>%
  inner_join(rra, by = c("gene" = "Name")) %>%
  mutate(
  freq_weight = log2(n),
  rank_pval = Score,
  composite_score = -log10(meta_pval) + -log10(rank_pval) + freq_weight)
str(composite_scores)

# save scores across methods for DEGs + composite scores
write_csv(composite_scores, paste0(filepath, "full_scored_degs.csv"))

# top 200 by composite score
composite_signature <- composite_scores %>%
  slice_max(composite_score, n = 200)
composite_signature

# save composite score signature
write_csv(composite_signature, paste0(filepath, "consensus_signature.csv"))

 
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

png(paste0(filepath, "upset_plot_compare_sigs.png"), width = 30, height = 20, units = "cm", res = 300)
upset(compare_sigs_up, sets = colnames(compare_sigs_up), order.by = "freq")
dev.off()

# up vs down DEGs in signature
composite_up <- composite_signature %>% filter(meta_logFC > 0)
composite_up
composite_down <- composite_signature %>% filter(meta_logFC < 0)
composite_down

# # compare to Tanaka et al DEGk
# tanaka_genes <- gorth(query = tanaka_degs$gene, source_organism = "mmusculus", target_organism = "hsapiens")
#
# compare to SETBP1 targets
setbp1_genes <- gconvert(query = setbp1_targets$X1)

verification_degs <- tanaka_degs %>%
  inner_join(tanaka_genes, by = c("gene" = "input_ensg")) %>%
  select(name = ortholog_name, target = ortholog_ensg, logFC, adj.P.Val) %>%
  distinct() %>%
  bind_rows(setbp1_genes[,4:5]) %>%
  mutate(
    Composite = ifelse(target %in% composite_signature$gene, TRUE, FALSE),
    RRA = ifelse(target %in% rra_results$gene, TRUE, FALSE),
    Meta = ifelse(target %in% meta_p$gene, TRUE, FALSE),
    Freq = ifelse(target %in% top_freq_degs$gene, TRUE, FALSE)
  ) %>%
  filter(if_any(Composite:Freq, ~ . == TRUE))

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
    # print all terms
    print(ora_res$term_name)

    # filter and select top terms
    plot_data <- ora_res %>%
      filter(term_size < 750) %>%
      slice_max(recall, n = 40)

    # create plot
    plot <- plot_data %>%
      ggplot(aes(
        x = recall, y = term_name, size = intersection_size, colour = precision # nolint
      )) +
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
ora_list <- lapply(seq_along(sig_list), function(i) {
  run_ora_plot(sig_list[[i]]$gene, sig_list[[i]]$method)
})

# ora_list <- list(ora_sig, ora_sig_up, ora_sig_down, ora_rank, ora_meta, ora_meta_up, ora_meta_down)
names(ora_list) <- c("composite", "RRA", "meta", "freq")
res_list <- lapply(seq_along(ora_list), function(i) {
                   ora_list[[i]]$result %>%
                   mutate(method = paste0(names(ora_list)[i]))
                   })
str(res_list)

# simplify results for term names in > 1 sig ORA results
res_df <- bind_rows(res_list)
str(res_df)

terms_filtered <- res_df %>%
  group_by(term_name) %>%
  mutate(term_count = n_distinct(method)) %>%
  filter(term_count > 1) %>%
  ungroup()
# ora_plots <- c(ora_sig$plot, ora_sig_up$plot, ora_sig_down$plot, ora_rank$plot, ora_meta$plot, ora_meta_up$plot, ora_meta_down$plot)
# plot
terms_plot <- terms_filtered %>%
  # filter(term_count > 2) %>%
  filter(term_size < 1750) %>%
  ggplot(aes(y = reorder(term_name, term_size), x = -log10(p_value), size = recall, colour = method)) +
  geom_point(alpha = 0.6)
ggsave(paste0(filepath, "ORA_combined_terms.png"), width = 20, height = 28, units = "cm")

ora_plots <- lapply(seq_along(ora_list), function(i) {
                   ora_list[[i]]$plot
})
s# save all plots to one PDF
pdf(paste0(filepath, "ORA_all_plots.pdf"), width = 10, height = 8)
for (plot in ora_plots) {
  print(plot)
}
dev.off()

