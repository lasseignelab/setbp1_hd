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
read_deg <- function(dataset_name) {
  file <- file.path(degs, dataset_name, "significant_degs.csv")
  if (!file.exists(file)) stop("Missing file: ", file)
  read_csv(file) %>% mutate(dataset = dataset_name)
}

# start timer
ptm <- proc.time()

###### READ IN DATA ######
cardo_degs <- read_deg("cardo")
shaw_degs <- read_deg("shaw")

# validation data
tanaka_degs <- bind_rows(read_deg("tanaka1"), read_deg("tanaka2"))
wong_fibro <- read_excel(paste0(
  indir, "wong/SourceData/SourceData6_forFIg5c_DEG.xlsx"
), sheet = 3, skip = 1)

wong_neu <- read_excel(paste0(
  indir, "wong/SourceData/SourceData10_forFig7b_DEG.xlsx"
), sheet = 3, skip = 1) %>%
  rename("ensemblID" = `...1`)

###### FIND CONSENSUS DEGS ######
unique(cardo_degs$contrast)
unique(shaw_degs$contrast)

# exclude redundant/biased contrasts
contrasts_skip <- c(
  "KOvsWT_NPC", "KOvsWT_all", "NPC_HD", "PATH_vs_PATH"
)

# combine human data sets and filter
combined_degs <- bind_rows(cardo_degs, shaw_degs)

consensus_degs <- combined_degs %>%
  filter(!contrast %in% contrasts_skip) %>%
  add_count(gene) %>%
  filter(n != 1) %>%
  select(!n)
consensus_degs
table(consensus_degs$contrast)

# save consensus DEGs
write_csv(consensus_degs, paste0(outdir, "consensus_degs_unfiltered.csv"))

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
  add_count(gene) %>%
  filter(n > 1) %>%
  group_by(gene) %>%
  nest() %>%
  mutate(meta = map(data, ~ tryCatch(
    rma(
      yi = .x$logFC, sei = .x$se_logFC,
      method = "REML"
    ),
    error = function(e) {
      print(paste("Error for gene:", unique(.x$gene_name), e$message))
      NULL
    }
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
  left_join(consensus_degs[, 1:2], by = "gene") %>%
  distinct()
summary(meta_stats)

# save meta analysis
write_csv(meta_stats, paste0(outdir, "metaanalysis_degs.csv"))

### opt 3) robust rang aggregate -- create ranked lists per contrast ###
print(paste0("Finding robust rank aggregate (RRA) DEGs..."))
ranked_lists <- lapply(unique(consensus_degs$contrast), function(c) {
  consensus_degs %>%
    filter(contrast == c) %>%
    arrange(abs(adj.P.Val)) %>%
    pull(gene)
})
# run and save RRA
rra <- aggregateRanks(ranked_lists, exact = TRUE)

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

png(
  paste0(
    outdir, "upset_plot_compare_sigs.png"
  ),
  width = 30, height = 20, units = "cm", res = 300
)
upset(compare_sigs, sets = colnames(compare_sigs), order.by = "freq")
dev.off()

# compare to Tanaka et al DEG
tanaka_genes <- gorth(
  query = tanaka_degs$gene,
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)

valid_degs <- tanaka_degs %>%
  inner_join(tanaka_genes, by = c("gene" = "input_ensg")) %>%
  select(
    symbol = ortholog_name,
    ensemblID = ortholog_ensg,
    log2FoldChange = logFC,
    padj = adj.P.Val, contrast, dataset
  ) %>%
  distinct() %>%
  bind_rows(wong_fibro) %>%
  mutate(
    contrast = ifelse(is.na(contrast), "fibroblasts", contrast),
  ) %>%
  bind_rows(wong_neu) %>%
  mutate(
    contrast = ifelse(is.na(contrast), "neurons", contrast),
    dataset = ifelse(is.na(dataset), "wong", dataset)
  ) %>%
  filter(!is.na(ensemblID))
write_csv(valid_degs, paste0(outdir, "validation_datasets_DEGs.csv"))

# plot DEGs across data sets for signature genes
dataset_sig_degs <- consensus_wide %>%
  select(gene, contains("_logFC")) %>%
  inner_join(composite_scores, by = "gene") %>%
  select(gene, contains("_logFC")) %>%
  select(!contains("se_logFC")) %>%
  distinct()
dataset_sig_degs

valid_degs_wide <- valid_degs %>%
  select(
    gene = ensemblID,
    gene_name = symbol,
    logFC = log2FoldChange,
    contrast, dataset
  ) %>%
  unite("contrast_dataset", contrast, dataset) %>%
  distinct() %>%
  pivot_wider(
    names_from = contrast_dataset, values_from = logFC, values_fill = 0
  ) %>%
  right_join(dataset_sig_degs, by = "gene") %>%
  mutate(
    across(3:length(.), ~ replace_na(.x, 0)),
    gene_name = ifelse(is.na(gene_name), "", gene_name)
  )
print(n = 20, valid_degs_wide %>% select(!gene) %>% arrange(neurons_wong))
write_csv(
  valid_degs_wide, paste0(outdir, "validation_datasets_consensus_DEGs.csv")
)

# save
pdf(paste0(outdir, "validation_degs_full.pdf"), width = 10, height = 22)
Heatmap(
  as.matrix(
    valid_degs_wide[, 3:length(valid_degs_wide)]
  ),
  row_labels = valid_degs_wide$gene_name
)
dev.off()

### gene functional pathway annotations ###
target_ids <- c(
  "GO:0030900", # forebrain dev
  "REAC:R-HSA-5653656", # Reactome: PI3K/AKT Signaling
  "KEGG:04064", # NF-kappa B signaling
  "KEGG:03320", # PPAR signaling
  "KEGG:04010", # MAPK signaling
  "KEGG:04310" # KEGG: Wnt signaling pathway
)

pathway_df <- gconvert(
  query = target_ids,
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
  select(
    !c(
      PATH3vsWT_shaw_logFC,
      IPSC_HD_shaw_logFC,
      HDvsWT_tanaka1,
      HDvsWT_tanaka2
    )
  ) %>%
  filter(rowSums(select(., 3:4)) != 0) %>%
  filter(gene_name %in% pathways$GeneSymbol)

### prepare Heatmap parameters ###

# define row annotation based on established pathways
pathway_map <- list(
  "Forebrain_Dev" = "GO:0030900", # forebrain dev
  "PI3K_AKT" = "REAC:R-HSA-5653656", # Reactome: PI3K/AKT Signaling
  "NF-kappa_B" = "KEGG:04064", # NF-kappa B signaling
  "PPAR" = "KEGG:03320", # PPAR signaling
  "MAPK" = "KEGG:04010", # MAPK signaling
  "Wnt" = "KEGG:04310" # KEGG: Wnt signaling pathway
)

target_pathway_ids <- unlist(pathway_map)

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
  mutate(across(everything(), ~ replace_na(., 0))) %>%
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

# define the correct binary column names
binary_col_names <- unlist(names(pathway_map)) # Use the IDs as column names

# set colors
unique_pathways <- unique(pathways_wide$Pathway)
pathway_colors <- setNames(
  colorRampPalette(
    RColorBrewer::brewer.pal(12, "Set3")
  )(length(unique_pathways)),
  unique_pathways
)

binary_colors <- c("0" = "white", "1" = "gray70")
binary_color_list <- rep(list(binary_colors), length(binary_col_names))
names(binary_color_list) <- binary_col_names

# creat gene row annotation
gene_ha <- rowAnnotation(
  df = pathways_wide %>% filter(row.names(.) %in% valid_degs_min$gene),
  col = binary_color_list,
  width = unit(7, "cm"),
  annotation_name_side = "top"
)

# color scale for logFC
max_lfc <- max(abs(meta_stats$meta_logFC))
col_fun <- circlize::colorRamp2(
  c(-max_lfc, 0, max_lfc), c("blue", "white", "red")
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
pdf(paste0(outdir, "validation_degs_min.pdf"), width = 10, height = 22)
draw(heatmap_obj, column_title_gp = gpar(fontsize = 14, fontface = "bold"))
dev.off()

heatmap_by_method <- function(full_deg_df, genes, gene_data, title) {
  # filter and select sig degs
  plot_data <- full_deg_df %>%
    filter(gene %in% genes) %>%
    column_to_rownames("gene")
  print(head(plot_data))

  annotation_df <- gene_data %>%
    filter(row.names(.) %in% genes) %>%
    select(all_of(names(target_pathway_ids)))

  gene_ha <- rowAnnotation(
    df = annotation_df %>% select(all_of(binary_col_names)),
    col = binary_color_list,
    width = unit(7, "cm"),
    annotation_name_side = "top"
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
pdf(
  paste0(outdir, "heatmap_all_sigs_validation.pdf"),
  width = 10, height = 18
)
for (plot in hm_list) {
  print(plot)
}
dev.off()

### COMPARE DEGs ACROSS CONTRASTS + DATA SETS ###
compare_contrasts <- consensus_long %>%
  select(gene_name, contrast) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(
    names_from = contrast, values_from = value, values_fill = 0
  ) %>%
  distinct() %>%
  column_to_rownames(var = "gene_name")

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
    outdir, "upset_plot_compare_datasets.png"
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

### RRA and Meta Analysis Volcano Plot ###
lfc_th <- 0.1 # logFC threshold
pval_th <- -log10(0.05) # p-value threshold

plot_data <- rra %>%
  inner_join(meta_stats, by = c("Name" = "gene")) %>%
  mutate(
    RRA_pval = -log10(Score),
    Direction = case_when(
      meta_logFC > lfc_th & RRA_pval > pval_th ~ "Upregulated",
      meta_logFC < -lfc_th & RRA_pval > pval_th ~ "Downregulated",
      TRUE ~ "Non-Significant/Inconsistent"
    )
  )
table(plot_data$Direction)
# save as meta-analysis labeled RRA signature
write_csv(
  plot_data, paste0(outdir, "RRA_meta_consensus_signature.csv")
)
core_genes_label <- plot_data %>%
  slice_max(RRA_pval, n = 10)

volcano_plot <- ggplot(plot_data, aes(x = meta_logFC, y = RRA_pval)) +
  geom_point(aes(color = Direction), alpha = 0.6, size = 2) +
  scale_color_manual(values = c(
    "Upregulated" = "firebrick3",
    "Downregulated" = "dodgerblue3",
    "Non-Significant/Inconsistent" = "gray50"
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
    title = "Consensus SETBP1-HD Signature: MA logFC vs. RRA Consistency",
    x = "MA Aggregate logFC (Magnitude)",
    y = "-log10(RobustRankAggregate P-value)",
    color = "Consensus Status"
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# add text labels for top genes
labeled_volcano_plot <- volcano_plot +
  ggrepel::geom_text_repel(
    data = core_genes_label,
    aes(label = gene_name),
    size = 3.5,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = "grey50"
  )
ggsave(paste0(outdir, "RRAp_metalogFC_volcano.png"),
  dpi = 300, height = 30, width = 26, units = "cm"
)

###### ORA OF CONSENSUS SIGNATURES ######
# run ORA for each signature
ora_all <- gost(sig_list, sources = c("GO", "KEGG", "WP"), evcodes = TRUE)
results_df <- ora_all$result

# prep data for plotting
plot_df <- results_df %>%
  filter(significant == TRUE) %>%
  filter(term_size < 400) %>%
  mutate(
    neg_log10_p = -log10(p_value),
    # Simplify term names for better readability
    term_name = gsub("Homo sapiens ", "", term_name)
  )

# Create the comparative dot plot
ggplot(plot_df, aes(x = query, y = reorder(term_name, neg_log10_p))) +
  geom_point(aes(size = intersection_size, color = neg_log10_p)) +
  scale_color_viridis_c(name = "-log10(p-value)") +
  labs(
    x = "Query Condition",
    y = "Enriched Term",
    size = "Gene Overlap Size"
  ) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  facet_grid(source ~ ., scales = "free_y", space = "free_y")
ggsave(paste0(outdir, "ORA_comparison_dotplot.png"),
  dpi = 300, height = 30, width = 26, units = "cm"
)

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
  theme_classic()
ggsave(paste0(outdir, "ORA_term_counts_barplot.png"),
  dpi = 300, height = 22, width = 18, units = "cm"
)

#### END ####

# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
