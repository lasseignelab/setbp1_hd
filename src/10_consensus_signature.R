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
  metavar = "PATH",
  help = "Path to the input data directory to access validation data"
)
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
# indir <- args$indir
# degs <- args$degs
# outdir <- args$outdir
indir <- "./data/"
outdir <- "./results/consensus_signature/"
degs <- "./results/dea/"

# create output dir if needed
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# helper function to read in DEGs
read_deg <- function(dataset_name) {
  file <- file.path(degs, dataset_name, "significant_degs.csv")
  if (!file.exists(file)) stop("Missing file: ", file)
  read_csv(file) %>% mutate(dataset = dataset_name)
}

# helper function to read in MSigDB tsv files
target_ids <- c(
  "GO:0030900", # forebrain dev
  "GO:0021537", # telencephalon dev
  "REAC:R-HSA-5653656", # Reactome: PI3K/AKT Signaling
  "KEGG:04392", # Hippo signaling
  "KEGG:04064", # NF-kappa B signaling
  "KEGG:03320", # PPAR signaling
  "KEGG:04010", # MAPK signaling
  "KEGG:04151", # PI3K-Akt signaling
  "KEGG:04340", # SHH signaling
  "KEGG:04310"          # KEGG: Wnt signaling pathway
  )

pathway_df <- gconvert(query = target_ids, organism = "hsapiens", target = "HGNC", mthreshold = Inf, filter_na = TRUE)

pathways <- pathway_df %>%
  dplyr::select(input, target) %>%
  # rename(input = "STANDARD_NAME", target = "GeneSymbol") %>%
  rename("STANDARD_NAME" = input , "GeneSymbol" = target ) %>%
  filter(!is.na(GeneSymbol))
str(pathways)
# forebrain_dev <- pathways %>%
  # filter(STANDARD_NAME == "GO:0030900")

# read_msigdb <- function(file_path) {
#   lines <- readLines(file_path)
#
#   # filter for pathway name and pathway genes
#   std_name_line <- str_subset(lines, "^STANDARD_NAME")
#   gene_symbols_line <- str_subset(lines, "^GENE_SYMBOLS")
#
#   # 3. Extract the values by removing the key and the tab delimiter
#   # The '?' makes the removal non-greedy, though not strictly necessary here.
#   std_name <- str_remove(std_name_line, "^STANDARD_NAME\t")
#   gene_symbols <- str_remove(gene_symbols_line, "^GENE_SYMBOLS\t")
#
#   # 4. Create a single-row data frame
#   data.frame(
#     STANDARD_NAME = std_name,
#     GENE_SYMBOLS = gene_symbols,
#     stringsAsFactors = FALSE
#   )
# }
#

# start timer
ptm <- proc.time()

###### READ IN DATA ######
# setbp1_targets <- read_csv("./doc/241203_setbp1_targets_prioritized.csv")
cardo_degs <- read_deg("cardo")
shaw_degs <- read_deg("shaw")

# validation data
tanaka_degs <- bind_rows(read_deg("tanaka1"), read_deg("tanaka2"))
wong_fibro <- read_excel(paste0(
  indir, "wong/SourceData/SourceData6_forFIg5c_DEG.xlsx"), sheet = 3, skip = 1)

wong_neu <- read_excel(paste0(
  indir, "wong/SourceData/SourceData10_forFig7b_DEG.xlsx"), sheet = 3, skip = 1) %>%
  rename("ensemblID" = `...1`)
  # rename(`...1` = "ensemblID")

# list all pathway tsv files to read in
# file_list <- list.files(path = "./doc/MSigDB_v2024.1.Hs", pattern = "*.tsv", full.names = TRUE)
# # read in pathway standard names and gene symbols
# pathways_df <- bind_rows(lapply(file_list, read_msigdb)) %>%
#   pivot_wider(names_from = STANDARD_NAME, values_from = GENE_SYMBOLS)

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

# specifically highlighted genes in SHD
wnt_genes <- grep("WNT", consensus_degs$gene_name, value = TRUE)
core_genes <- c("MYC", "CCND1", "AXIN2", "GATA2", "DCX", "ASCL1", "NEUROD1", "HES1", "FGF8", "FOXG1", wnt_genes)
core_degs <- consensus_degs %>%
  filter(gene_name %in% core_genes) %>% arrange(gene_name)
print(n = 31, core_degs %>%
  select(gene_name, logFC, adj.P.Val, contrast))
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

top_freq_degs <- freq_degs %>%
  mutate(method = "majority") %>%
  filter(n > 4)
print(paste0(length(unique(top_freq_degs$gene)), " DEGs with freq > 3"))

# save majority/high freq signature
write_csv(
  top_freq_degs, paste0(outdir, "majorityfreq_consensus_signature.csv")
)

### opt 2) meta-analysis per gene ###
print(paste0("Finding effect size meta analysis DEGs..."))
meta_results <- consensus_degs %>%
  add_count(gene) %>%
  filter(n > 1) %>%
  group_by(gene) %>%
  nest() %>%

  mutate(meta = map(data, ~ tryCatch(
    rma(yi = .x$logFC, sei = .x$se_logFC,
      method = "REML"),
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
  left_join(consensus_degs[,1:2], by = "gene") %>%
  distinct()
summary(meta_stats)

meta_p <- meta_stats %>%
  slice_min(meta_pval, n = 300) %>%
  mutate(method = "meta_p")

# save meta analysis & signature
# meta_stats <- read_csv(paste0(outdir, "metaanalysis_degs.csv"))
write_csv(meta_stats, paste0(outdir, "metaanalysis_degs.csv"))
write_csv(meta_p, paste0(outdir, "metaanalysis_consensus_signature.csv"))

### opt 3) robust rang aggregate -- create ranked lists per contrast ###
print(paste0("Finding robust rank aggregate (RRA) DEGs..."))
ranked_lists <- lapply(unique(consensus_degs$contrast), function(c) {
  consensus_degs %>%
    filter(contrast == c) %>%
    arrange(abs(adj.P.Val)) %>%
    pull(gene)
})
# run RRA
rra <- aggregateRanks(ranked_lists, exact = TRUE)

rra_results <- rra %>%
  slice_min(order_by = Score, n = 300) %>%
  # add meta analysis for directionality
  left_join(meta_stats, by = c("Name" = "gene")) %>%
  distinct() %>%
  mutate(method = "RRA") %>%
  select(gene = Name, gene_name, Score, meta_logFC, method)

# save RRA & signature
write_csv(rra, paste0(outdir, "rankaggregate_degs.csv"))
write_csv(rra_results, paste0(outdir, "rankaggregate_consensus_signature.csv"))
# rra_results <- read_csv(paste0(outdir, "rankaggregate_consensus_signature.csv"))

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
  select(gene, gene_name, n) %>%
  left_join(entrez_degs, by = "gene") %>%
  distinct() %>%
  inner_join(meta_stats, by = "gene") %>%
  inner_join(rra, by = c("gene" = "Name")) %>%
  # inner_join(rra, by = "gene") %>%
  mutate(
    freq_weight = log2(n),
    rank_pval = Score,
    composite_score = -log10(meta_pval) + -log10(rank_pval) + freq_weight
  )

# save scores across methods for DEGs + composite scores
write_csv(
  composite_scores, paste0(outdir, "compositescore_full_table_degs.csv")
)

# top 300 by composite score
composite_signature <- composite_scores %>%
  mutate(method = "CompositeScore") %>%
  slice_max(composite_score, n = 300)
composite_signature

# save composite score signature
write_csv(
  composite_signature, paste0(
    outdir, "compositescore_consensus_signature.csv"
  )
)
###### VALIDATION ######
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

png(
  paste0(
    outdir, "upset_plot_compare_sigs.png"
  ),
  width = 30, height = 20, units = "cm", res = 300
)
upset(compare_sigs_up, sets = colnames(compare_sigs_up), order.by = "freq")
dev.off()

# compare to Tanaka et al DEG
tanaka_genes <- gorth(
  query = tanaka_degs$gene,
  source_organism = "mmusculus",
  target_organism = "hsapiens"
)

tanaka_anno <- tanaka_genes %>%
  select(gene_name = ortholog_name) %>%
  mutate(contrast = "Tanaka_DEG")

valid_degs <- tanaka_degs  %>%
  inner_join(tanaka_genes, by = c("gene" = "input_ensg")) %>%
  select(symbol = ortholog_name, ensemblID = ortholog_ensg, log2FoldChange = logFC, padj = adj.P.Val, contrast, dataset) %>%
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
  select(gene = ensemblID, gene_name = symbol, logFC = log2FoldChange, contrast, dataset) %>%
  unite("contrast_dataset", contrast, dataset) %>%
  distinct() %>%
  pivot_wider(
    names_from = contrast_dataset, values_from = logFC, values_fill = 0) %>%
  right_join(dataset_sig_degs, by = "gene") %>%
  mutate(
    across(3:length(.), ~ replace_na(.x, 0)),
    gene_name = ifelse(is.na(gene_name), "", gene_name)
  )
print(n = 50, valid_degs_wide %>% select(!gene) %>% arrange(neurons_wong))
write_csv(valid_degs_wide, paste0(outdir, "validation_datasets_consensus_DEGs.csv"))

# save
pdf(paste0(outdir, "validation_degs_full.pdf"), width = 10, height = 22)
Heatmap(as.matrix(valid_degs_wide[,3:length(valid_degs_wide)]), row_labels = valid_degs_wide$gene_name)
dev.off()

# filter for DEGs that are also in at least one validation set
valid_degs_min <- valid_degs_wide %>%
  select(!c(PATH3vsWT_shaw_logFC, IPSC_HD_shaw_logFC, HDvsWT_tanaka1, HDvsWT_tanaka2)) %>%
  filter(rowSums(select(., 3:4)) != 0) %>%
  filter(gene_name %in% pathways$GeneSymbol)

### prepare Heatmap parameters ###
colnames(valid_degs_min)
# TODO: FIX -- NOT SHOWING UP RIGHT IN HEATMAP
column_data <- data.frame(
  Study = c(rep("Wong", 2), rep("Cardo", 3), rep("Shaw", 2), "Consensus")
)

ha = HeatmapAnnotation(
  Dataset = column_data$Study,
  col = list(Study_Type = c("Cardo" = "#FF7F50", "Shaw" = "#4682B4", 
                            "Wong" = "#3CB371", "Consensus" = "#6A5ACD")),
  annotation_name_side = "left"
)

# define row annotation based on established pathways
# 1. Define the Pathway IDs and their Desired Labels
pathway_map <- list(
    "Forebrain_Dev" = "GO:0030900",       # forebrain dev
    # "Telencephalon_Dev" = "GO:0021537",   # telencephalon dev
    "PI3K_AKT" = "REAC:R-HSA-5653656", # Reactome: PI3K/AKT Signaling
    # "Hippo" = "KEGG:04392",               # Hippo signaling
    "NF-kappa_B" = "KEGG:04064",          # NF-kappa B signaling
    "PPAR" = "KEGG:03320",                # PPAR signaling
    "MAPK" = "KEGG:04010",                # MAPK signaling
    # "KEGG_PI3K_Akt" = "KEGG:04151",       # PI3K-Akt signaling
    # "SHH" = "KEGG:04340",                 # SHH signaling
    "Wnt" = "KEGG:04310"                  # KEGG: Wnt signaling pathway
)

target_pathway_ids <- unlist(pathway_map)

pathways_pivoted <- pathways %>%
    distinct() %>%
    filter(STANDARD_NAME %in% target_pathway_ids) %>% # Optional: Filter only relevant pathways
    mutate(value = 1) %>%
    pivot_wider(
        names_from = STANDARD_NAME, values_from = value, values_fill = 0
    )

# Identify missing target columns and add them filled with 0
missing_cols <- setdiff(target_pathway_ids, names(pathways_pivoted))
if(length(missing_cols) > 0) {
    for (col in missing_cols) {
        pathways_pivoted[[col]] <- 0
    }
}

# 3. Join with DEGs and use case_when for clean prioritization
pathways_wide <- valid_degs_wide %>%
    # rename(gene_name = "GeneSymbol") %>% # Assuming valid_degs_min has gene_name as the symbol
    rename("GeneSymbol" = gene_name) %>% # Assuming valid_degs_min has gene_name as the symbol
    left_join(pathways_pivoted, by = "GeneSymbol") %>%
    mutate(across(everything(), ~replace_na(., 0))) %>%
    # Assign a single, prioritized pathway using case_when
    mutate(
        Pathway = case_when(
            # Highest Priority
            `KEGG:04310` == 1 ~ "Wnt",
            `GO:0030900` == 1 ~ "Forebrain_Dev",
            # `GO:0021537` == 1 ~ "Telencephalon_Dev", # all genes in forebrain dev
            `REAC:R-HSA-5653656` == 1 ~ "PI3K_AKT",
            # `REAC:R-HSA-5653656` == 1 ~ "Reactome_PI3K_AKT",
            # `KEGG:04392` == 1 ~ "Hippo", # no genes annotated
            `KEGG:04064` == 1 ~ "NF-kappa_B",
            `KEGG:04010` == 1 ~ "MAPK",
            # `KEGG:04340` == 1 ~ "SHH", # no genes annotated
            `KEGG:03320` == 1 ~ "PPAR",
            # Lowest Priority / Default: if none of the above are 1
            TRUE ~ "Other"
        )
    ) %>%
    # Prepare for ComplexHeatmap annotation
    column_to_rownames(var = "gene") %>%
    # Select only the columns needed for the annotation dataframe (GeneSymbol, Pathway, and the binary indicators)
    select(GeneSymbol, Pathway, all_of(target_pathway_ids))
str(pathways_wide)


# define the correct binary column names
binary_col_names <- unlist(names(pathway_map)) # Use the IDs as column names

# set colors
unique_pathways <- unique(pathways_wide$Pathway)
pathway_colors <- setNames(
    colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(length(unique_pathways)),
    unique_pathways
)

binary_colors <- c("0" = "white", "1" = "gray70")
binary_color_list <- rep(list(binary_colors), length(binary_col_names))
names(binary_color_list) <- binary_col_names

# creat gene row annotation
gene_ha = rowAnnotation(
    df = pathways_wide %>% filter(row.names(.) %in% valid_degs_min$gene),#%>% select(all_of(binary_col_names)),
    col = binary_color_list, #list(
    width = unit(7, "cm"),
    annotation_name_side = "top"#,
  )

# color scale for logFC
max_logFC <- max(abs(meta_stats$meta_logFC))
col_fun = circlize::colorRamp2(c(-max_logFC, 0, max_logFC), c("blue", "white", "red"))

### plot heatmap ###
heatmap_obj = Heatmap(
    as.matrix(valid_degs_min[,3:length(valid_degs_min)]), 
  name = "Log2(Fold Change)",
  col = col_fun,
  row_labels = valid_degs_min$gene_name,
  row_names_gp = gpar(fontsize = 6), 
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  top_annotation = ha,
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
    column_to_rownames("gene")# %>%
    # select(!gene_name) %>%
    # as.matrix()
  print(head(plot_data))

  annotation_df <- gene_data %>%
  filter(row.names(.) %in% genes) %>%
    select(all_of(names(target_pathway_ids)))

gene_ha = rowAnnotation(
    df = annotation_df %>% select(all_of(binary_col_names)),
    col = binary_color_list, #list(
    width = unit(7, "cm"),
    annotation_name_side = "top"#,
  )


    # create plot
  hm <- Heatmap(
    as.matrix(plot_data[,2:length(plot_data)]), 
  name = "Log2(Fold Change)",
  col = col_fun,
  row_labels = plot_data$gene_name,
  row_names_gp = gpar(fontsize = 6), 
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  # top_annotation = ha,
  right_annotation = gene_ha,
  border = TRUE,
  column_title = paste0(title, "-Based SETBP1-HD Consensus Signature")
    # as.matrix(plot_data[,3:length(plot_data)]), 
    # row_labels = plot_data$gene_name,
    # column_title = title
  )
     list(plot = hm)
}

# plot heatmap for each method
sig_list <- list(
  Composite = c(composite_signature)$gene,
  RRA = c(rra_results)$gene,
  Meta = c(meta_p)$gene,
  Freq = c(top_freq_degs)$gene
)
# sig_list <- list(composite_signature, rra_results, meta_p, top_freq_degs, core_sig)
hm_list <- lapply(seq_along(sig_list), function(sig) {
  heatmap_by_method(
    full_deg_df = valid_degs_wide,
    genes = sig_list[[sig]], 
    gene_data = pathways_wide,
    title = names(sig_list)[[sig]])
})

# save all plots to one PDF
pdf(paste0(outdir, "heatmap_all_sigs_validation.pdf"), width = 10, height = 18)
for (plot in hm_list) {
  print(plot)
}
dev.off()

### COMPARE DEGs ACROSS CONTRASTS + DATA SETS ###
compare_contrasts <- consensus_long %>%
  select(gene_name, contrast) %>%
  # bind_rows(tanaka_anno) %>%
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
    outdir, "upset_plot_compare_shaw_cardo.png"
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
    outdir, "upset_plot_compare_contrasts.png"
  ),
  width = 30, height = 20, units = "cm", res = 300
)
upset(compare_contrasts, sets = colnames(compare_contrasts), order.by = "freq")
dev.off()

### RRA and Meta Analysis Volcano Plot ###
# set logFC and p-val thresholds
logFC_threshold <- 0.1
pval_threshold <- -log10(0.05)

plot_data <- rra %>%
  inner_join(meta_stats, by = c("Name" = "gene")) %>%
  mutate(
    # calculate -log10(p-value RRA)
    RRA_pval = -log10(Score),
    Direction = case_when(
      meta_logFC > logFC_threshold & RRA_pval > pval_threshold ~ "Upregulated",
      meta_logFC < -logFC_threshold & RRA_pval > pval_threshold ~ "Downregulated",
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
  # filter(gene_name %in% core_genes)

volcano_plot <- ggplot(plot_data, aes(x = meta_logFC, y = RRA_pval)) +
  geom_point(aes(color = Direction), alpha = 0.6, size = 2) +
  scale_color_manual(values = c("Upregulated" = "firebrick3",
                                "Downregulated" = "dodgerblue3",
                                "Non-Significant/Inconsistent" = "gray50")) +
  # horizontal threshold line for RRA significance
  geom_hline(yintercept = pval_threshold, linetype = "dashed", color = "black") +
  # vertical threshold lines for logFC magnitude
  geom_vline(xintercept = c(-logFC_threshold, logFC_threshold), linetype = "dashed", color = "black") +
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
    segment.color = 'grey50'
  )
ggsave(paste0(outdir, "RRA_meta_volcano.png"), dpi = 300, height = 30, width = 26, units = "cm")

###### ORA OF CONSENSUS SIGNATURES ######
# function for ORA plotting
run_ora_plot <- function(genes, title) {
  ora_res <- tryCatch(
    {
      gost(genes, sources = c("GO", "KEGG", "WP"), evcodes = TRUE)$result
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
        aes(
          x = recall, y = term_name,
          size = intersection_size, colour = precision
        )
      ) +
      geom_point() +
      labs(title = title, x = "Recall", y = NULL) +
      theme_minimal()

    list(result = ora_res, plot = plot)
  }
}

# run ORA for each signature
list_sig_genes <- list(
  Composite = c(composite_signature$gene),
  RRA = c(rra_results$gene),
  Meta = c(meta_p$gene),
  Freq = c(top_freq_degs$gene)
)
ora_all <- gost(list_sig_genes, sources = c("GO", "KEGG", "WP"), evcodes = TRUE)
results_df <- ora_all$result

# Tidy the data for plotting: Select only significant terms and calculate -log10(p_value)
plot_df <- results_df %>%
  filter(significant == TRUE) %>%
  # filter(query != "Freq") %>%
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
  facet_grid(source ~ ., scales = "free_y", space = "free_y") # Separate by database
ggsave(paste0(outdir, "ORA_comparison_dotplot.png"), dpi = 300, height = 30, width = 26, units = "cm")

plot_df %>%
  group_by(query, source) %>%
  summarise(
    n_significant_terms = n(),
    .groups = 'drop'
  ) %>%
  ggplot(aes(x = query, y = n_significant_terms, fill = source)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  labs(
    x = "Query Condition",
    y = "Count of Significant Terms",
    fill = "Data Source"
  ) +
  theme_classic()
ggsave(paste0(outdir, "ORA_term_counts.png"), dpi = 300, height = 22, width = 18, units = "cm")

### Signature Plotting ###
# select wnt genes
wnt_terms <- grep("Wnt", results_df$term_name, value = TRUE)
wnt_terms
wnt_genes <- results_df %>%
  filter(term_name %in% wnt_terms) %>%
  select(term_name, intersection) %>%
  distinct() %>%
  separate_longer_delim(intersection, delim = ",")
wnt_genes

composite_signature_wnt <- composite_scores %>%
  filter(
    gene %in% unique(wnt_genes$intersection) | gene_name.x %in% core_genes
  ) %>%
  select(gene, gene_name = gene_name.x, meta_logFC, meta_pval) %>%
  distinct() %>%
  column_to_rownames(var = "gene")
length(unique(composite_signature_wnt$gene_name))
intersect(composite_signature$gene_name.x, composite_signature_wnt$gene_name)
intersect(rownames(composite_signature_wnt), rra_results$gene)


# save
pdf(paste0(outdir, "composite_signature_wnt_genes.pdf"), width = 10, height = 18)
Heatmap(as.matrix(composite_signature_wnt[,2:3]), row_labels = composite_signature_wnt$gene_name)
dev.off()

# plot DEGs of signature
wnt_degs <- composite_signature_wnt %>%
  inner_join(consensus_wide[,-1], by = "gene_name") %>%
  select(gene_name, contains("_logFC")) %>%
  select(!c(contains("_se_logFC")))
head(wnt_degs)

wnt_degs[is.na(wnt_degs)] <- 0
# save
pdf(paste0(outdir, "composite_wnt_degs.pdf"), width = 10, height = 18)
Heatmap(as.matrix(wnt_degs[,2:length(wnt_degs)]), row_labels = wnt_degs$gene_name)
dev.off()
write_csv(wnt_degs, paste0(outdir, "composite_wnt_degs.csv"))


#### END ####

# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
