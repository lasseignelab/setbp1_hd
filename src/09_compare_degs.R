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
meta <- rse@metadata
sra <- args$sra

# create output dir if needed
filepath <- paste0(args$output_dir, sra, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
# save labeled and sig-filtered results
# degs_comb <- read_csv(, file = paste0(filepath, "combined_contrasts_labeled_degs.csv"))
filepath <- "./results/dea/"
cardo_degs <- read_csv(paste0(filepath, "cardo/combined_contrasts_sig_degs.csv")) %>%
  mutate(dataset = "cardo")
shaw_degs <- read_csv(paste0(filepath, "shaw/combined_contrasts_sig_degs.csv")) %>%
  mutate(dataset = "shaw")
setbp1_targets <- read_csv("./data/241203_setbp1_targets_names_only.csv", col_names = FALSE)

unique(cardo_up$contrast)
unique(shaw_up$contrast)

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
write_csv(consensus_degs, paste0(filepath, "combined/consensus_degs_unfiltered.csv"))

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

### PRIORTIZE/NARROW DEGS TO MAKE SIGNATURE ###
# opt 1) high frequency DEGs
top_freq_degs <- consensus_long %>%
  add_count(gene) %>%
  mutate(method = "majority") %>%
  filter(n > 2)

# opt 2) top DEGs by logFC
lfc_degs <- consensus_long %>%
  filter(metric == "logFC") %>%
  group_by(gene) %>%
  summarise(median_lfc = median(metric), mean_lfc = mean(metric))
str(lfc_up)

lfc_down <- consen_down_long %>%
  group_by(gene) %>%
  summarise(median_lfc = median(logFC), mean_lfc = mean(logFC))
str(lfc_down)

# by median
med_lfc_up <- lfc_up %>%
  mutate(method = "median_lfc") %>%
  slice_max(median_lfc, n = 100)

med_lfc_down <- lfc_down %>%
  mutate(method = "median_lfc") %>%
  slice_min(median_lfc, n = 100)

# by mean
mean_lfc_up <- lfc_up %>%
  mutate(method = "mean_lfc") %>%
  slice_max(mean_lfc, n = 100)
str(mean_lfc_up)

mean_lfc_down <- lfc_down %>%
  mutate(method = "mean_lfc") %>%
  slice_min(mean_lfc, n = 100)
str(mean_lfc_down)

# opt 4) top median p-val
p_up <- consen_up %>%
  pivot_longer(
    cols = starts_with("adj.P.Val_"),
    names_to = "contrast",
    names_prefix = "adj.P.Val_",
    values_to = "adj_p_val",
    values_drop_na = TRUE
  ) %>%
  group_by(gene) %>%
  summarise(median_p = median(adj_p_val), mean_p = mean(adj_p_val))
str(p_up)

p_down <- consen_down %>%
  pivot_longer(
    cols = starts_with("adj.P.Val_"),
    names_to = "contrast",
    names_prefix = "adj.P.Val_",
    values_to = "adj_p_val",
    values_drop_na = TRUE
  ) %>%
  group_by(gene) %>%
  summarise(median_p = median(adj_p_val), mean_p = mean(adj_p_val))
str(p_down)

# by median
med_p_up <- p_up %>%
  mutate(method = "median_p") %>%
  slice_min(median_p, n = 100)

med_p_down <- p_down %>%
  mutate(method = "median_p") %>%
  slice_min(median_p, n = 100)

 # by mean
mean_p_up <- p_up %>%
  mutate(method = "mean_p") %>%
  slice_min(mean_p, n = 100)

mean_p_down <- p_down %>%
  mutate(method = "mean_p") %>%
  slice_min(mean_p, n = 100)

# opt 4) robust rang aggregate -- create ranked lists per contrast
ranked_lists_up <- lapply(unique(consen_up_long$contrast), function(c) {
  consen_up_long %>%
                          filter(contrast == c) %>%
                          arrange(desc(abs(logFC))) %>%
                          pull(gene)
})

ranked_lists_down <- lapply(unique(consen_down_long$contrast), function(c) {
  consen_down_long %>%
                          filter(contrast == c) %>%
                          arrange(desc(abs(logFC))) %>%
                          pull(gene)
})

# run RRA
rra_up <- aggregateRanks(ranked_lists_up, exact = TRUE)
str(rra_up)
summary(rra_up)
rra_down <- aggregateRanks(ranked_lists_down, exact = TRUE)
str(rra_down)
summary(rra_down)

rra_results <- rra_up %>%
  mutate(direction = "up") %>%
  bind_rows(rra_down) %>%
  group_by(direction) %>%
  slice_min(order_by = Score, n = 100) %>%
  mutate(method = "RRA") %>%
  select(gene = Name, method)
 
## compare overlap between options for finding signature ##
compare_sigs_up <- majority_up %>%
  bind_rows(
    majority_down, rra_results,
    med_p_up, med_p_down, mean_p_up, mean_p_down, 
    med_lfc_up, med_lfc_down, mean_lfc_up, mean_lfc_down) %>%
  select(gene, method) %>%
    distinct() %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = method, values_from = value, values_fill = 0) %>%
  column_to_rownames(var = "gene")

str(compare_sigs_up)

png(paste0(filepath, "combined/upset_plot_compare_sigs.png"), width = 30, height = 20, units = "cm", res = 300)
upset(compare_sigs_up, sets = colnames(compare_sigs_up), order.by = "freq")
dev.off()



# save consensus signature
write_csv(consen_up_long, paste0(filepath, "combined/consensus_up_degs.csv"))
write_csv(consen_down_long, paste0(filepath, "combined/consensus_down_degs.csv"))

# compare to SETBP1 targets
# intersect(shaw_down$gene_name, setbp1_targets$X1)
intersect(consen_up_long$gene_name, setbp1_targets$X1)
intersect(consen_down_long$gene_name, setbp1_targets$X1)

consen_up_long %>%
  filter(n > 3) %>%
  ggplot(aes(x = logFC, y = reorder(gene_name, n), colour = contrast)) +
  geom_point(size = 2.5)
ggsave(paste0(filepath, "combined/consensus_up_degs.png"), height = 22, width = 18, units = "cm")

consen_down_long %>%
  filter(n > 3) %>%
  ggplot(aes(x = -logFC, y = reorder(gene_name, n), colour = contrast)) +
  geom_point(size = 2.5)
ggsave(paste0(filepath, "combined/consensus_down_degs.png"), height = 32, width = 18, units = "cm")

# consensus up DEGs upset
upset_df_up <- consen_up_long %>%
    select(gene, contrast) %>%
    distinct() %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = contrast, values_from = value, values_fill = 0) %>%
  column_to_rownames(var = "gene")
str(upset_df_up)

png(paste0(filepath, "combined/upset_plot_upregulated.png"), width = 30, height = 20, units = "cm", res = 300)
upset(upset_df_up, sets = colnames(upset_df_up), order.by = "freq")
dev.off()

# consensus down DEGs upset
upset_df_down <- consen_down_long %>%
    select(gene, contrast) %>%
    distinct() %>%
    mutate(value = 1) %>%
    pivot_wider(names_from = contrast, values_from = value, values_fill = 0) %>%
  column_to_rownames(var = "gene")
str(upset_df_down)

png(paste0(filepath, "combined/upset_plot_downregulated.png"), width = 30, height = 20, units = "cm", res = 300)
upset(upset_df_down, sets = colnames(upset_df_down), order.by = "freq")
dev.off()

#### ORA of Consensus DEGs ####
ora_up <- gost(unique(consen_up_long$gene), sources = c("GO", "KEGG", "WP"))
ora_up$result$term_name

ora_down <- gost(unique(consen_down_long$gene), sources = c("GO", "KEGG", "WP"))
ora_down$result$term_name

# plot pathways
  ora_up$result %>%
      filter(term_size < 750) %>%
      slice_max(recall, n = 40) %>%
      ggplot(aes(x = recall, y = term_name, size = intersection_size, colour = precision)) +
      geom_point()
  ggsave(paste0(filepath, "combined/up_top_ora.png"))

# plot pathways
  ora_down$result %>%
      filter(term_size < 750) %>%
      slice_max(recall, n = 40) %>%
      ggplot(aes(x = recall, y = term_name, size = intersection_size, colour = precision)) +
      geom_point()
  ggsave(paste0(filepath, "combined/down_top_ora.png"))
