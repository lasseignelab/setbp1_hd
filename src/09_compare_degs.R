suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(UpSetR)
  library(ComplexHeatmap)
  library(statmod)
  library(gprofiler2)
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
cardo_up <- read_csv(file = paste0(filepath, "cardo/upregulated_degs_combined.csv"))
cardo_down <- read_csv(file = paste0(filepath, "cardo/downregulated_degs_combined.csv"))
shaw_up <- read_csv(file = paste0(filepath, "shaw/upregulated_degs_combined.csv"))
shaw_down <- read_csv(file = paste0(filepath, "shaw/downregulated_degs_combined.csv"))

setbp1_targets <- read_csv("./data/241203_setbp1_targets_names_only.csv", col_names = FALSE)

unique(cardo_up$contrast)
unique(shaw_up$contrast)

# function to find consensus between cardo/shaw DEGs
find_consensus <- function(df1, df2, contrasts) {
  df2_wide <- df2 %>%
    filter(contrast %in% paste0(contrasts)) %>%
    select(gene, contrast, logFC, adj.P.Val) %>%
    pivot_wider(
      names_from = contrast,
      values_from = c(logFC, adj.P.Val))

  consensus_df <- df1 %>%
    filter(contrast %in% paste0(contrasts)) %>%
    select(!c(AveExpr, t, P.Value, B)) %>%
    pivot_wider(
      names_from = contrast,
      values_from = c(logFC, adj.P.Val)) %>%
    inner_join(df2_wide, by = "gene")

  return(consensus_df)
}

# contrasts to keep
contrasts_keep <- c("KOvsWT_15", "KOvsWT_21", "KOvsWT_34", "VUSvsWT", "PATH2vsWT", "PATH3vsWT")
# consensus up DEGs
consen_up <- find_consensus(cardo_up, shaw_up, contrasts_keep)
# consen_up <- cardo_up %>%
#   filter(contrast %in% c("KOvsWT_15", "KOvsWT_21", "KOvsWT_34")) %>%
#   pivot_wider(
#     names_from = contrast,
#     values_from = c(logFC, adj.P.Val))
#   # inner_join(shaw_up, by = "gene") %>%
  # filter(contrast %in% c("IPSC_HD", "NPC_HD"))
str(consen_up)

# consensus down DEGs
consen_down <- find_consensus(cardo_down, shaw_down, contrasts_keep)
str(consen_down)

# compare to SETBP1 targets
# intersect(shaw_down$gene_name, setbp1_targets$X1)
intersect(consen_up$gene_name, setbp1_targets$X1)
intersect(consen_down$gene_name, setbp1_targets$X1)

# save consensus signature
write_csv(consen_up, paste0(filepath, "combined/consensus_up_degs_wide_unfiltered.csv"))
write_csv(consen_down, paste0(filepath, "combined/consensus_down_degs_wide_unfiltered.csv"))

# long format for plotting
consen_up_long <- consen_up %>%
  pivot_longer(
    cols = starts_with("logFC"),
    names_to = "contrast",
    names_prefix = "logFC_",
    values_to = "logFC",
    values_drop_na = TRUE
  ) %>%
  add_count(gene)
length(unique(consen_up_long$gene))
str(consen_up_long)
table(consen_up_long$contrast)

consen_down_long <- consen_down %>%
  pivot_longer(
    cols = starts_with("logFC"),
    names_to = "contrast",
    names_prefix = "logFC_",
    values_to = "logFC",
    values_drop_na = TRUE
  ) %>%
  add_count(gene)
length(unique(consen_down_long$gene))
str(consen_down_long)
table(consen_down_long$contrast)

### PRIORTIZE/NARROW DEGS TO MAKE SIGNATURE ###
# opt 1) majority vote DEGs
majority_up <- consen_up_long %>%
  filter(n > 2)

majority_down <- consen_down_long %>%
  filter(n > 2)
str(majority_down)

# opt 2) top median logFC
med_lfc_up <- consen_up_long %>%
  group_by(gene) %>%
  summarise(median_lfc = median(logFC)) %>%
  slice_max(median_lfc, n = 100)
str(med_lfc_up)

med_lfc_down <- consen_down_long %>%
  group_by(gene) %>%
  summarise(median_lfc = median(logFC)) %>%
  slice_min(median_lfc, n = 100)
str(med_lfc_down)

# opt 3) top mean logFC
mean_lfc_up <- consen_up_long %>%
  group_by(gene) %>%
  summarise(mean_lfc = mean(logFC)) %>%
  slice_max(mean_lfc, n = 100)
str(mean_lfc_up)

mean_lfc_down <- consen_down_long %>%
  group_by(gene) %>%
  summarise(mean_lfc = mean(logFC)) %>%
  slice_min(mean_lfc, n = 100)
str(mean_lfc_down)

# opt 4) top median p-val
med_p_up <- consen_up %>%
  pivot_longer(
    cols = starts_with("adj.P.Val_"),
    names_to = "contrast",
    names_prefix = "adj.P.Val_",
    values_to = "adj_p_val",
    values_drop_na = TRUE
  ) %>%
  group_by(gene) %>%
  summarise(median_p = median(adj_p_val), mean_p = mean(adj_p_val)) %>%
  slice_min(median_p, n = 100)
str(med_p_up)

med_p_down <- consen_down %>%
  pivot_longer(
    cols = starts_with("adj.P.Val_"),
    names_to = "contrast",
    names_prefix = "adj.P.Val_",
    values_to = "adj_p_val",
    values_drop_na = TRUE
  ) %>%
  group_by(gene) %>%
  summarise(median_p = median(adj_p_val), mean_p = mean(adj_p_val)) %>%
  slice_min(median_p, n = 100)
str(med_p_down)

## overlap between options for finding signature ##
length(intersect(majority_up$gene, med_lfc_up$gene)) #28
length(intersect(majority_down$gene, med_lfc_down$gene)) # 18

length(intersect(mean_lfc_up$gene, med_lfc_up$gene)) #92
length(intersect(mean_lfc_down$gene, med_lfc_down$gene)) #92

length(intersect(majority_up$gene, mean_lfc_up$gene)) # 28
length(intersect(majority_down$gene, mean_lfc_down$gene)) #21

length(intersect(majority_up$gene, med_p_up$gene)) # 91
length(intersect(majority_down$gene, med_p_down$gene)) # 100


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
