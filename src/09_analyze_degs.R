suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(ggplot2)
  library(UpSetR)
  library(circlize)
  library(gprofiler2)
})

# create parser object
parser <- ArgumentParser()

parser$add_argument("-d", "--degs",
  type = "character",
  metavar = "PATH",
  help = "Significant DEGS as a CSV file (rds format) from limma-voom analysis with contrasts combined as 'contrast' column" # nolint
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

# create output dir if needed
filepath <- paste0(args$output_dir, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

# filepath <- "./results/archive/"
# sig_degs_comb <- read_csv(file = "./results/dea/cardo/combined_contrasts_sig_degs.csv")

#### Compare Up/Down DEGs Across Conditions ####
# split significant DEGs by direction
up_degs_comb <- sig_degs_comb %>% filter(logFC > 0)
down_degs_comb <- sig_degs_comb %>% filter(logFC < 0)

# gene-contrast matrix for upregulated DEGs
up_gene_contrast_matrix <- up_degs_comb %>%
  select(gene, contrast) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = contrast, values_from = value, values_fill = 0)

# gene-contrast matrix for downregulated DEGs
down_gene_contrast_matrix <- down_degs_comb %>%
  select(gene, contrast) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = contrast, values_from = value, values_fill = 0)

# prep dataframes for UpSetR
up_gene_contrast_matrix <- as.data.frame(up_gene_contrast_matrix)
rownames(up_gene_contrast_matrix) <- up_gene_contrast_matrix$gene
up_gene_contrast_matrix$gene <- NULL

down_gene_contrast_matrix <- as.data.frame(down_gene_contrast_matrix)
rownames(down_gene_contrast_matrix) <- down_gene_contrast_matrix$gene
down_gene_contrast_matrix$gene <- NULL

# up DEGs upset plot
png(paste0(filepath, "upset_plot_upregulated.png"), width = 30, height = 20, units = "cm", res = 300)
upset(up_gene_contrast_matrix, sets = colnames(up_gene_contrast_matrix), order.by = "freq")
dev.off()

# down DEGs upset plot
png(paste0(filepath, "upset_plot_downregulated.png"), width = 30, height = 20, units = "cm", res = 300)
upset(down_gene_contrast_matrix, sets = colnames(down_gene_contrast_matrix), order.by = "freq")
dev.off()

write_csv(up_degs_comb, paste0(filepath, "upregulated_degs_combined.csv"))
write_csv(down_degs_comb, paste0(filepath, "downregulated_degs_combined.csv"))

# plot shared DEGs
sig_degs_freq <- sig_degs_comb %>%
  add_count(gene)
print(paste0("Frequency of shared DEGs across contrasts (1=only in one contrast)"))
table(sig_degs_freq$n)

sig_degs_freq %>%
  filter(contrast != "PATH_vs_PATH") %>%
  filter(n > 2) %>%
  ggplot(aes(x = logFC, y = gene_name, colour = contrast)) +
  geom_point()
ggsave(paste0(filepath, "sig_degs_across_conditions_dotplot.png"))

#### ORA ####
consensus_up <- sig_degs_freq %>%
  # rm comparisons that are not SETBP1 HD to WT comparisons
  filter(contrast != "PATH_vs_PATH") %>%
  filter(n > 1 & logFC > 0) 
consensus_ora_up <- gost(consensus_up$gene, sources = "GO")
write_csv(consensus_ora_up$result, paste0(filepath, "consensus_up_degs_ora.csv"))

consensus_down <- sig_degs_freq %>%
  # rm comparisons that are not SETBP1 HD to WT comparisons
  filter(contrast != "PATH_vs_PATH") %>%
  filter(n > 1 & logFC < 0) 
consensus_ora_down <- gost(consensus_down$gene, sources = "GO")
write_csv(consensus_ora_down$result, paste0(filepath, "consensus_down_degs_ora.csv"))

conditions <- unique(sig_degs_comb$contrast)

 for(i in seq_along(colnames(contr))) {
  print(paste0(colnames(contr)[i]))
  contrast <- colnames(contr)[i]

  up_degs <- sig_degs_comb %>%
    filter(contrast == paste0(contrast)) %>%
    filter(logFC > 0)

  ora_up <- gost(up_degs$gene, sources = "GO")
  write_csv(ora_up$result, paste0(filepath, contrast, "_up_degs_ora.csv"))

  down_degs <- sig_degs_comb %>%
    filter(contrast == paste0(contrast)) %>%
    filter(logFC < 0)

  ora_down <- gost(down_degs$gene, sources = "GO")
  write_csv(ora_down$result, paste0(filepath, contrast, "_down_degs_ora.csv"))

  # Plot Pathways
  ora_up$result %>%
    filter(term_size < 750) %>%
    slice_max(recall, n = 30) %>%
    ggplot(aes(x = recall, y = term_name, size = intersection_size, colour = precision)) +
    geom_point()
  ggsave(paste0(filepath, contrast, "_up_top_ora.png"))

  ora_down$result %>%
    filter(term_size < 750) %>%
    slice_max(recall, n = 30) %>%
    ggplot(aes(x = recall, y = term_name, size = intersection_size, colour = precision)) +
    geom_point()
  ggsave(paste0(filepath, contrast, "_down_top_ora.png"))


}

# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
