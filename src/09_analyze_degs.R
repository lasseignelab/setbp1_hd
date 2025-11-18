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

# start timer
ptm <- proc.time()

#### Set variables for dataset ####
print(paste0("Reading in DEGs from: ", args$degs))
sig_degs_comb <- read_csv(args$degs)
filepath <- args$output_dir

# create output dir if needed
filepath <- paste0(args$output_dir, "/")
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

# filepath <- "./results/archive/"
# sig_degs_comb <- read_csv(file = "./results/dea/shaw/combined_contrasts_sig_degs.csv")

#### Compare Up/Down DEGs Across Conditions ####
# gene-contrast matrix for upregulated DEGs
up_gene_contrast_matrix <- sig_degs_comb %>%
    filter(logFC < 0) %>%
  select(gene, contrast) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = contrast, values_from = value, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene")

# gene-contrast matrix for downregulated DEGs
down_gene_contrast_matrix <- sig_degs_comb %>%
    filter(logFC > 0) %>%
  select(gene, contrast) %>%
  distinct() %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = contrast, values_from = value, values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames(var = "gene")

# DEGs upset plots by direction (up/down)
upset_up <- upset(
  up_gene_contrast_matrix, 
  sets = colnames(up_gene_contrast_matrix), 
  order.by = "degree",
  mainbar.y.label = "Upregulated Intersection Size",
sets.bar.color = "steelblue")

upset_down <- upset(
  down_gene_contrast_matrix, 
  sets = colnames(down_gene_contrast_matrix), 
  order.by = "degree",
  mainbar.y.label = "Downregulated Intersection Size",
sets.bar.color = "firebrick")

pdf(paste0(filepath, "upset_plots_across_contrasts.pdf"))
upset_up
upset_down
dev.off()
print(
  paste0(
    "Upset plots saved to: ", filepath, "upset_plots_across_contrasts.pdf"))


#### ORA ####
# function for ORA plotting
run_ora_plot <- function(genes, title) {
  ora_res <- tryCatch({
    gost(genes, sources = "GO")$result
  }, error = function(e) {
    message(paste("ORA failed for", title, ":", e$message))
    return(NULL)
  })
  # skip plotting for no ORA results
  if (is.null(ora_res) || nrow(ora_res) == 0) {
    message(paste("No ORA results for", title, "-- skipping plot"))
    return(NULL)
  }
  # filter and select top terms
  plot_data <- ora_res %>%
    filter(term_size < 750) %>%
    slice_max(recall, n = 30)
  
  # create plot
  plot <- ggplot(plot_data, aes(x = recall, y = term_name, size = intersection_size, colour = precision)) +
    geom_point() +
    labs(title = title, x = "Recall", y = NULL) +
    theme_minimal()
  
  return(list(result = ora_res, plot = plot))
}

# consensus DEGs ORA for each direction (up/down)
sig_degs_freq <- sig_degs_comb %>%
  add_count(gene) %>%
  filter(contrast != "PATH_vs_PATH", n > 1)

consensus_up <- sig_degs_freq %>% filter(logFC > 0)
consensus_down <- sig_degs_freq %>% filter(logFC < 0)

consensus_up_ora <- run_ora_plot(consensus_up$gene, "Consensus Upregulated")
consensus_down_ora <- run_ora_plot(consensus_down$gene, "Consensus Downregulated")

write_csv(consensus_up_ora$result, paste0(filepath, "consensus_up_degs_ora.csv"))
write_csv(consensus_down_ora$result, paste0(filepath, "consensus_down_degs_ora.csv"))
print(
  paste0(
    "Consensus ORA results saved as: ", "consensus_[up/down]_degs_ora.csv"))

contrasts <- unique(sig_degs_comb$contrast)
plots <- list(consensus_up_ora$plot, consensus_down_ora$plot)

# run ORA for DEGs from each contrast
 for(contr in seq_along(contrasts)) {
  print(paste0("Running ORA for ", contrasts[contr]))
  # split by upregulated genes, run ORA
  up_degs <- sig_degs_comb %>%
    filter(contrast == contrasts[contr], logFC > 0)
 
  ora_up <- run_ora_plot(up_degs$gene, paste0(contrasts[contr], " Upregulated"))
  if (!is.null(ora_up)) {
    write_csv(ora_up$result, paste0(filepath, contrasts[contr], "_up_degs_ora.csv"))
    plots <- append(plots, list(ora_up$plot))
  }

  # split by downregulated genes, run ORA
  down_degs <- sig_degs_comb %>%
    filter(contrast == contrasts[contr], logFC < 0)

   ora_down <- run_ora_plot(down_degs$gene, paste0(contrasts[contr], " Downregulated"))
  if (!is.null(ora_up)) {
    write_csv(ora_down$result, paste0(filepath, contrasts[contr], "_down_degs_ora.csv"))
    plots <- append(plots, list(ora_up$plot))
  }
}

# save all plots to one PDF
pdf(paste0(filepath, "ORA_all_plots.pdf"), width = 10, height = 8)
for (plot in plots) {
  print(plot)
}
dev.off()

# end timer
fptm <- proc.time()
(fptm[3] / 60)

#  session info
print(sessionInfo())
