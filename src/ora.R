suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(recount3)
  library(ggplot2)
  library(ComplexHeatmap)
  library(limma)
  library(edgeR)
  # source(paste0(getwd(), "/functions.R"))
  # source(paste0(getwd(), "/src/functions.R"))
})

sig_degs <- read_csv("results/dea/shaw/shaw/combined_contrasts_sig_degs.csv")
combined_contrasts_labeled_degs <- read_csv("results/dea/shaw/shaw/combined_contrasts_labeled_degs.csv")

# compare to neuronal markers used in original study
markers <- c("POU5F1", "PAX6", "NES")
marker_degs <- combined_contrasts_labeled_degs %>% 
  filter(gene_name %in% markers) %>%
  arrange(gene_name)

write.csv(marker_degs, "./results/dea/shaw/shaw/neuronal_marker_degs.csv")

ipsc <- sig_degs %>%
  filter(contrast %in% c("IPSC_PATH", "IPSC_HD"))


sig_degs_wide <- sig_degs %>%
  select(gene, gene_name, logFC, contrast) %>%
  pivot_wider( names_from = contrast, values_from = logFC)

npc_path <- sig_degs %>% filter(contrast == "NPC_HD")
npc_path_ora <- gost(npc_path$gene, sources = c("GO"), evcodes = TRUE)

filtered_ora <- npc_path_ora[["result"]] %>%
  filter(term_size < 1000)

filtered_ora %>%
  ggplot(aes(x = intersection_size, y = recall, color = source)) +
  geom_point(aes(size = -log(p_value))) +
  geom_label(
    data = subset(filtered_ora, recall > 0.18),
    nudge_x = 6,
    nudge_y = 0.02,
    check_overlap = TRUE,
    aes(label = term_name)
  )
ggsave("./temp.png")



