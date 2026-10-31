# script for plotting signature reversion drug figure
set.seed(42)

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(ggplot2)
  library(viridis)
  library(ggrepel)
  library(UpSetR)
  library(circlize)
  library(gprofiler2)
  library(cowplot)
})

# create parser object
parser <- ArgumentParser()

parser$add_argument("-i", "--input_dir",
  type = "character", metavar = "PATH",
  help = "Path to the directory with input files"
)
parser$add_argument("-d", "--data_dir",
  type = "character", metavar = "PATH",
  help = "Path to the data directory"
)
parser$add_argument("-o", "--output_dir",
  type = "character", metavar = "PATH",
  help = "Path to the directory to write results"
)

args <- parser$parse_args()

# start timer
ptm <- proc.time()

#### Set variables for dataset ####
indir <- args$input_dir
datadir <- args$data_dir
filepath <- args$output_dir

if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)
print(paste0("Results will be saved to ", filepath))

# ggplot theme
mytheme <- theme(
  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  panel.background = element_blank(),
  text = element_text(family = "Helvetica", face = "bold", size = 14)
)

#### READ IN DRUG DATA ####
approved_drugs <- read_csv(
  paste0(indir, "approved_drugs_res.csv"), show_col_types = FALSE
)
drug_targets_setbp1 <- read_csv(
  paste0(indir, "setbp1_drug_targets.csv"), show_col_types = FALSE
)
drug_safety <- read_csv(
  paste0(datadir, "drug_safety_summary.csv"), show_col_types = FALSE
)


# ==============================================================================
# SAFETY ANNOTATION
# ==============================================================================
# FIX: the old code built `bbw_drugs` by joining safety onto approved_drugs,
# selecting columns, then joining that BACK onto approved_drugs. That round trip
# can multiply rows and made `bbw_drugs` a data frame, which then broke
# `filter(!pref_name %in% bbw_drugs)` silently further down. Join once.

# Guard against duplicate safety rows multiplying the drug table
drug_safety_1 <- drug_safety %>%
  distinct(OT_Name, .keep_all = TRUE)

if (nrow(drug_safety_1) < nrow(drug_safety)) {
  message("Note: dropped ", nrow(drug_safety) - nrow(drug_safety_1),
          " duplicate OT_Name rows from drug_safety before joining.")
}

n_before <- nrow(approved_drugs)
approved_drugs <- approved_drugs %>%
  left_join(drug_safety_1, by = c("pref_name" = "OT_Name"))
stopifnot(nrow(approved_drugs) == n_before)

# Report drugs that found no safety record, rather than letting them sit as NA
unmatched <- approved_drugs %>%
  filter(is.na(Safety_Flag)) %>%
  distinct(pref_name) %>%
  pull(pref_name)
if (length(unmatched)) {
  message("Drugs with no safety record (", length(unmatched), "): ",
          paste(head(unmatched, 20), collapse = ", "))
}


# ==============================================================================
# SETBP1 REGULATORY TARGET FLAG  --  computed ONCE, early
# ==============================================================================
# FIX: this was the fatal error. The old script used SETBP1_Target inside the
# drug_level block but only created it ~300 lines later. It is also the source
# of the mislabelling: SETBP1_Target means "acts on any of the 361 curated
# targets", NOT "acts on one of the prioritized 34" (that is SETBP1_top_target
# in drug_targets_setbp1). Every legend below is labelled accordingly.

approved_drugs <- approved_drugs %>%
  mutate(SETBP1_Target = pert %in% drug_targets_setbp1$pert)

message("Rows: ", nrow(approved_drugs),
        " | unique drugs: ", n_distinct(approved_drugs$pref_name),
        " | rows acting on a SETBP1 regulatory target: ",
        sum(approved_drugs$SETBP1_Target))


# ------------------------------------------------------------------------------
# SAFETY FLAG -> THREE TIERS
# ------------------------------------------------------------------------------
#   Excluded      Boxed warning | Withdrawn
#   Needs review  Boxed warning (combination label only)
#                 Boxed warning (Open Targets record only)
#                 No US label found
#   Cleared       Label, no boxed warning
#                 Combination labels only, no boxed warning

safety_tier <- function(x) {
  x <- trimws(as.character(x))
  dplyr::case_when(
    x %in% c("Boxed warning", "Withdrawn")                   ~ "Excluded",
    x %in% c("Boxed warning (combination label only)",
             "Boxed warning (Open Targets record only)",
             "No US label found")                            ~ "Needs review",
    x %in% c("Label, no boxed warning",
             "Combination labels only, no boxed warning")    ~ "Cleared",
    is.na(x)                                                 ~ "Needs review",
    TRUE                                                     ~ "Needs review"
  )
}

known_flags <- c("Boxed warning", "Withdrawn",
                 "Boxed warning (combination label only)",
                 "Boxed warning (Open Targets record only)",
                 "No US label found", "Label, no boxed warning",
                 "Combination labels only, no boxed warning")
unknown <- setdiff(unique(trimws(approved_drugs$Safety_Flag)), c(known_flags, NA))
if (length(unknown)) {
  warning("Unrecognised Safety_Flag values (treated as Needs review): ",
          paste(unknown, collapse = " | "))
}

approved_drugs <- approved_drugs %>%
  mutate(Safety = factor(safety_tier(Safety_Flag),
                         levels = c("Cleared", "Needs review", "Excluded")))

print(paste0("Total number of unique drugs: ",
             n_distinct(approved_drugs$pref_name)))
approved_drugs %>% distinct(pref_name, Safety, Safety_Flag) %>%
  count(Safety, Safety_Flag) %>% print(n = 20)

write_csv(approved_drugs, paste0(filepath, "approved_drugs_safetyannotated.csv"))

# Names only, as a character vector -- this is what the old `%in% bbw_drugs`
# needed and never got.
excluded_names <- approved_drugs %>%
  filter(Safety == "Excluded") %>%
  distinct(pref_name) %>%
  pull(pref_name)


# ==============================================================================
# COLLAPSE TO ONE ROW PER DRUG
# ==============================================================================
drug_level <- approved_drugs %>%
  group_by(pref_name) %>%
  summarise(
    mean_WTCS   = mean(WTCS, na.rm = TRUE),
    median_WTCS = median(WTCS, na.rm = TRUE),
    n_cell      = n_distinct(cell),
    cells       = paste(sort(unique(cell)), collapse = ", "),
    best_pval   = min(WTCS_Pval, na.rm = TRUE),
    reg_target  = any(SETBP1_Target),
    targets     = dplyr::first(mergeTargets),
    n_annot_tar = suppressWarnings(as.numeric(dplyr::first(Ntar))),
    moa_raw     = dplyr::first(MOAclue),
    safety_raw  = dplyr::first(Safety_Flag),
    Safety      = dplyr::first(Safety),
    .groups = "drop"
  ) %>%
  mutate(
    # missing Ntar -> smallest size rather than dropped from the plot
    n_annot_tar = if_else(is.na(n_annot_tar) | n_annot_tar < 1, 1, n_annot_tar)
  )

message("Drugs: ", nrow(drug_level),
        " | acting on a SETBP1 regulatory target: ", sum(drug_level$reg_target))
drug_level %>% count(Safety) %>% print()


# ------------------------------------------------------------------------------
# MOA HIGHLIGHT PALETTE
# ------------------------------------------------------------------------------
first_moa <- function(x) {
  x <- trimws(sub(";.*$", "", as.character(x)))
  x[is.na(x) | x == "" | toupper(x) == "NA"] <- "Unannotated"
  x
}

# CHECK THIS LIST against the printout below -- spellings must match exactly or
# the class silently falls into grey.
highlight <- c("Calcium channel blocker",
               "Cyclooxygenase inhibitor",
               "HDAC inhibitor",
               "Dopamine receptor antagonist")

moa_cols <- c(
  "Calcium channel blocker"      = viridis(10)[2],
  "Cyclooxygenase inhibitor"     = viridis(10)[5],
  "HDAC inhibitor"               = viridis(10)[7],
  "Dopamine receptor antagonist" = viridis(10)[9],
  "Other / unannotated"          = "grey72"
)

message("Top MOA classes present (check your highlight spellings against these):")
drug_level %>% mutate(m = first_moa(moa_raw)) %>%
  count(m, sort = TRUE) %>% head(15) %>% print()

drug_level <- drug_level %>%
  mutate(
    MOA = first_moa(moa_raw),
    MOA = if_else(MOA %in% highlight, MOA, "Other / unannotated"),
    MOA = factor(MOA, levels = names(moa_cols))
  )
drug_level %>% count(MOA) %>% print()


# ------------------------------------------------------------------------------
# AXIS LABEL COLOURS  --  must be in factor-level order
# ------------------------------------------------------------------------------
safety_cols <- c("Cleared"      = "black",
                 "Needs review" = "grey30",
                 "Excluded"     = "grey50")

plot_dat <- drug_level %>%
  mutate(drug_f = fct_reorder(pref_name, -mean_WTCS))

# ggplot draws discrete-axis level 1 at the BOTTOM and applies an axis.text
# colour vector in level order, so arrange by the factor to match.
lab_cols <- unname(
  safety_cols[as.character(dplyr::arrange(plot_dat, drug_f)$Safety)]
)
stopifnot(length(lab_cols) == nlevels(plot_dat$drug_f), !any(is.na(lab_cols)))


# ------------------------------------------------------------------------------
# THE NEW DRUG-LEVEL PLOT
# ------------------------------------------------------------------------------
# The zero-size layer exists only to give the axis-label colours a legend entry;
# ggplot will not legend a theme element on its own.

safety_moa_plot <- ggplot(plot_dat, aes(x = -mean_WTCS, y = drug_f)) +
  geom_point(aes(fill = Safety), shape = 22, size = 0, stroke = 0,
             colour = NA, na.rm = TRUE, show.legend = TRUE) +
  geom_point(aes(colour = MOA, shape = reg_target, size = n_annot_tar),
             alpha = 0.9) +
  scale_colour_manual(values = moa_cols, drop = FALSE) +
  scale_fill_manual(values = safety_cols, drop = FALSE) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "No", `TRUE` = "Yes")) +
  scale_size(range = c(7.5, 2)) + 
  labs(
    x = "-(Mean Weighted Connectivity Score)",
    y = "Drug Name",
    colour = "Mechanism of Action",
    shape  = "Acts on a SETBP1\nregulatory target",
    size   = "Annotated\ndrug targets",
    fill   = "Safety screening\n(drug name colour)"
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 3.5)),
    shape  = guide_legend(order = 2, override.aes = list(size = 3.5)),
    size   = guide_legend(order = 3),
    fill   = guide_legend(order = 4,
                          override.aes = list(size = 4, shape = 22, stroke = 0.3))
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", face = "bold"),
    axis.text.y = element_text(size = 7.5, colour = lab_cols),
    panel.grid.major.y = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.key.height = unit(0.9, "lines")
  )

# FIX: was paste0(outdir, ...) -- `outdir` does not exist in this script
ggsave(paste0(filepath, "approved_drugs_bydrug_safety_MOA.png"),
       plot = safety_moa_plot, width = 28, height = 46, units = "cm",
       dpi = 300, limitsize = FALSE)

# Talk-sized version: top 30 by reversion strength
safety_moa_plot_top <- safety_moa_plot %+%
  (plot_dat %>% slice_min(mean_WTCS, n = 30) %>%
     mutate(drug_f = fct_drop(fct_reorder(pref_name, -mean_WTCS))))
lab_cols_top <- unname(safety_cols[as.character(
  dplyr::arrange(plot_dat %>% slice_min(mean_WTCS, n = 30) %>%
                   mutate(drug_f = fct_drop(fct_reorder(pref_name, -mean_WTCS))),
                 drug_f)$Safety)])
safety_moa_plot_top <- safety_moa_plot_top +
  theme(axis.text.y = element_text(size = 9, colour = lab_cols_top))

ggsave(paste0(filepath, "approved_drugs_bydrug_safety_MOA_top30.png"),
       plot = safety_moa_plot_top, width = 26, height = 20, units = "cm", dpi = 300)

# ==============================================================================
# ALT FULL APPROVED-DRUG PLOT: cell line + individual WTCS
# ==============================================================================
 
# One row per drug x cell line. Ntar and Safety are drug-level and ride along.
pt_dat <- approved_drugs %>%
  filter(!is.na(WTCS), !is.na(pref_name)) %>%
  mutate(
    Safety = factor(Safety, levels = names(safety_cols)),
    # level 1 is drawn at the BOTTOM, so order by the drug's BEST -WTCS
    drug_f = fct_reorder(pref_name, -WTCS, .fun = max)
  )
 
# ---- the part that breaks when you move off mean WTCS -----------------------
# axis.text.y takes one colour per LEVEL, in level order. pt_dat has one row
# per drug x cell, so indexing it directly gives 104 colours for 97 levels and
# silently shifts every label. Build the vector from a drug-level frame.
lab_dat <- pt_dat %>%
  distinct(pref_name, drug_f, Safety) %>%
  arrange(drug_f)
 
stopifnot(
  nrow(lab_dat) == nlevels(pt_dat$drug_f),            # exactly one row per level
  identical(as.character(lab_dat$drug_f),             # and in level order
            levels(pt_dat$drug_f)),
  !any(is.na(lab_dat$Safety))
)
lab_cols <- unname(safety_cols[as.character(lab_dat$Safety)])
# -----------------------------------------------------------------------------
 
# Size channel: Ntar keeps it comparable to the MOA version.
# Swap to `-log10(WTCS_Pval)` here and in labs() if you prefer the old one.
cell_wtcs_plot <- ggplot(pt_dat, aes(x = -WTCS, y = drug_f)) +
  # zero-size layer: gives the axis-label colours a legend entry, since ggplot
  # will not legend a theme element on its own
  geom_point(aes(fill = Safety), shape = 22, size = 0, stroke = 0,
             colour = NA, na.rm = TRUE, show.legend = TRUE) +
  geom_point(aes(colour = cell, shape = SETBP1_Target, size = Ntar),
             alpha = 0.85) +
  scale_colour_viridis_d(end = 0.9) +
  scale_fill_manual(values = safety_cols, drop = FALSE) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "No", `TRUE` = "Yes")) +
  scale_size(range = c(7.5, 2)) + 
  # scale_size_continuous(range = c(2, 7.5), trans = "sqrt",
  #                       breaks = c(1, 5, 20, 50, 80)) +
  labs(
    x = "-(Weighted Connectivity Score)",
    y = "Drug Name",
    colour = "Cell Line",
    shape  = "Acts on a SETBP1\nregulatory target",
    size   = "Annotated\ndrug targets",
    fill   = "Safety screening\n(drug name colour)"
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 3.5)),
    shape  = guide_legend(order = 2, override.aes = list(size = 3.5)),
    size   = guide_legend(order = 3),
    fill   = guide_legend(order = 4,
                          override.aes = list(size = 4, shape = 22, stroke = 0.3))
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", face = "bold"),
    axis.text.y = element_text(size = 7.5, colour = lab_cols),
    panel.grid.major.y = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.key.height = unit(0.9, "lines")
  )
 
ggsave(paste0(filepath, "approved_drugs_bycell_safety.png"),
       plot = cell_wtcs_plot, width = 28, height = 46, units = "cm",
       dpi = 300, limitsize = FALSE)
 
# Talk-sized: top 30 drugs by their best WTCS, all their cell lines kept
top30 <- pt_dat %>%
  group_by(pref_name) %>% mutate(best = min(WTCS)) %>% ungroup() %>%
  filter(pref_name %in% (distinct(., pref_name, best) %>%
                           slice_min(best, n = 30) %>% pull(pref_name))) %>%
  mutate(drug_f = fct_drop(fct_reorder(pref_name, -WTCS, .fun = max)))
 
lab_cols_top <- unname(safety_cols[as.character(
  (top30 %>% distinct(pref_name, drug_f, Safety) %>% arrange(drug_f))$Safety)])
stopifnot(length(lab_cols_top) == nlevels(top30$drug_f))
 
ggsave(paste0(filepath, "approved_drugs_bycell_safety_top30.png"),
       plot = (cell_wtcs_plot %+% top30) +
         theme(axis.text.y = element_text(size = 9, colour = lab_cols_top)),
       width = 26, height = 20, units = "cm", dpi = 300)
 

# ------------------------------------------------------------------------------
# COUNTS FOR THE CAPTION
# ------------------------------------------------------------------------------
cap <- drug_level %>% count(Safety, .drop = FALSE) %>% deframe()
getn <- function(k) if (k %in% names(cap)) as.integer(cap[[k]]) else 0L
message(sprintf(
  "Caption: of %d approved drugs reversing the signature, %d cleared screening, %d carry a boxed warning or withdrawal, and %d need manual review.",
  nrow(drug_level), getn("Cleared"), getn("Excluded"), getn("Needs review")
))

# Convergent drugs by safety tier -- the number that decides Act 3
convergent_by_safety <- drug_level %>%
  filter(reg_target) %>%
  select(pref_name, mean_WTCS, n_cell, n_annot_tar, targets, Safety, safety_raw) %>%
  arrange(mean_WTCS)
print(convergent_by_safety, n = 40)
write_csv(convergent_by_safety,
          paste0(filepath, "convergent_drugs_by_safety.csv"))


# ==============================================================================
# SETBP1 DRUG TARGETS
# ==============================================================================
message("Checking NA values that drop rows from the target plots...")
drug_targets_setbp1 %>%
  filter(is.na(max_phase) | is.na(WTCS) | is.na(mergeTargets)) %>%
  select(pert, pref_name, max_phase, WTCS, mergeTargets) %>%
  print()

message("Plotting shared drug/SETBP1 targets plots...")

# Annotate rather than filter, and keep safety on these rows too
drug_targets_setbp1 <- drug_targets_setbp1 %>%
  left_join(drug_safety_1 %>% select(OT_Name, Safety_Flag),
            by = c("pref_name" = "OT_Name")) %>%
  mutate(Safety = factor(safety_tier(Safety_Flag),
                         levels = c("Cleared", "Needs review", "Excluded")))

# FIX: `bbw_drugs` was a data frame, so this filter silently did nothing.
# Now a character vector. Set KEEP_EXCLUDED <- TRUE to annotate instead of drop.
KEEP_EXCLUDED <- TRUE

approved_drugs_setbp1 <- drug_targets_setbp1 %>%
  filter(as.integer(max_phase) == 4) %>%
  { if (KEEP_EXCLUDED) . else filter(., !pref_name %in% excluded_names) }

message("Convergent approved drugs: ", n_distinct(approved_drugs_setbp1$pref_name),
        " (", sum(approved_drugs_setbp1$Safety == "Excluded", na.rm = TRUE),
        " rows flagged Excluded)")

write_csv(approved_drugs_setbp1,
          paste0(filepath, "setbp1_converge_targets_safetyannotated.csv"))

simple_setbp1_drug_filtered <- approved_drugs_setbp1 %>%
  select(Drug = pref_name, Cell = cell, WTCS, Targets = mergeTargets,
         SETBP1_top_target, Safety, Safety_Flag)
write_csv(simple_setbp1_drug_filtered,
          paste0(filepath, "setbp1_converge_targets_safetyannotated_simplified.csv"))


# ==============================================================================
# ALT CONVERGENT DRUGS: target colour, target labels
# ==============================================================================
 
# ---- MOA classes ------------------------------------------------------------
# Regex over the WHOLE mergeMOA string, not the first token. first_moa() reads
# only the text before the first ";", which puts DIFLUNISAL ("Prostanoid
# receptor antagonist; Cyclooxygenase inhibitor") into grey instead of with the
# other two NSAIDs. Order matters: first match wins.
moa_class <- function(x) {
  x <- as.character(x); x[is.na(x)] <- ""
  dplyr::case_when(
    grepl("cyclooxygenase|prostanoid",   x, ignore.case = TRUE) ~ "COX inhibitor (NSAID)",
    grepl("calcium channel",             x, ignore.case = TRUE) ~ "Calcium channel blocker",
    grepl("hdac|histone deacetylase",    x, ignore.case = TRUE) ~ "HDAC inhibitor",
    grepl("kinase|\\bALK\\b|\\bABL\\b|\\bSRC\\b",
                                         x, ignore.case = TRUE) ~ "Tyrosine kinase inhibitor",
    grepl("phosphodiesterase|adenosine", x, ignore.case = TRUE) ~ "PDE / adenosine antagonist",
    grepl("bile acid",                   x, ignore.case = TRUE) ~ "Bile acid",
    TRUE                                                        ~ "Other / unannotated"
  )
}
 
conv_moa_cols <- c(
  "Calcium channel blocker"    = viridis(10)[2],
  "COX inhibitor (NSAID)"      = viridis(10)[4],
  "HDAC inhibitor"             = viridis(10)[6],
  "Tyrosine kinase inhibitor"  = viridis(10)[8],
  "PDE / adenosine antagonist" = viridis(10)[9],
  "Bile acid"                  = "#8C6D31",
  "Other / unannotated"        = "grey72"
)
 
# ---- build the plotting frame ----------------------------------------------
# MOA and Safety are drug-level: pull them once per pert, or the join multiplies
# rows by cell line.
drug_meta <- approved_drugs %>%
  distinct(pert, .keep_all = TRUE) %>%
  select(pert, mergeMOA, Safety, Ntar)
 
n_before <- nrow(approved_drugs_setbp1)
conv <- approved_drugs_setbp1 %>%
  select(pert, pref_name, cell, WTCS, Target = mergeTargets,
         Prioritized = SETBP1_top_target) %>%
  left_join(drug_meta, by = "pert")
stopifnot(nrow(conv) == n_before)
 
# One point per drug x cell. A drug hitting two targets (bosutinib: ERBB4 and
# MAP4K1) would otherwise plot two points on top of each other at the same x.

conv_pts <- conv %>%
# TODO: ADD MEAN WTCS HERE
  # group_by(pert, pref_name, cell, WTCS, mergeMOA, Safety, Ntar) %>%
  group_by(pert, pref_name, mergeMOA, Safety, Ntar) %>%
  summarise(
    Target      = paste(sort(unique(Target)), collapse = " · "),
    Prioritized = any(as.logical(Prioritized)),
    mean_WTCS   = mean(WTCS, na.rm = TRUE),
    median_WTCS = median(WTCS, na.rm = TRUE),
    n_cell      = n_distinct(cell),
    cells       = paste(sort(unique(cell)), collapse = ", "),
    .groups = "drop"
  ) %>%
  mutate(
    MOA    = factor(moa_class(mergeMOA), levels = names(conv_moa_cols)),
    Safety = factor(Safety, levels = names(safety_cols)),
    drug_f = fct_reorder(pref_name, -mean_WTCS, .fun = max),
    n_annot_tar = if_else(is.na(Ntar) | Ntar < 1, 1, Ntar)
  )
 
message("Convergent drugs: ", n_distinct(conv_pts$pref_name),
        " | points (drug x cell): ", nrow(conv_pts))
message("Check the MOA assignment -- anything in 'Other / unannotated' is a miss:")
conv_pts %>% distinct(pref_name, MOA, mergeMOA) %>% arrange(MOA) %>%
  select(pref_name, MOA) %>% print(n = 30)
 
# ---- one label per drug, at its strongest point ----------------------------
# Labelling every point double-prints the target for drugs seen in two cell
# lines (belinostat, cholic acid), and their points can sit almost on top of
# each other, so repel fights itself. One label per row instead.
conv_labs <- conv_pts %>%
  group_by(drug_f, Target) %>%
  slice_min(mean_WTCS, n = 1, with_ties = FALSE) %>%
  ungroup()
 
lab_cols_conv <- unname(safety_cols[as.character(
  (conv_pts %>% distinct(pref_name, drug_f, Safety) %>% arrange(drug_f))$Safety)])
stopifnot(length(lab_cols_conv) == nlevels(conv_pts$drug_f),
          !any(is.na(lab_cols_conv)))
 
convergent_plot <- ggplot(conv_pts, aes(x = -mean_WTCS, y = drug_f)) +
  geom_point(aes(fill = Safety), shape = 22, size = 0, stroke = 0,
             colour = NA, na.rm = TRUE, show.legend = TRUE) +
  geom_point(aes(colour = Target, shape = Prioritized, size = n_annot_tar), #size = 4.2, 
alpha = 0.9) +
  geom_text(data = conv_labs, aes(label = Target),
            hjust = 0, nudge_x = 0.004, size = 3, fontface = "bold",
            colour = "grey25") +
  # scale_colour_manual(values = conv_moa_cols, drop = FALSE) +
  scale_size(range = c(7.5, 2)) + 
  scale_colour_viridis_d(end = 0.7) +
  scale_fill_manual(values = safety_cols, drop = FALSE) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "No", `TRUE` = "Yes")) +
  # room on the right for the target labels
  scale_x_continuous(expand = expansion(mult = c(0.04, 0.22))) +
  labs(
    x = "-(Mean Weighted Connectivity Score)",
    y = "Drug Name",
    colour = "SETBP1 Regulatory Target",
    shape  = "Prioritized SETBP1 target",
    fill   = "Safety screening\n(drug name colour)"
  ) +
  guides(
    colour = guide_legend(order = 1, override.aes = list(size = 3.5)),
    shape  = guide_legend(order = 2, override.aes = list(size = 3.5)),
    fill   = guide_legend(order = 3,
                          override.aes = list(size = 4, shape = 22, stroke = 0.3))
  ) +
  theme_minimal() +
  theme(
    text = element_text(family = "Helvetica", face = "bold"),
    axis.text.y = element_text(size = 10, colour = lab_cols_conv),
    panel.grid.major.y = element_line(linewidth = 0.2),
    panel.grid.minor = element_blank(),
    legend.key.height = unit(0.9, "lines")
  )
 
ggsave(paste0(filepath, "convergent_drugs_by_target.png"),
       plot = convergent_plot, width = 26, height = 15, units = "cm", dpi = 300)
 

# FIX: assignment captured the whole pipe including ggplot(), so the object was
# a ggplot and could not be piped into ggplot() again. Data and plot split.
p_targets_hist <- approved_drugs_setbp1 %>%
  distinct(mergeTargets, n, SETBP1_target) %>%
  ggplot(aes(y = n, fill = SETBP1_target)) +
  geom_bar() +
  scale_fill_viridis_d(end = 0.7) +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "drug_targets_histogram.png"),
       plot = p_targets_hist, width = 19, height = 24, units = "cm")

# FIX: size and alpha are constants -- outside aes(), so no guides() workaround
p_trialphase <- drug_targets_setbp1 %>%
  ggplot(aes(x = -WTCS, y = max_phase, colour = mergeTargets, shape = cell)) +
  geom_point(aes(size = -log10(WTCS_Pval)), alpha = 0.8) +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "setbp1_drug_shared_targets_trialphase.png"),
       plot = p_trialphase, width = 19, height = 19, units = "cm")

p_shared_phase <- drug_targets_setbp1 %>%
  ggplot(aes(y = mergeTargets, x = -WTCS)) +
  geom_point(aes(colour = SETBP1_top_target, shape = as.factor(max_phase)),
             size = 4, alpha = 0.7) +
  scale_colour_viridis_d(end = 0.7) +
  labs(
    x = "-(Weighted Connectivity Score)",
    y = "SETBP1 & Drug Target",
    colour = "Prioritized SETBP1 Target",
    shape = "Max Drug Approval Phase"
  ) +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "setbp1_drugs_shared_phase.png"),
       plot = p_shared_phase, width = 19, height = 16, units = "cm")

# FIX: legend said "Priority SETBP1 Target" while the shape mapped to
# SETBP1_top_target; now shape = prioritized, alpha dropped (it encoded the
# same variable twice), and safety carried as the point outline.
setbp1_targets_plot <- approved_drugs_setbp1 %>%
  ggplot(aes(y = mergeTargets, x = -WTCS)) +
  geom_point(aes(colour = pref_name, shape = SETBP1_top_target),
             size = 4, alpha = 0.85) +
  scale_colour_viridis_d(end = 0.7) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "No", `TRUE` = "Yes")) +
  labs(
    x = "-(Weighted Connectivity Score)",
    y = "SETBP1 & Drug Target",
    colour = "Drug",
    shape = "Prioritized SETBP1 Target"
  ) +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "setbp1_drugs_shared.png"),
       plot = setbp1_targets_plot, width = 19, height = 16, units = "cm")


# ==============================================================================
# APPROVED DRUGS
# ==============================================================================
message("Plotting approved drug plots...")

p_multi <- approved_drugs %>%
  filter(n > 1) %>%
  ggplot(aes(x = -WTCS, y = reorder(pref_name, -drug_top_WTCS),
             colour = cell)) +
  geom_point(aes(size = -log10(WTCS_Pval)), alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = "-(Weighted Connectivity Score)", y = "Drug Name",
       colour = "Cell Line", size = "-log10(P-value)") +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "approved_drugs_braincells_multiplecells.png"),
       plot = p_multi, width = 16, height = 25, units = "cm")

p_full <- approved_drugs %>%
  ggplot(aes(x = -WTCS, y = reorder(pref_name, -drug_top_WTCS),
             colour = cell)) +
  geom_point(aes(size = -log10(WTCS_Pval)), alpha = 0.7) +
  scale_colour_viridis_d(end = 0.9) +
  labs(x = "-(Weighted Connectivity Score)", y = "Drug Name",
       colour = "Cell Line", size = "-log10(P-value)") +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "approved_drugs_braincells_full.png"),
       plot = p_full, width = 16, height = 29, units = "cm")


### TOP WTCS APPROVED DRUGS, SETBP1 TARGET ANNOTATED ###
message("Plotting top WTCS-approved, SETBP1-annotated drug plots...")

# FIX: the legend read "Priority SETBP1 Target" but SETBP1_Target means ANY
# regulatory target. Label corrected. SETBP1_Target is already on the table.
top_wtcs_plot <- approved_drugs %>%
  ggplot(aes(x = -WTCS, y = reorder(pref_name, -drug_top_WTCS))) +
  geom_point(aes(alpha = SETBP1_Target, colour = cell,
                 size = -log10(WTCS_Pval), shape = SETBP1_Target)) +
  scale_size(range = c(4, 7.5)) +
  scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.45)) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                     labels = c(`FALSE` = "No", `TRUE` = "Yes")) +
  scale_colour_viridis_d(end = 0.9) +
  labs(
    x = "-(Weighted Connectivity Score)",
    y = "Drug Name",
    shape = "Acts on a SETBP1\nregulatory target",
    size = "-log10(P-value)",
    colour = "Cell Line"
  ) +
  guides(size = "none", alpha = "none") +
  theme_minimal() +
  theme(text = element_text(face = "bold"))

ggsave(paste0(filepath, "approved_drugs_braincells_topWTCS.png"),
       plot = top_wtcs_plot, width = 19, height = 25, units = "cm")


### ANNOTATED DRUG TABLE ###
message("Writing drug tables with SETBP1 + safety annotations...")
converg_drug_anno <- approved_drugs %>% filter(SETBP1_Target)

write_csv(converg_drug_anno,
          paste0(filepath, "setbp1_converge_drug_anno_table.csv"))

# Safety is already joined, so no second join needed
write_csv(converg_drug_anno,
          paste0(filepath, "drug_full_anno_safetyannotated_table.csv"))


### MOA ###
message("Plotting approved drug MOA plots...")

# mergeMOA is the only column with no missing values (MOA is NA for 15/97,
# MOAclue for 9/97), so it is the one to count from.
MOA_COL <- "mergeMOA"
 
# Merge classes that are the same concept under different words ACROSS drugs
# too, not just within one drug. FALSE is the conservative default; see the
# printout at the end for exactly which drugs each merge would move.
GLOBAL_MERGE <- FALSE
 
 
# ------------------------------------------------------------------------------
# NORMALISATION
# ------------------------------------------------------------------------------
# Different words for the same thing. Applied to the lowercased, punctuation-
# stripped key, so "Histone deacetylase inhibitor" and "HDAC inhibitor" collide.
moa_synonyms <- c(
  "histone deacetylase"                        = "hdac",
  "peroxisome proliferator activated receptor" = "ppar",
  "retinoic acid receptor"                     = "retinoid receptor",
  "cyclooxygenase 2"                           = "cyclooxygenase",
  # the only two drugs tagged "Prostanoid receptor antagonist" are meclofenamate
  # and diflunisal, which also carry the COX token; drop this line to keep them apart
  "prostanoid receptor antagonist"             = "cyclooxygenase inhibitor"
)
 
# A "...receptor ligand" annotation is a vaguer restatement of a "...receptor
# agonist/antagonist" one on the same drug (loxapine carries both).
GENERIC_HEADS <- c("ligand", "modulator")
 
moa_key <- function(x) {
  k <- tolower(as.character(x))
  k <- gsub("[^a-z0-9]+", " ", k)
  k <- trimws(gsub("\\s+", " ", k))
  for (i in seq_along(moa_synonyms)) {
    k <- gsub(names(moa_synonyms)[i], moa_synonyms[[i]], k, fixed = TRUE)
  }
  trimws(gsub("\\s+", " ", k))
}
 
# ------------------------------------------------------------------------------
# WITHIN-DRUG HIERARCHY COLLAPSE
# ------------------------------------------------------------------------------
# Keep the MOST GENERAL wording of each concept a drug has. Class B is a more
# general form of A when B's words are a strict subset of A's:
#   {dopamine, receptor, antagonist} subset of {d2, like, dopamine, receptor, antagonist}
#   {abl, inhibitor}                 subset of {bcr, abl, kinase, inhibitor}
# Genuinely different mechanisms share no such nesting and are both kept, which
# is why clozapine still counts as dopamine AND serotonin antagonist.
drop_nested <- function(surface) {
  surface <- trimws(as.character(surface))
  surface <- surface[nzchar(surface) & !is.na(surface)]
  if (!length(surface)) return(character(0))
 
  k    <- moa_key(surface)
  keep <- nzchar(k) & !duplicated(k)        # kills case variants: SRC vs Src
  surface <- surface[keep]; k <- k[keep]
  if (length(k) <= 1) return(surface)
 
  w <- lapply(strsplit(k, " ", fixed = TRUE), unique)
  drop <- rep(FALSE, length(k))
  for (i in seq_along(k)) {
    for (j in seq_along(k)) {
      if (i == j || drop[j]) next
      if (length(w[[j]]) < length(w[[i]]) && all(w[[j]] %in% w[[i]])) {
        drop[i] <- TRUE; break
      }
      last <- w[[i]][length(w[[i]])]
      if (last %in% GENERIC_HEADS && all(setdiff(w[[i]], last) %in% w[[j]])) {
        drop[i] <- TRUE; break
      }
    }
  }
  surface[!drop]
}
 
# ------------------------------------------------------------------------------
# BUILD THE DRUG x CLASS TABLE  (no safety filtering anywhere)
# ------------------------------------------------------------------------------
moa_by_drug <- approved_drugs %>%
  distinct(pref_name, .keep_all = TRUE) %>%          # 104 rows -> 97 drugs
  select(pref_name, moa_src = all_of(MOA_COL)) %>%
  rowwise() %>%
  mutate(moa = list(drop_nested(str_split(
    if_else(is.na(moa_src) | str_trim(moa_src) %in% c("", "NA"), "Unannotated", moa_src),
    ";")[[1]]))) %>%
  ungroup() %>%
  select(pref_name, moa) %>%
  unnest(moa) %>%
  # "Unknown" is a real token (cholic acid = "Bile acid; Unknown"); keep it only
  # when it is all the drug has
  group_by(pref_name) %>%
  filter(!(moa_key(moa) == "unknown" & n() > 1)) %>%
  ungroup() %>%
  mutate(moa = if_else(moa_key(moa) == "unknown", "Unannotated", moa),
         key = moa_key(moa))
 
# one display name per concept: the most common spelling of it
moa_by_drug <- moa_by_drug %>%
  group_by(key) %>%
  mutate(MOA = names(sort(table(moa), decreasing = TRUE))[1]) %>%
  ungroup() %>%
  distinct(pref_name, key, MOA)
 
stopifnot(n_distinct(moa_by_drug$pref_name) == n_distinct(approved_drugs$pref_name))
 
# ------------------------------------------------------------------------------
# OPTIONAL: collapse the same nesting ACROSS drugs
# ------------------------------------------------------------------------------
keys <- sort(unique(moa_by_drug$key))
kw   <- lapply(strsplit(keys, " ", fixed = TRUE), unique); names(kw) <- keys
gmap <- setNames(keys, keys)
for (a in keys) {
  subs <- keys[vapply(keys, function(b)
    b != a && length(kw[[b]]) < length(kw[[a]]) && all(kw[[b]] %in% kw[[a]]),
    logical(1))]
  if (length(subs)) gmap[a] <- subs[order(lengths(kw[subs]), subs)][1]
}
changed <- names(gmap)[gmap != names(gmap)]
 
message("\nClasses still split across drugs (GLOBAL_MERGE would fold these in):")
for (a in changed) {
  added <- setdiff(moa_by_drug$pref_name[moa_by_drug$key == a],
                   moa_by_drug$pref_name[moa_by_drug$key == gmap[a]])
  if (length(added))
    message("  ", moa_by_drug$MOA[match(gmap[a], moa_by_drug$key)],
            "  +=  ", paste(sort(added), collapse = ", "))
}
message("  NOTE: 'L-type calcium channel blocker' -> 'Calcium channel blocker' adds ",
        "CARBINOXAMINE, an antihistamine whose L-type call looks like an ",
        "off-target annotation. The sodium-channel and PPAR merges are sound.\n")
 
if (GLOBAL_MERGE) {
  moa_by_drug <- moa_by_drug %>%
    mutate(key = unname(gmap[key])) %>%
    group_by(key) %>%
    mutate(MOA = names(sort(table(MOA), decreasing = TRUE))[1]) %>%
    ungroup() %>%
    distinct(pref_name, key, MOA)
}
 
# ------------------------------------------------------------------------------
# THE NUMBERS TO QUOTE
# ------------------------------------------------------------------------------
moa_counts <- moa_by_drug %>%
  count(MOA, name = "Drugs") %>%
  arrange(desc(Drugs), MOA)
 
message("Drugs: ", n_distinct(moa_by_drug$pref_name),
        " | drug x class rows: ", nrow(moa_by_drug),
        " | distinct classes: ", n_distinct(moa_by_drug$key),
        " | mean classes per drug: ",
        round(nrow(moa_by_drug) / n_distinct(moa_by_drug$pref_name), 2))
message("A drug can hold several classes, so Drugs does NOT sum to 97.")
print(moa_counts %>% filter(Drugs >= 3), n = 30)
 
write_csv(moa_counts,  paste0(filepath, "moa_counts_by_drug.csv"))
write_csv(moa_by_drug, paste0(filepath, "moa_by_drug_long.csv"))
 
# sanity: no drug counted twice in nested classes
dupes <- moa_by_drug %>% count(pref_name, key) %>% filter(n > 1)
stopifnot(nrow(dupes) == 0)
 
# ------------------------------------------------------------------------------
# HISTOGRAM  --  bars are drugs, so the axis number is the number you can say
# ------------------------------------------------------------------------------
moa_hist <- moa_by_drug %>%
  group_by(MOA) %>% filter(n() >= 3) %>% ungroup() %>%
  ggplot(aes(y = reorder(MOA, MOA, length))) +
  geom_bar(fill = "grey35", width = 0.75) +
  labs(y = "Mechanism of action", x = "Approved drugs") +
  theme_minimal() +
  theme(text = element_text(face = "bold"),
        panel.grid.major.y = element_blank())
 
ggsave(paste0(filepath, "moa_counts_by_drug.png"),
       plot = moa_hist, width = 22, height = 14, units = "cm", dpi = 300)
 
# for safety split later, join it on here -- it does not touch the
# counts above:
#   moa_by_drug %>%
#     left_join(distinct(approved_drugs, pref_name, Safety), by = "pref_name") %>%
#     count(MOA, Safety) %>% pivot_wider(names_from = Safety, values_from = n,
#                                        values_fill = 0)


#### FIGURE 5 ASSEMBLY ####
message("Plotting figure 5...")
png(paste0(filepath, "figure_5.png"),
    width = 350, height = 425, units = "mm", res = 300)
print(cowplot::plot_grid(
  moa_hist, top_wtcs_plot, setbp1_targets_plot,
  labels = "AUTO", rel_heights = c(1, 1.5, 1), ncol = 1
))
dev.off()

# Poster figure 4
message("Plotting poster figure 4...")
png(paste0(filepath, "poster_figure_4.png"),
    width = 250, height = 600, units = "mm", res = 300)
print(cowplot::plot_grid(
  top_wtcs_plot, setbp1_targets_plot,
  labels = "AUTO", rel_heights = c(2.5, 1), ncol = 1
))
dev.off()

# end timer
fptm <- proc.time()
(fptm[3] / 60)

print(sessionInfo())

