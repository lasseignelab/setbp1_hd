suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
})

# create parser object
parser <- ArgumentParser()

# specify desired options
parser$add_argument("-f", "--functions",
  type = "character",
  metavar = "PATH",
  default = NULL,
  help = "Path to functions file to source custom functions" # nolint
)
parser$add_argument("-c", "--counts",
  type = "character",
  metavar = "PATH",
  help = "Counts summarized experiment file (rds format) from nf-core outputs"
)
parser$add_argument("-m", "--meta",
  type = "character",
  metavar = "PATH",
  help = "Metadata for dataset used in --counts"
)
parser$add_argument("-n", "--name",
  type = "character",
  metavar = "NAME",
  help = "Dataset name  for labeling purposes" # nolint
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

rse <- readRDS(args$counts)
meta <- read.csv(args$meta, row.names = 1, stringsAsFactors = TRUE)
name <- args$name
# rse <- readRDS("./data/cardo/star_salmon/salmon.merged.gene_counts.SummarizedExperiment.rds")
# sra <- "cardo"
# meta <- read.csv("./doc/cardo_metadata.csv")
# create output dir if needed
filepath <- paste0(args$output_dir)
if (!dir.exists(filepath)) dir.create(filepath, recursive = TRUE)

#### Metadata ####
# dataset-specific adjustments for easier downstream analyses

if (name == "cardo") {
  meta_filt <- meta %>%
    select(
      !starts_with("DATASTORE") | c(
        Bytes, Consent, ReleaseDate, create_date
      )
    ) %>%
    mutate(
      # rm whitespace
      across(everything(), ~ str_replace_all(., "\\s+", "_")),
      cells = ifelse(cell_type == "hESC-derived_neurons", "Neuron", "NPC"),
      sample_type = paste0(cells, "_", genotype),
      sample = paste0(genotype, "_", differentiation_days),
      # broad identifier of SETBP1-KOs vs control
      condition = ifelse(
        grepl("WT", sample, ignore.case = TRUE), "Control", "SETBP1-HD"
      )
    ) %>%
    rownames_to_column() %>%
    group_by(sample) %>%
    mutate(rep = seq_along(sample)) %>%
    column_to_rownames()

  print(meta_filt)
}

if (name == "shaw") {
  meta_filt <- meta %>%
    select(
      !starts_with("DATASTORE") | c(
        Bytes, Collection_Date, Consent, ReleaseDate, create_date
      )
    ) %>%
    filter(rownames(.) %in% colnames(rse)) %>%
    mutate(
      # rm whitespace
      across(everything(), ~ str_replace_all(., "\\s+", "_")),
      cells = ifelse(cell_type == "neural_progenitor_cells", "NPC", "IPSC"),
      sample = paste0(cells, "_", genotype),
      ACMG = ifelse(
        genotype == "WT", "WT", ifelse(
          genotype == "VUS2", "VUS", "PATH"
        )
      ),
      # broad identifier of SETBP1-KOs vs control
      condition = ifelse(
        grepl("WT", sample, ignore.case = TRUE), "Control", "SETBP1-HD"
      )
    ) %>%
    rownames_to_column() %>%
    group_by(sample) %>%
    mutate(rep = seq_along(sample)) %>%
    column_to_rownames()

  print(meta_filt)
}

if (name == "tanaka") {
  meta_filt <- meta %>%
    select(
      !starts_with("DATASTORE") | c(
        Bytes, Collection_Date, Consent, ReleaseDate, create_date
      )
    ) %>%
    mutate(
      # rm whitespace
      across(everything(), ~ str_replace_all(., "\\s+", "_")),
      # rm prefix from samples
      sample = str_replace(`Sample.Name`, "^RNA-seq,_LSK,_", ""),
      rep = as.integer(str_extract(sample, "[0-9]+$")),
      source_name = ifelse(
        grepl("^RNA-seq,_LSK,_", `Sample.Name`), "LSK", "LS"
      ),
      # split sample names and rep numbers
      sample = str_replace(sample, "-", ""),
      sample_name = str_extract(sample, ".*(?=_[0-9]+$|[0-9]+$)"),
      sample_type = ifelse(
        sample_name == "KO", "Setbp1_KO",
        ifelse(sample_name == "control", "Control", sample_name)
      ),
      cells = "HSC",
      # broad identifier of SETBP1-KOs vs control
      condition = ifelse(
        grepl("control", sample, ignore.case = TRUE), "Control", "SETBP1-HD"
      )
    )

  print(meta_filt)
}

# save filtered/refined metadata
write_csv(meta_filt, paste0(filepath, name, "_filtered_meta.csv"))

##### Remove Ensembl Versions from Gene IDs #####
print(paste0("Removing Ensembl versions from ID's... ", head(rownames(rse))))
rownames(rse) <- gsub("\\.\\d+$", "", rownames(rse))

# update rownames in individual assays as well
for (i in seq_along(assays(rse))) {
  rownames(rse@assays@data@listData[[i]]) <- rownames(rse)
}

print(paste0("Adding metadata to RSE...", str(rse)))
rse@metadata <- meta_filt
print(
  paste0(
    "Saving full prepped (with metadata and non-versioned gene IDs) RSE as: ",
    name, "_prepped_rse.rds"
  )
)
write_rds(rse, paste0(filepath, name, "_prepped_rse.rds"))

##### Detect/Filter Gene Counts #####
# filtering for 0 variance genes
cts_filtered <- calc_nonzero(
  rse@assays@data@listData[["salmon.merged.gene_tpm"]]
)
rse_filtered <- rse[rownames(cts_filtered), ]

print(
  paste0(
    "Saving prepped RSE filtered for nonzero-variance genes as: ",
    name, "_filtered_rse.rds"
  )
)
write_rds(rse_filtered, paste0(filepath, name, "_filtered_rse.rds"))
