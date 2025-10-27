suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(magrittr)
  library(e1071)
  library(gbm) # for gbm
})

##### calc_nonzero #####
# fxn removes genes with zero variance across samples
# counts = counts dataframe/matrix
calc_nonzero <- function(counts) {
  nonzero_cts <- counts %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rowwise() %>%
    filter(var(c_across(-rowname)) != 0) %>%
    ungroup() %>%
    column_to_rownames(var = "rowname")
}


##### get_imp_feats #####
# pcameta = df of PC's/patterns + metadata aggrgated,
# modeltype = model type to use for finding variable importance (gbm, svm, rf)
# suffix = suffix to look for for PC/pattern columns
# numpc = number of PC's or patterns,
# featlist = list of 3 metadata columns to use for model recipe, and
# filepath = filepath for saving PDF of feature importance plots for each PC
# function returns a list of feature importance info for each PC
get_imp_feats <- function(pcameta, modeltype, suffix, numpc, sra, filepath) {
  vi_list <- NULL

  for (i in 1:numpc) {
    temp <- pcameta %>%
      select(PC = i, !starts_with(suffix))

    # make model fit obj
    set.seed(42)
    # Support Vector Machine Model
    if (modeltype == "svm") {
      mod <- svm(PC ~ .,
        data = temp,
        probability = TRUE
      )
    } else if (modeltype == "gbm") {
    # Generalized Boosted Model
      mod <- gbm(PC ~ .,
        data = temp,
        n.trees = 500,
        interaction.depth = 1,
        n.cores = 2
      )
    } else if (modeltype == "rf") {
    # Random Forest Model
      mod <- ranger(PC ~ .,
        data = temp,
        classification = FALSE
      )
    } else {
      stop("Model type not available! Please use 'svm', 'gbm', or 'rf'")
    }
    # create explainer for PC/pattern
    explained <- explain(mod,
      data = temp[, -1], # all metadata cols as predictors
      y = temp[, 1], # PC as outcome to pred in col 1
      label = paste0(modeltype, " model for ", sra, suffix, i)
    )
    vi <- model_parts(explained)

    # variable importance
    set.seed(1)
    vi <- model_parts(explained)
    vi_list[[i]] <- vi
  }
  # set filename for plots and object
  savename <- paste0(filepath, "feature_importance_", modeltype)

  # plot variable importance
  pdf(paste0(savename, "_plots.pdf"))
  print(plot(vi_list))
  dev.off()
  # file save message
  print(paste0("Plot saved to: ", savename, "_plots.pdf"))

  # save results
  saveRDS(vi_list, paste0(savename, ".rds"))
  # file save message
  print(paste0("Explainer object saved to: ", savename, ".rds"))
  return(vi_list)
}

##### rse_pca accepts #####
# pcameta = dataframe (with PCs and metadata)
# var_explained = PCs meeting variance explained cutoff
# color_by = desired pcameta metadata variable to color PC's by
# sra = SRA number or some other identifier for the dataset for plot title
# filepath = filepath for saving PCA plot
rse_pca <- function(
    pcameta, var_explained, color_by, shape = 16, sra, filepath) {
  pdf(filepath)
  for (i in 1:(nrow(var_explained) - 1)) {
    set.seed(42)

    pca_plot <- pcameta %>%
      # color by discrete values
      ggplot(
        aes(
          x = pcameta[, i],
          y = pcameta[, i + 1],
          color = pcameta[, color_by],
          shape = pcameta[, shape]
        )
      ) + {
        if (class(pcameta[, color_by]) == "factor") {
          scale_colour_viridis_d(
            name = paste0(color_by), # fix names and wrap long names
            labels = function(x) str_wrap(x, width = 15)
          )
        }
      } + {
        # color by continuous values, check for discrete vs continuous
        if (class(pcameta[, color_by]) != "factor") {
          scale_colour_viridis_c(
            # wrap long names
            name = paste0(color_by),
            labels = function(x) str_wrap(x, width = 15)
          )
        }
      } +
      geom_point(size = 6, alpha = 0.6) +
      xlab(paste0(
        "PC", i, " ",
        # % variance explained by PC
        round(var_explained[i, 2] * 100, digits = 2),
        "% Variance"
      )) +
      ylab(paste0(
        "PC", i + 1, " ",
        # % variance explained by PC
        round(var_explained[i + 1, 2] * 100, digits = 2),
        "% Variance"
      )) +
      labs(title = str_wrap(paste0("PCA of ", sra), width = 23)) +
      theme(
        panel.background = element_blank(),
        text = element_text(face = "bold", size = 12)
      )

    print(pca_plot)
  }
  dev.off()
}

##### pca_topvar_plot #####
# imp_feats_res = DALEX object of feature importance from get_imp_feats
# pcameta = df of PC's + metadata aggrgated,
# var_explained = a df with columns "PC" and "var_explained" values
# sra = sra ID for data set being analyzed
# study = study title for data set being analyzed
# modeltype = model used with explainer
# filepath = filepath to save PCA's
# function plots and saves PCA's colored by variable importance
pca_topvar_plot <- function(
    imp_feats_res, pcameta, var_explained, sra, modeltype, filepath) {
  # plot
  savefile <- paste0(filepath, "/PCA_coloredby_imp_feat_", modeltype, ".pdf")
  pdf(savefile)

  for (i in 1:(length(imp_feats_res) - 1)) {
    # find top variable by median dropout loss
    top_v <- imp_feats_res[[i]] %>%
      group_by(variable) %>%
      # find highest median
      summarise(Median_DL = median(dropout_loss)) %>%
      filter(variable %in% colnames(pcameta)) %>%
      arrange(desc(Median_DL))

    top_v <- top_v$variable[1] # name of top variable for PC

    pca_plot <- pcameta %>%
      # color by discrete values
      ggplot(
        aes(
          x = pcameta[, i],
          y = pcameta[, i + 1],
          color = pcameta[, top_v]
        )
      ) + {
      if (class(pcameta[, top_v]) == "factor") {
          scale_color_discrete(
            name = paste0(top_v), # fix names and wrap long names
            labels = function(x) str_wrap(x, width = 15)
          )
        }
      } + {
        # color by continuous values
        if (class(pcameta[, top_v]) != "factor") {
          # check for discrete vs continuous
          scale_color_continuous(
            name = paste0(top_v), # fix names and wrap long names
            labels = function(x) str_wrap(x, width = 15)
          )
        }
      } +
      geom_point(size = 4) +
      xlab(paste0(
        "PC", i, " ",
        # % variance explained by PC
        round(var_explained[i, 2] * 100, digits = 2),
        "% Variance"
      )) +
      ylab(paste0(
        "PC", i + 1, " ",
        # % variance explained by PC
        round(var_explained[i + 1, 2] * 100, digits = 2),
        "% Variance"
      )) +
      labs(title = str_wrap(paste0("PCA of ", sra, ": ")), width = 23)

    print(pca_plot)
    # PC "colored by" message
    print(paste0("PC", i, " Colored by most important variable: ", top_v))
  }
  dev.off()
  # file save message
  print(paste0("Plot saved to: ", savefile))
}

##### simplify_cts #####
# reformats RSE w/ nonzero var genes into a simplified dataframe with
# rownames as columns
# checks for adj_cts assay and uses it be default, otherewise uses TPM
simplify_cts <- function(rse) {
  # check if there are cts adjustments -- use TPM if not
  if (is.null(assays(rse)$adj_cts)) {
    cts <- assays(rse)$TPM
  } else {
    cts <- assays(rse)$adj_cts
  }
  cts <- cts %>%
    as.data.frame() %>%
    rownames_to_column()
}

##### top_quartile #####
# accepts RSE and returns dataframe of top quartile variable genes
# checks for adj_cts assay and uses it be default, otherewise uses TPM
top_quartile <- function(rse) {
  # check if there are cts adjustments -- use TPM if not
  if (is.null(assays(rse)$adj_cts)) {
    cts <- assays(rse)$TPM
  } else {
    cts <- assays(rse)$adj_cts
  }
  var_genes <- cts %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    rowwise() %>%
    mutate(Variance = var(c_across(-rowname))) %>%
    ungroup()

  var_quart <- quantile(var_genes$Variance)
  top_quart <- var_genes %>%
    filter(Variance >= var_quart[4]) %>%
    select(!Variance)
}
