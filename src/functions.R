suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(magrittr)
  library(e1071)
  library(gbm)
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
    }
    # Generalized Boosted Model
    else if (modeltype == "gbm") {
      mod <- gbm(PC ~ .,
        data = temp,
        n.trees = 500,
        interaction.depth = 1,
        n.cores = 2
      )
    }
    # Random Forest Model
    else if (modeltype == "rf") {
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
# color_by = desired metadata variable to color PC's by (must be pcameta colname)
# sra = SRA number or some other identifier for the dataset for plot title
# filepath = filepath for saving PCA plot
rse_pca <- function(pcameta, var_explained, color_by, shape = 16, sra, filepath) {
  pdf(filepath)
  for (i in 1:(nrow(var_explained) - 1)) {
    set.seed(42)

    pca_plot <- pcameta %>%
      # color by discrete values
      ggplot(
        aes(
          x = pcameta[, i], y = pcameta[, i + 1], color = pcameta[, color_by], shape = pcameta[, shape]
        )
      ) +
      {
        if (class(pcameta[, color_by]) == "factor") {
          scale_colour_viridis_d(
            name = paste0(color_by), # fix names and wrap long names
            labels = function(x) str_wrap(x, width = 15)
          )
        }
      } +
      # color by continuous values
      {
        if (class(pcameta[, color_by]) != "factor") { # check for discrete vs continuous
          scale_colour_viridis_c(
            name = paste0(color_by), # fix names and wrap long names
            labels = function(x) str_wrap(x, width = 15)
          )
        }
      } +
      geom_point(size = 6, alpha = 0.6) +
      xlab(paste0(
        "PC", i, " ",
        round(var_explained[i, 2] * 100, digits = 2), # % variance explained by PC
        "% Variance"
      )) +
      ylab(paste0(
        "PC", i + 1, " ",
        round(var_explained[i + 1, 2] * 100, digits = 2), # % variance explained by PC
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
# var_explained = var_explained_cutoff, a df with columns "PC" and "var_explained" values
# sra = sra ID for data set being analyzed
# study = study title for data set being analyzed
# modeltype = model used with explainer
# filepath = filepath to save PCA's
# function plots and saves PCA's colored by variable importance
pca_topvar_plot <- function(imp_feats_res, pcameta, var_explained, sra, modeltype, filepath) {
  # plot
  savefile <- paste0(filepath, "/PCA_coloredby_imp_feat_", modeltype, ".pdf")
  pdf(savefile)

  for (i in 1:(length(imp_feats_res) - 1)) {
    # find top variable by median dropout loss
    top_v <- imp_feats_res[[i]] %>%
      group_by(variable) %>%
      summarise(Median_DL = median(dropout_loss)) %>% # find highest median
      filter(variable %in% colnames(pcameta)) %>% # rm model info
      arrange(desc(Median_DL)) # reorder variables by median

    top_v <- top_v$variable[1] # name of top variable for PC

    pca_plot <- pcameta %>%
      # rse_pca(pcameta, var_explained, color_by = top_v, sra, study, filepath)
      # color by discrete values
      ggplot(aes(x = pcameta[, i], y = pcameta[, i + 1], color = pcameta[, top_v])) +
      {
        if (class(pcameta[, top_v]) == "factor") {
          scale_color_discrete(
            name = paste0(top_v), # fix names and wrap long names
            labels = function(x) str_wrap(x, width = 15)
          )
        }
      } +
      # color by continuous values
      {
        if (class(pcameta[, top_v]) != "factor") { # check for discrete vs continuous
          scale_color_continuous(
            name = paste0(top_v), # fix names and wrap long names
            labels = function(x) str_wrap(x, width = 15)
          )
        }
      } +
      geom_point(size = 4) +
      xlab(paste0(
        "PC", i, " ",
        round(var_explained[i, 2] * 100, digits = 2), # % variance explained by PC
        "% Variance"
      )) +
      ylab(paste0(
        "PC", i + 1, " ",
        round(var_explained[i + 1, 2] * 100, digits = 2), # % variance explained by PC
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


##### prep_exerc_data #####
# fxn renames exercise column, corrects for known batch effect (if any),
# removes zero var genes,  adds log transformations, and adds lab name to rse
# rse = recount3 RSE object
# batch = colname indicating known batch effect to correct for (if any)
# lab = adds lab name to RSE to help distinguish datasets (if)
prep_exerc_data <- function(rse, exercise_col, exercise_pre_term = NULL, interact_cov = NULL, batch = NULL, lab = NULL) {
  # rename exercise column to "exercise" for consistency across data sets
  names(rse@colData)[names(rse@colData) == exercise_col] <- "sra_attribute.exercise"

  # remove sample_attribute column -- aggregate info across the other cols
  rse@colData@listData$sra.sample_attributes <- NULL

  # rm zero var genes and adjust rse
  cts <- calc_nonzero(assays(rse)$TPM)
  rse <- rse[rownames(rse) %in% rownames(cts), ]

  # LIMMA correct for batch affects if provided
  if (!is.null(batch)) {
    print(paste0("Adjusting for ", batch, " and log transforming counts"))
    print(paste0("Ensure non 'exercise_pre_term' options are post-exercise"))
    print(paste0(unique(rse@colData@listData$sra_attribute.exercise)))
    # coerce pre/post terminology to be consistent across data sets
    prepost <- as.data.frame(rse@colData@listData) %>%
      transmute(
        group = ifelse(
          sra_attribute.exercise == paste0(exercise_pre_term), "pre", "post"
        )
      )
    exercise <- prepost$group

    # make design
    des <- model.matrix(~ 0 + exercise)

    # batch correction
    if (is.null(interact_cov)) {
      assay(rse, "adj_cts") <- removeBatchEffect(
        assays(rse)$TPM,
        batch = rse@colData[[paste0(batch)]],
        design = des
      )
    } else {
      assay(rse, "adj_cts") <- removeBatchEffect(
        assays(rse)$TPM,
        batch = rse@colData[[paste0(batch)]],
        covariates = rse@colData[[paste0(interact_cov)]],
        design = des
      )
    }

    print(paste0("Corrected counts stored in 'adj_cts' of assay(rse)"))
    # log transformed adjusted counts
    assay(rse, "log_adj") <- log((1 + assays(rse)$adj_cts), base = 2)
    print(
      paste0(
        "Log-transformed adjusted counts stored in 'log_adj' of assay(rse)"
      )
    )
  }
  # log-transform counts as recommended in projectR paper
  assay(rse, "logTPM") <- log((1 + assays(rse)$TPM), base = 2)
  print(paste0("Log-transformed TPM counts stored in 'logTPM' of assay(rse)"))
  # add lab name for easier referencing
  rse$lab <- paste0(lab)

  return(rse)
}

##### simplify_cts #####
# reformats RSE w/ nonzero var genes into a simplified dataframe with rownames as columns
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
# accepts RSE and returns dataframe of top quartile variable genes across dataset
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

##### dea_exercise #####
# accepts RSE, count type to use from RSE, and pre-exercise term to run LIMMA on pre/post exercise samples and convert ens genes to symbols
# exercise_pre_term = how "pre" exercise metadata is listed in sra_attribute.exercise of RSE
dea_exercise <- function(rse, cts, exercise_pre_term) {
  # check options
  print(paste0(
    "Ensure non 'exercise_pre_term' options are referring to post-exercise"
  ))
  print(paste0(nique(rse@colData@listData$sra_attribute.exercise)))
  # coerce pre/post terminology to be conistent acros data sets
  prepost <- as.data.frame(rse@colData@listData) %>%
    # column_to_rownames(var = "external_id") %>%
    transmute(
      group = ifelse(
        sra_attribute.exercise == paste0(exercise_pre_term), "pre", "post"
      )
    )
  exercise <- prepost$group

  # make design
  des <- model.matrix(~ 0 + exercise)

  # contrast pre vs post exercise
  contr <- makeContrasts(
    preVSpost = exercisepre - exercisepost, levels = colnames(des)
  )

  # LIMMA and empirical bayes, take top 100 DEGs
  fit <- lmFit(rse@assays@data@listData[[paste0(cts)]], des)
  fit <- contrasts.fit(fit, contrasts = contr)
  efit <- eBayes(fit)
  degs <- topTreat(efit, n = 500, p.value = 0.05)

  # ens to symbol conversion
  t <- gconvert(rownames(degs))
  degs <- degs %>%
    rownames_to_column(var = "input") %>%
    inner_join(t, c = "input")
}
