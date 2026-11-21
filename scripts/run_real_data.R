#!/usr/bin/env Rscript
# Usage: Rscript scripts/run_real_data.R cpu 1
# Produces one 10-fold CV repeat across 9 historical knot budget settings.
source("R/common.R")
source("R/data_prep.R")
root <- project_root()
check_hal()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript scripts/run_real_data.R <dataset> <repeat 1..10>")
name <- args[[1L]]
repeat_id <- as.integer(args[[2L]])
datasets <- c("cpu", "automobile", "communities_crime", "pima", "oecdpanel",
              "bike_sharing", "liver_disorders", "laheart", "real_estate")
if (!(name %in% datasets) || is.na(repeat_id) || repeat_id < 1L) stop("Invalid dataset or repeat")
group_id <- match(name, datasets)
dat <- load_dataset(name)
stopifnot(nrow(dat$X) == length(dat$Y), !anyNA(dat$X), !anyNA(dat$Y))
# Same 10 folds are used for all methods and knot settings within each repeat.
fold_seed <- seed_for(group_id, repeat_id, 1L)
fold_ids <- make_balanced_folds(nrow(dat$X), 10L, seed = fold_seed)
out_path <- file.path(root, "results", "real_data", sprintf("%s_repeat_%02d.csv", name, repeat_id))
rows <- load_previous(out_path)
assert_checkpoint_version(rows, out_path)
all_knot_settings <- seq(20L, 100L, by = 10L)
# Optional smoke test: e.g. HAL_KNOT_STARTS=20 (otherwise all 9 values).
if (nzchar(Sys.getenv("HAL_KNOT_STARTS", ""))) {
  settings <- as.integer(strsplit(Sys.getenv("HAL_KNOT_STARTS"), ",", fixed = TRUE)[[1]])
  if (!all(settings %in% all_knot_settings)) stop("Invalid HAL_KNOT_STARTS")
} else settings <- all_knot_settings

for (fold in 1:10) {
  train_idx <- which(fold_ids != fold)
  valid_idx <- which(fold_ids == fold)
  train_X <- dat$X[train_idx, , drop = FALSE]
  train_Y <- dat$Y[train_idx]
  selected <- select_predictors(train_X, train_Y, dat$selection,
                                seed = seed_for(group_id, repeat_id, fold, 1L))
  train_X <- train_X[, selected, drop = FALSE]
  valid_X <- dat$X[valid_idx, selected, drop = FALSE]
  valid_Y <- dat$Y[valid_idx]
  n_cont <- sum(vapply(seq_len(ncol(train_X)), function(j) {
    length(unique(train_X[, j])) > 2L
  }, logical(1)))

  for (k0 in settings) {
    knot_label <- sprintf("%d,10,5", k0)
    if (completed_key(rows, c("fold", "knot_start"), c(fold, k0))) next
    set.seed(seed_for(group_id, repeat_id, fold, k0))
    message(sprintf("%s repeat=%d fold=%d knots=%s", name, repeat_id, fold, knot_label))
    q <- fit_hal_screened(train_X, train_Y, valid_X, "quantile", 3L, c(k0, 10, 5))
    budget <- clustering_budget(q$basis_list, n_continuous = n_cont,
                                n_columns = ncol(train_X))
    k <- fit_hal_screened(train_X, train_Y, valid_X, "kmeans", 3L, budget)
    p <- fit_hal_screened(train_X, train_Y, valid_X, "pam", 3L, budget)
    fits <- list(quantile = q, kmeans = k, pam = p)
    checkpoint <- do.call(rbind, lapply(names(fits), function(method) {
      f <- fits[[method]]
      data.frame(analysis_version = analysis_version, dataset = name, repeat_id = repeat_id, fold = fold,
                 knot_start = k0, knot_setting = knot_label, method = method,
                 mse = mean((valid_Y - f$prediction)^2),
                 total_bases = f$total_bases, active_bases = f$active_bases,
                 selected_predictors = ncol(train_X),
                 n_train = nrow(train_X), n_valid = length(valid_Y),
                 seed = seed_for(group_id, repeat_id, fold, k0),
                 stringsAsFactors = FALSE)
    }))
    if (!all(is.finite(checkpoint$mse))) stop("Non-finite MSE; inspect predictions")
    if (!is.null(rows)) {
      rows <- rows[!(rows$fold == fold & rows$knot_start == k0), , drop = FALSE]
    }
    rows <- append_checkpoint(rows, checkpoint, out_path)
    rm(q, k, p, fits)
    invisible(gc(verbose = FALSE))
  }
}
message("Completed ", out_path)
