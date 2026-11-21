#!/usr/bin/env Rscript
# Usage: Rscript scripts/run_simulation.R regular 1
# Produces ONE independent Monte Carlo replicate (9 DGPs x 3 methods).
source("R/common.R")
source("R/simulation_config.R")
root <- project_root()
check_hal()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript scripts/run_simulation.R <group> <repeat 1..10>")
group <- args[[1L]]
repeat_id <- as.integer(args[[2L]])
if (is.na(repeat_id) || repeat_id < 1L) stop("Invalid repeat ID")
group_id <- match(group, simulation_groups$group)
if (is.na(group_id)) stop("Unknown group ", group)
spec <- simulation_groups[group_id, , drop = FALSE]
source(file.path("R", "generators", spec$generator))

out_path <- file.path(root, "results", "simulations", sprintf("%s_repeat_%02d.csv", group, repeat_id))
rows <- load_previous(out_path)
assert_checkpoint_version(rows, out_path)
# Historical settings: n=1000, max interaction degree=3, quantile knots 50/15/5.
for (d in c(1L, 3L, 5L)) {
  for (shape in c("smooth", "jump", "sin")) {
    # Use an explicit formula independent of the progress of prior tasks.
    setting_id <- match(d, c(1L, 3L, 5L)) * 3L + match(shape, c("smooth", "jump", "sin")) - 3L
    if (completed_key(rows, c("dimension", "scenario"), c(d, shape))) next
    this_seed <- seed_for(group_id, repeat_id, setting_id)
    set.seed(this_seed)
    message(sprintf("%s repeat=%d d=%d shape=%s", group, repeat_id, d, shape))
    train <- get_simulation_dataset(spec, d, shape, 1000L)
    test <- get_simulation_dataset(spec, d, shape, spec$test_n)
    stopifnot("Y" %in% names(train), "Y" %in% names(test))
    X <- train[, setdiff(names(train), "Y"), drop = FALSE]
    Xtest <- test[, colnames(X), drop = FALSE]
    Y <- train$Y
    Ytest <- test$Y

    q <- fit_hal_screened(X, Y, Xtest, "quantile", 3L, c(50, 15, 5))
    k_budget <- clustering_budget(q$basis_list, n_continuous = d, n_columns = d)
    k <- fit_hal_screened(X, Y, Xtest, "kmeans", 3L, k_budget)
    p <- fit_hal_screened(X, Y, Xtest, "pam", 3L, k_budget)
    fits <- list(quantile = q, kmeans = k, pam = p)
    one_setting <- do.call(rbind, lapply(names(fits), function(method) {
      fit <- fits[[method]]
      data.frame(analysis_version = analysis_version, group = group, repeat_id = repeat_id, dimension = d,
                 scenario = shape, method = method,
                 r2 = 1 - sum((Ytest - fit$prediction)^2) /
                     sum((Ytest - mean(Ytest))^2),
                 active_bases = fit$active_bases,
                 total_bases = fit$total_bases,
                 seed = this_seed, n_train = nrow(X), n_test = nrow(Xtest),
                 stringsAsFactors = FALSE)
    }))
    if (!all(is.finite(one_setting$r2))) stop("Non-finite R-squared; examine DGP and predictions")
    # Append exactly one complete setting; if interrupted, the last successful
    # setting is available and can be resumed without changing random seeds.
    if (!is.null(rows)) {
      rows <- rows[!(rows$dimension == d & rows$scenario == shape), , drop = FALSE]
    }
    rows <- append_checkpoint(rows, one_setting, out_path)
    rm(q, k, p, fits, train, test, X, Xtest, Y, Ytest)
    invisible(gc(verbose = FALSE))
  }
}
message("Completed ", out_path)
