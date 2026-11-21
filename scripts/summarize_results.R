#!/usr/bin/env Rscript
# Creates tidy mean/SD/CI tables from COMPLETE replicate-level CSV files.
source("R/common.R")
source("R/simulation_config.R")
root <- project_root()
args <- commandArgs(trailingOnly = TRUE)
expected <- if (length(args)) as.integer(args[[1L]]) else 10L
if (is.na(expected) || expected < 2L) stop("Expected number of repeats must be >=2")

group_ci <- function(df, group_keys, metric, interval = TRUE) {
  if (nrow(df) == 0L) return(data.frame())
  groups <- interaction(df[, group_keys, drop = FALSE], drop = TRUE, lex.order = TRUE)
  pieces <- split(df, groups)
  result <- lapply(pieces, function(z) {
    values <- z[[metric]]
    n <- length(values)
    est <- mean(values)
    sdv <- if (n > 1L) sd(values) else NA_real_
    lower <- if (interval && n > 1L) est - qt(.975, n - 1) * sdv / sqrt(n) else NA_real_
    upper <- if (interval && n > 1L) est + qt(.975, n - 1) * sdv / sqrt(n) else NA_real_
    cbind(z[1L, group_keys, drop = FALSE],
          data.frame(metric = metric, n_repeats = n, mean = est, sd = sdv,
                     lower_95 = lower, upper_95 = upper,
                     stringsAsFactors = FALSE))
  })
  do.call(rbind, result)
}

# SIMULATIONS: one row per group x repeat x DGP x method
simfiles <- list.files("results/simulations", "_repeat_[0-9]+\\.csv$", full.names = TRUE)
if (length(simfiles)) {
  sim <- do.call(rbind, lapply(simfiles, read.csv, stringsAsFactors = FALSE))
  assert_checkpoint_version(sim, "loaded result files")
  basekeys <- c("group", "dimension", "scenario", "method")
  if (!setequal(unique(sim$group), simulation_groups$group)) {
    stop("Missing simulation GROUPS: ", paste(setdiff(simulation_groups$group, unique(sim$group)), collapse=", "))
  }
  if (nrow(sim) != nrow(simulation_groups) * expected * 9L * 3L) {
    stop("Expected ", nrow(simulation_groups) * expected * 9L * 3L, " complete simulation rows; found ", nrow(sim))
  }
  duplicates <- duplicated(sim[c("group", "repeat_id", "dimension", "scenario", "method")])
  if (any(duplicates)) stop("Duplicate simulation rows; inspect checkpoint files")
  # Check completeness before reporting any confidence intervals.
  tab <- aggregate(sim$repeat_id, sim[basekeys], function(x) length(unique(x)))
  if (any(tab$x != expected)) stop("Missing simulation repeats; inspect results/simulations before summarizing")
  summary_sim <- do.call(rbind, lapply(c("r2", "total_bases", "active_bases"), function(m) {
    group_ci(sim, basekeys, m, interval = TRUE)
  }))
  atomic_write_csv(summary_sim, "results/tables/simulation_summary.csv")
  # Paired differences (same generated train/test data for each algorithm).
  wide <- reshape(sim[c("group", "repeat_id", "dimension", "scenario", "method", "r2")],
                  idvar = c("group", "repeat_id", "dimension", "scenario"),
                  timevar = "method", direction = "wide")
  wide$kmeans_minus_quantile <- wide$r2.kmeans - wide$r2.quantile
  wide$pam_minus_quantile <- wide$r2.pam - wide$r2.quantile
  paired <- do.call(rbind, lapply(c("kmeans_minus_quantile", "pam_minus_quantile"), function(m) {
    group_ci(wide, c("group", "dimension", "scenario"), m, interval = TRUE)
  }))
  atomic_write_csv(paired, "results/tables/simulation_paired_differences.csv")
  message("Simulation summaries written")
} else message("No simulation data found")

# REAL DATA: mean each 10-fold CV first; aggregate SD across 10 CV repeats.
realfiles <- list.files("results/real_data", "_repeat_[0-9]+\\.csv$", full.names = TRUE)
if (length(realfiles)) {
  raw <- do.call(rbind, lapply(realfiles, read.csv, stringsAsFactors = FALSE))
  assert_checkpoint_version(raw, "loaded result files")
  expected_datasets <- c("cpu", "automobile", "communities_crime", "pima", "oecdpanel",
                         "bike_sharing", "liver_disorders", "laheart", "real_estate")
  if (!setequal(unique(raw$dataset), expected_datasets)) {
    stop("Missing real DATASETS: ", paste(setdiff(expected_datasets, unique(raw$dataset)), collapse=", "))
  }
  if (nrow(raw) != 9L * expected * 9L * 10L * 3L) {
    stop("Expected ", 9L * expected * 9L * 10L * 3L, " complete real-data fold rows; found ", nrow(raw))
  }
  foldkeys <- c("dataset", "repeat_id", "knot_start", "method")
  dup <- duplicated(raw[c(foldkeys, "fold")])
  if (any(dup)) stop("Duplicate real-data fold records")
  foldcount <- aggregate(raw$fold, raw[foldkeys], function(x) length(unique(x)))
  if (any(foldcount$x != 10L)) stop("Incomplete 10-fold CV repeat; do not report partial averages")
  # Equal-weight folds, matching the historical average of 10 fold MSE values.
  avg <- aggregate(raw[c("mse", "total_bases", "active_bases")], raw[foldkeys], mean)
  key <- c("dataset", "knot_start", "method")
  count <- aggregate(avg$repeat_id, avg[key], function(x) length(unique(x)))
  if (any(count$x != expected)) stop("Missing real-data repeats; inspect results/real_data")
  atomic_write_csv(avg, "results/tables/real_data_repeat_means.csv")
  real_summary <- do.call(rbind, lapply(c("mse", "total_bases", "active_bases"), function(m) {
    group_ci(avg, key, m, interval = FALSE)
  }))
  # SD quantifies variability across CV splits of ONE observed dataset;
  # it is NOT an iid sampling confidence interval for population generalization.
  atomic_write_csv(real_summary, "results/tables/real_data_summary.csv")
  message("Real data summaries written (mean and SD; no false 95% CI)")
} else message("No real data found")
