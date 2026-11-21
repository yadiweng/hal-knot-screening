#!/usr/bin/env Rscript
# Are 10 repeats enough?  Run AFTER scripts/summarize_results.R.
#   Rscript scripts/check_precision.R [target_r2] [target_rel_mse]
# target_r2      : acceptable 95% CI half-width for a paired R^2 difference (default 0.01)
# target_rel_mse : acceptable half-width for a paired real-data MSE difference,
#                  relative to the quantile MSE (default 0.02 = 2%)
# Distinguish (a) difference from zero, (b) interval entirely inside a
# PRACTICAL EQUIVALENCE margin, and (c) inconclusive. A narrow interval
# straddling zero is not sufficient to establish equivalence.
# Real-data repeats quantify split sensitivity, NOT population-level uncertainty.
source("R/common.R")
args <- commandArgs(trailingOnly = TRUE)
target_r2 <- if (length(args) >= 1L) as.numeric(args[[1L]]) else 0.01
target_rel <- if (length(args) >= 2L) as.numeric(args[[2L]]) else 0.02

repeats_needed <- function(sd, h) {
  if (!is.finite(sd) || sd == 0) return(NA_real_)
  R <- max(2, ceiling((qnorm(0.975) * sd / h)^2))
  while (qt(0.975, R - 1) * sd / sqrt(R) > h) R <- R + 1
  R
}

# ---- simulations: paired R^2 differences ----------------------------------------
f <- "results/tables/simulation_paired_differences.csv"
if (file.exists(f)) {
  pd <- read.csv(f, stringsAsFactors = FALSE)
  pd$half_width <- (pd$upper_95 - pd$lower_95) / 2
  pd$ci_excludes_0 <- pd$lower_95 > 0 | pd$upper_95 < 0
  pd$within_equivalence_margin <- pd$lower_95 >= -target_r2 & pd$upper_95 <= target_r2
  pd$classification <- ifelse(pd$ci_excludes_0, "direction_supported",
                              ifelse(pd$within_equivalence_margin,
                                     "within_practical_margin", "inconclusive"))
  pd$resolved <- pd$classification != "inconclusive"
  pd$repeats_for_target <- vapply(pd$sd, repeats_needed, numeric(1), h = target_r2)
  atomic_write_csv(pd, "results/tables/precision_simulation.csv")
  un <- pd[!pd$resolved, ]
  message(sprintf("Simulation: %d of %d paired comparisons resolved with %d repeats (target half-width %.3f).",
                  sum(pd$resolved), nrow(pd), max(pd$n_repeats), target_r2))
  if (nrow(un)) message(sprintf("  Unresolved: %d; repeats needed for target: median %g, max %g.",
                                nrow(un), median(un$repeats_for_target, na.rm = TRUE),
                                max(un$repeats_for_target, na.rm = TRUE)))
}

# ---- real data: paired CV-MSE differences across CV repeats ---------------------
f <- "results/tables/real_data_repeat_means.csv"
if (file.exists(f)) {
  rm <- read.csv(f, stringsAsFactors = FALSE)
  wide <- reshape(rm[c("dataset", "repeat_id", "knot_start", "method", "mse")],
                  idvar = c("dataset", "repeat_id", "knot_start"),
                  timevar = "method", direction = "wide")
  out <- do.call(rbind, lapply(c("kmeans", "pam"), function(m) {
    wide$rel_diff <- (wide[[paste0("mse.", m)]] - wide$mse.quantile) / wide$mse.quantile
    do.call(rbind, lapply(split(wide, list(wide$dataset, wide$knot_start), drop = TRUE), function(z) {
      n <- nrow(z); s <- sd(z$rel_diff); h <- qt(0.975, n - 1) * s / sqrt(n)
      data.frame(dataset = z$dataset[1], knot_start = z$knot_start[1],
                 comparison = paste0(m, "_vs_quantile"), n_repeats = n,
                 mean_rel_diff = mean(z$rel_diff), sd_rel_diff = s, half_width = h,
                 split_ci_excludes_zero = abs(mean(z$rel_diff)) > h,
                 split_interval_within_margin = abs(mean(z$rel_diff)) + h <= target_rel,
                 repeats_for_target = repeats_needed(s, target_rel))
    }))
  }))
  out$classification <- ifelse(out$split_ci_excludes_zero,
                               "split_direction_consistent",
                               ifelse(out$split_interval_within_margin,
                                      "within_split_margin", "inconclusive"))
  atomic_write_csv(out, "results/tables/precision_real_data.csv")
  message(sprintf("Real data: %d of %d split-sensitivity comparisons inconclusive (margin %.0f%%).",
                  sum(out$classification == "inconclusive"), nrow(out), 100 * target_rel))
  message("  NOTE: repeated-CV t intervals are NOT population 95% confidence intervals.")
}
