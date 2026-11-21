#!/usr/bin/env Rscript
# What do the simulation generators ACTUALLY produce?  (no HAL needed, ~1 min)
#   Rscript scripts/dgp_distribution_check.R [n = 200000]
# For every simulation group, reports the empirical skewness, spread and
# correlation of the generated covariates (d = 5 generators; X1 also for d = 1),
# so Methods/Results can describe the realised distributions rather than the
# nominal `skewness_factor` / `correlation` parameters.
# Output: results/tables/dgp_distribution_check.csv
source("R/common.R")
source("R/simulation_config.R")
project_root()
args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args)) as.integer(args[[1L]]) else 200000L

skew <- function(x) { x <- x - mean(x); mean(x^3) / mean(x^2)^1.5 }   # = e1071 type 1
out <- list()
for (i in seq_len(nrow(simulation_groups))) {
  spec <- simulation_groups[i, , drop = FALSE]
  source(file.path("R", "generators", spec$generator))
  set.seed(seed_for(i, 999L, 1L))
  x <- get_simulation_dataset(spec, 5L, "smooth", n)
  u <- get_simulation_dataset(spec, 1L, "smooth", n)
  cm <- cor(x[, paste0("X", 1:5)])
  out[[i]] <- data.frame(
    group = spec$group, latent_corr = spec$correlation, skewness_factor = spec$skewness,
    skew_X1_d1 = skew(u$X1), sd_X1_d1 = sd(u$X1),
    skew_X1 = skew(x$X1), skew_X2 = skew(x$X2), skew_X4 = skew(x$X4), skew_X5 = skew(x$X5),
    sd_X1 = sd(x$X1), cor_X1_X2 = cm[1, 2], mean_abs_cor = mean(abs(cm[upper.tri(cm)])),
    X3_binary = length(unique(x$X3)) == 2L)
}
res <- do.call(rbind, out)
num <- vapply(res, is.numeric, logical(1))
res[num] <- lapply(res[num], round, 3)
atomic_write_csv(res, "results/tables/dgp_distribution_check.csv")
print(res, row.names = FALSE)
