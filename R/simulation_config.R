# Simulation groups. Groups 1-5 are unchanged from audited v2 (same order, so
# their seeds and any results already computed stay valid).
# Groups 6-17 (2026-11-16b): skewness x correlation factorial with the Gaussian
# copula / standardised-Gamma generator `skew_copula.R`, which controls the
# covariate skewness (0.5, 1, 2, 3) and the latent correlation (0.1, 0.5, 0.9)
# separately. They replace the normal-mixture groups of skewness.R, which
# produced |skewness| <= 0.33 and destroyed the correlation (see
# scripts/dgp_distribution_check.R). skewness.R is kept for reference only.
.skew_grid <- expand.grid(corr = c(0.1, 0.5, 0.9), skew = c(0.5, 1, 2, 3))
simulation_groups <- data.frame(
  group = c("regular", "corr_0.1", "corr_0.5", "corr_0.9", "skewed",
            sprintf("gskew_%s_corr_%s", .skew_grid$skew, .skew_grid$corr)),
  generator = c("sim_data.R", rep("correlation.R", 3), "skewness_sim_nofactor.R",
                rep("skew_copula.R", nrow(.skew_grid))),
  correlation = c(NA, 0.1, 0.5, 0.9, NA, .skew_grid$corr),
  skewness = c(rep(NA, 5), .skew_grid$skew),
  test_n = c(rep(10000L, 4), 1000L, rep(10000L, nrow(.skew_grid))),
  stringsAsFactors = FALSE
)

# For skew_copula.R, `skewness` is the target covariate skewness (passed as
# `skewness_factor` for compatibility with the generator signatures).
get_simulation_dataset <- function(spec, d, shape, n) {
  prefix <- switch(as.character(d), `1` = "uni", `3` = "tri", `5` = "five")
  f <- get(sprintf("sim_%s_%s", prefix, shape), mode = "function")
  args <- list(n = n)
  if (!is.na(spec$correlation)) args$corr <- spec$correlation
  if (!is.na(spec$skewness)) args$skewness_factor <- spec$skewness
  do.call(f, args)
}
