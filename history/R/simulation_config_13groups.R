# One group per historical manuscript simulation setting.
simulation_groups <- data.frame(
  group = c("regular", "corr_0.1", "corr_0.5", "corr_0.9",
            "skewed", "skew_0.5_corr_0.1", "skew_0.5_corr_0.5",
            "skew_0.5_corr_0.9", "skew_2_corr_0.1", "skew_2_corr_0.5",
            "skew_2_corr_0.9", "skew_3_corr_0.1", "skew_5_corr_0.1"),
  generator = c("sim_data.R", rep("correlation.R", 3),
                "skewness_sim_nofactor.R", rep("skewness.R", 8)),
  correlation = c(NA, 0.1, 0.5, 0.9, NA, 0.1, 0.5, 0.9, 0.1, 0.5, 0.9, 0.1, 0.1),
  skewness = c(rep(NA, 5), rep(0.5, 3), rep(2, 3), 3, 5),
  test_n = c(rep(10000L, 4), 1000L, rep(10000L, 8)),
  stringsAsFactors = FALSE
)

# NOTE: skew=3 and skew=5 original source scripts are absent in supplied ZIP.
# Those two groups are reconstructed using the supplied skewness.R generator,
# exactly as the accessible skew=0.5/2 scripts pass the parameters.
get_simulation_dataset <- function(spec, d, shape, n) {
  prefix <- switch(as.character(d), `1` = "uni", `3` = "tri", `5` = "five")
  f <- get(sprintf("sim_%s_%s", prefix, shape), mode = "function")
  args <- list(n = n)
  if (!is.na(spec$correlation)) args$corr <- spec$correlation
  if (!is.na(spec$skewness)) args$skewness_factor <- spec$skewness
  do.call(f, args)
}
