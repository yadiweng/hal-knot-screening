#!/usr/bin/env Rscript
# Optimal ("oracle") test R^2 of every simulation setting.
#   Rscript scripts/optimal_r2.R [N = 2e6]
# With Y = f0(X) + e, e ~ N(0, 1), the population R^2 of the true regression
# function is Var f0(X) / (Var f0(X) + 1) = 1 - 1 / Var(Y). It is estimated by
# Monte Carlo with N draws per setting, using the same generators as the
# simulations. No HAL fitting is needed (about 1-2 minutes).
# Output: results/tables/optimal_r2.csv  (quoted in Section 5 and Appendix A)
source("R/common.R")
source("R/simulation_config.R")
invisible(project_root())
args <- commandArgs(trailingOnly = TRUE)
N <- if (length(args)) as.numeric(args[[1L]]) else 2e6
set.seed(1)
out <- list()
for (gi in seq_len(nrow(simulation_groups))) {
  spec <- simulation_groups[gi, , drop = FALSE]
  source(file.path("R", "generators", spec$generator))
  for (d in c(1L, 3L, 5L)) for (sh in c("smooth", "jump", "sin")) {
    y <- get_simulation_dataset(spec, d, sh, N)$Y
    out[[length(out) + 1L]] <- data.frame(group = spec$group, dimension = d,
                                          scenario = sh, optimal_r2 = 1 - 1 / var(y))
  }
}
res <- do.call(rbind, out)
atomic_write_csv(res, "results/tables/optimal_r2.csv")
rng <- aggregate(optimal_r2 ~ group, res, function(x) round(range(x), 3))
print(rng, row.names = FALSE)
