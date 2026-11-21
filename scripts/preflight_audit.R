#!/usr/bin/env Rscript
# Fast data/DGP audit. Run on a compute node before launching full arrays.
# Does not refit the manuscript's heavy simulations.
source("R/common.R")
source("R/simulation_config.R")
source("R/data_prep.R")
project_root()

crime <- load_dataset("communities_crime")
cat(sprintf("Communities after dropping IDs and complete cases: n=%d; predictors=%d\n",
            nrow(crime$X), ncol(crime$X)))
if (nrow(crime$X) != 319L) stop("Expected 319 Communities rows in bundled input; inspect the preprocessing")

# Verify that the formerly ignored d=1 skewness factor actually changes X.
source("R/generators/skewness.R")
set.seed(2468)
a <- sim_uni_smooth(500L, corr=0.1, skewness_factor=0.5)
set.seed(2468)
b <- sim_uni_smooth(500L, corr=0.1, skewness_factor=5)
if (isTRUE(all.equal(a$X1, b$X1))) stop("d=1 skewness_factor does not affect X1")
cat("Univariate skewness factor changes X1: OK\n")

# The jump error was a lack of noise; this check can detect that regression.
step_mean <- function(x) -2.7*(x < -3) + 2.5*(x > -2) - 2*(x > 0) +
  4*(x > 2) - 3*(x > 3)
for (s in c("correlation.R", "skewness.R")) {
  source(file.path("R", "generators", s))
  set.seed(1001)
  z <- if (s == "correlation.R") sim_uni_jump(3000, corr=0.1) else
    sim_uni_jump(3000, corr=0.1, skewness_factor=2)
  noise_sd <- sd(z$Y - step_mean(z$X1))
  cat(sprintf("%s: jump residual noise SD=%.3f\n", s, noise_sd))
  if (!is.finite(noise_sd) || abs(noise_sd - 1) > 0.1)
    stop("Jump DGP residual noise SD inconsistent with N(0,1)")
}

# Copula generator: realised skewness and rank correlation match their targets,
# and X3 stays binary.
source("R/generators/skew_copula.R")
skew_stat <- function(x) { x <- x - mean(x); mean(x^3) / mean(x^2)^1.5 }
for (g in c(0.5, 1, 2, 3)) for (r in c(0.1, 0.9)) {
  set.seed(77)
  z <- sim_five_smooth(50000L, corr = r, skewness_factor = g)
  sk <- skew_stat(z$X1); rs <- cor(z$X1, z$X2, method = "spearman")
  cat(sprintf("skew_copula target skew=%.1f corr=%.1f: skew(X1)=%.2f spearman=%.3f (theory %.3f)\n",
              g, r, sk, rs, 6 / pi * asin(r / 2)))
  if (abs(sk - g) > 0.15 * max(1, g) || abs(rs - 6 / pi * asin(r / 2)) > 0.02 ||
      length(unique(z$X3)) != 2L)
    stop("skew_copula.R does not reproduce its target skewness/correlation")
}

out <- list()
for (i in seq_len(nrow(simulation_groups))) {
  spec <- simulation_groups[i,,drop=FALSE]
  source(file.path("R", "generators", spec$generator))
  for (d in c(1L,3L,5L)) {
    for (shape in c("smooth", "jump", "sin")) {
      set.seed(seed_for(i, 1L, d*3L + match(shape, c("smooth", "jump", "sin"))))
      z <- get_simulation_dataset(spec,d,shape,250L)
      if (nrow(z)!=250L || any(!is.finite(as.matrix(z))))
        stop("Invalid generated data: ", spec$group, " d=",d," shape=",shape)
      x <- z[,setdiff(names(z),"Y"),drop=FALSE]
      observed_abs_corr <- if (ncol(x)>1) {
        m <- suppressWarnings(cor(x))
        mean(abs(m[upper.tri(m)]),na.rm=TRUE)
      } else NA_real_
      out[[length(out)+1L]] <- data.frame(
        group=spec$group,d=d,scenario=shape,
        sd_x1=sd(z$X1), mean_x1=mean(z$X1),
        mean_absolute_x_correlation=observed_abs_corr)
    }
  }
}
atomic_write_csv(do.call(rbind,out),"results/tables/dgp_diagnostics.csv")
cat("All", nrow(simulation_groups), "x 9 generator probes returned finite data. Wrote results/tables/dgp_diagnostics.csv\n")

if (requireNamespace("hal9001",quietly=TRUE)) {
  # Small synthetic all-binary predictor case: does the custom fork accept it?
  set.seed(2026)
  x <- data.frame(a=rbinom(80,1,0.5),b=rbinom(80,1,0.5))
  y <- 0.4*x$a - 0.2*x$b + rnorm(80)
  q <- fit_hal_screened(x[1:60,],y[1:60],x[61:80,],"quantile",2L,c(5,3))
  kb <- clustering_budget(q$basis_list,0L,2L)
  cat("Binary-only knot budget: ",paste(kb,collapse=","),"\n")
  for (method in c("kmeans","pam")) {
    result <- try(fit_hal_screened(x[1:60,],y[1:60],x[61:80,],method,2L,kb),silent=TRUE)
    if (inherits(result,"try-error")) {
      cat("Binary-only ",method," FAILED: ",as.character(result),"\n",sep="")
      stop("Custom HAL fork binary-only case needs a scientifically defined fallback")
    }
    cat("Binary-only ",method," fit succeeded; basis count ",result$total_bases,"\n",sep="")
  }
} else warning("Skipping binary-only HAL smoke test; hal9001 unavailable")
cat("Preflight complete. Still inspect actual full-job logs and output summaries.\n")
