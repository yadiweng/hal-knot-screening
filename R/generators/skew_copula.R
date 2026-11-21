# Skewed + correlated covariates with SEPARATELY CONTROLLED skewness and
# correlation (replacement for the normal-mixture design in skewness.R).
#
# Same structure as correlation.R (Tables 2-4):
#   S ~ N_4(0, Sigma), Sigma = exchangeable with latent correlation `corr`
#   X3 ~ Bernoulli(plogis(0.25 X1 - 0.25 X2)),  X4 = Z3 + X1 * X3
# Only the marginal transform changes: instead of Uniform(-4, 4) margins, every
# continuous covariate is a standardised Gamma quantile of its latent normal,
#   Z_j = (qgamma(pnorm(S_j), shape = k) - k) / sqrt(k),   k = 4 / skewness^2,
# so Z_j has mean 0, SD 1 and population skewness exactly `skewness`
# (argument `skewness_factor`, kept for compatibility with the pipeline).
#
#   X1, X2 = 8/sqrt(12) * Z     same mean (0) and SD (2.31) as Uniform(-4, 4)
#   Z3     = Z                  same mean/SD as the N(0, 1) used for X4
#   X5     = 4 + sqrt(2) * Z    same mean/SD as Gamma(2, 1) + 2; skewness = sqrt(2)
#                               reproduces the original X5 exactly
#
# Consequences:
# * Location and scale match the correlated design, so changing `skewness`
#   changes the SHAPE of the covariates, not their range (the mixture design
#   widened X from SD 2.3 to 12.7, which drove most R^2 changes).
# * All margins are monotone transforms of S, so the rank (Spearman)
#   correlation of any two continuous covariates is exactly
#   (6/pi) * asin(corr / 2) = 0.096, 0.483, 0.891 for corr = 0.1, 0.5, 0.9,
#   whatever the skewness. Pearson correlation is slightly lower and is
#   reported by scripts/dgp_distribution_check.R.
# * X3 stays binary.
# Regression functions (Y) are identical to correlation.R / sim_data.R.
library(MASS)

skewed_margin <- function(u, skewness) {
  if (!is.finite(skewness) || skewness <= 0) stop("skewness must be > 0")
  k <- 4 / skewness^2
  (qgamma(u, shape = k, rate = 1) - k) / sqrt(k)
}

copula_covariates <- function(n, d, corr, skewness) {
  Sigma <- matrix(corr, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  U <- pnorm(mvrnorm(n, mu = rep(0, 4), Sigma = Sigma))
  sd_unif <- 8 / sqrt(12)
  X1 <- sd_unif * skewed_margin(U[, 1], skewness)
  if (d == 1) return(data.frame(X1 = X1))
  X2 <- sd_unif * skewed_margin(U[, 2], skewness)
  X3 <- rbinom(n, 1, plogis(0.25 * X1 - 0.25 * X2))
  if (d == 3) return(data.frame(X1 = X1, X2 = X2, X3 = X3))
  X4 <- skewed_margin(U[, 3], skewness) + X1 * X3
  X5 <- 4 + sqrt(2) * skewed_margin(U[, 4], skewness)
  data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5)
}

# univariate -------------------------------------------------------------------
sim_uni_smooth <- function(n, corr, skewness_factor) {
  d <- copula_covariates(n, 1, corr, skewness_factor); X1 <- d$X1
  d$Y <- 0.05*X1 - 0.42*X1^2 + rnorm(n)
  d
}
sim_uni_jump <- function(n, corr, skewness_factor) {
  d <- copula_covariates(n, 1, corr, skewness_factor); X1 <- d$X1
  d$Y <- -2.7*as.numeric(X1 < -3) + 2.5*as.numeric(X1 > -2) - 2*as.numeric(X1 > 0) +
    4*as.numeric(X1 > 2) - 3*as.numeric(X1 > 3) + rnorm(n)
  d
}
sim_uni_sin <- function(n, corr, skewness_factor) {
  d <- copula_covariates(n, 1, corr, skewness_factor); X1 <- d$X1
  d$Y <- 2*sin(0.5*pi*abs(X1)) + 2*cos(0.5*pi*abs(X1)) + rnorm(n)
  d
}

# trivariate -------------------------------------------------------------------
sim_tri_smooth <- function(n, corr, skewness_factor) {
  d <- copula_covariates(n, 3, corr, skewness_factor)
  X1 <- d$X1; X2 <- d$X2; X3 <- d$X3
  d$Y <- X1/15 - 0.28*X1^2 + 0.5*X2 + 0.25*X2*X3 + rnorm(n)
  d
}
sim_tri_jump <- function(n, corr, skewness_factor) {
  d <- copula_covariates(n, 3, corr, skewness_factor)
  X1 <- d$X1; X2 <- d$X2; X3 <- d$X3
  d$Y <- -2*as.numeric(X1 < -3)*X3 + 2.5*as.numeric(X1 > -2) -
    2*as.numeric(X1 > 0) + 2.5*as.numeric(X1 > 2)*X3 - 2.5*as.numeric(X1 > 3) +
    as.numeric(X2 > -1) - 4*as.numeric(X2 > 1)*X3 + 2*as.numeric(X2 > 3) + rnorm(n)
  d
}
sim_tri_sin <- function(n, corr, skewness_factor) {
  d <- copula_covariates(n, 3, corr, skewness_factor)
  X1 <- d$X1; X2 <- d$X2; X3 <- d$X3
  d$Y <- X3 * (4 * sin(pi / 2 * abs(X1)) * (X2 < 0) +
               4.1 * cos(pi / 2 * abs(X1)) * (X2 > 0)) + rnorm(n)
  d
}

# five-variate -----------------------------------------------------------------
sim_five_smooth <- function(n, corr, skewness_factor) {
  d <- copula_covariates(n, 5, corr, skewness_factor)
  X1 <- d$X1; X2 <- d$X2; X3 <- d$X3; X4 <- d$X4; X5 <- d$X5
  d$Y <- X1/10 - 0.3*X1^2 + 0.25*X2 + 0.5*X2*X3 - 0.5*X4 + 0.04*X5^2 - 0.1*X5 + rnorm(n)
  d
}
sim_five_jump <- function(n, corr, skewness_factor) {
  d <- copula_covariates(n, 5, corr, skewness_factor)
  X1 <- d$X1; X2 <- d$X2; X3 <- d$X3; X4 <- d$X4; X5 <- d$X5
  d$Y <- -as.numeric(X1 < -3)*X3 + 0.5*as.numeric(X1 > -2) - as.numeric(X1 > 0) +
    2*as.numeric(X1 > 2)*X3 - 3*as.numeric(X1 > 3) + 1.5*as.numeric(X2 > -1) -
    5*as.numeric(X2 > 1)*X3 + 2*as.numeric(X2 > 3) + 2*as.numeric(X4 < 0) -
    as.numeric(X5 > 5) - as.numeric(X4 < 0)*as.numeric(X1 < 0) + 2*X3 + rnorm(n)
  d
}
sim_five_sin <- function(n, corr, skewness_factor) {
  d <- copula_covariates(n, 5, corr, skewness_factor)
  X1 <- d$X1; X2 <- d$X2; X3 <- d$X3; X4 <- d$X4; X5 <- d$X5
  d$Y <- X3*(3.75*sin(pi/2*abs(X1))*(X2 < 0) + 4*cos(pi/2*abs(X1))*(X2 > 0)) +
    sin(X4*pi)*X5/10 + cos(abs(X4-X5))*X3 + rnorm(n)
  d
}
