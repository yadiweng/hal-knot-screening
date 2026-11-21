library(MASS)
library(rlang)

# univariate -------------------------------------------------------------------
# smooth
sim_uni_smooth <- function(n, corr) {
  mu <- c(0, 0, 0, 0)
  Sigma <- matrix(corr, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  samples <- mvrnorm(n, mu = mu, Sigma = Sigma)
  X1 <- 8 * (pnorm(samples[,1]) - 0.5)
  Y <- 0.05*X1 - 0.42*X1^2 + rnorm(n)
  return(data.frame(X1 = X1, Y = Y))
}

#jumps
sim_uni_jump <- function(n, corr){
  mu <- c(0, 0, 0, 0)
  Sigma <- matrix(corr, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  samples <- mvrnorm(n, mu = mu, Sigma = Sigma)
  X1 <- 8 * (pnorm(samples[,1]) - 0.5)
  # FIXED (2026-11): the original line break made R drop the last two jumps
  # and the N(0,1) noise; now identical to sim_uni_jump() in sim_data.R.
  Y <- -2.7*as.numeric(X1 < -3) + 2.5*as.numeric(X1 > -2) - 2*as.numeric(X1 > 0) +
    4*as.numeric(X1 > 2) - 3*as.numeric(X1 > 3) + rnorm(n)
  return(data.frame(X1 = X1, Y = Y))
}

#sinusoidal
sim_uni_sin <- function(n, corr){
  mu <- c(0, 0, 0, 0)
  Sigma <- matrix(corr, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  samples <- mvrnorm(n, mu = mu, Sigma = Sigma)
  X1 <- 8 * (pnorm(samples[,1]) - 0.5)
  Y <- 2*sin(0.5*pi*abs(X1)) + 2*cos(0.5*pi*abs(X1)) + rnorm(n)
  return(data.frame(X1 = X1, Y = Y))
}

# trivariate -------------------------------------------------------------------
# smooth
sim_tri_smooth <- function(n, corr) {
  mu <- c(0, 0, 0, 0)
  Sigma <- matrix(corr, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  samples <- mvrnorm(n, mu = mu, Sigma = Sigma)
  X1 <- 8 * (pnorm(samples[,1]) - 0.5)
  X2 <- 8 * (pnorm(samples[,2]) - 0.5)
  X3 <- rbinom(n, 1, plogis(0.25 * X1 - 0.25 * X2))
  Y <- X1/15-0.28*X1^2+0.5*X2+0.25*X2*X3+rnorm(n)
  data <- data.frame(X1 = X1, X2 = X2, X3 = X3, Y = Y)
  return(data)
}


# jumps
sim_tri_jump <- function(n, corr){
  mu <- c(0, 0, 0, 0)
  Sigma <- matrix(corr, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  samples <- mvrnorm(n, mu = mu, Sigma = Sigma)
  X1 <- 8 * (pnorm(samples[,1]) - 0.5)
  X2 <- 8 * (pnorm(samples[,2]) - 0.5)
  X3 <- rbinom(n, 1, plogis(0.25 * X1 - 0.25 * X2))
  Y <- -2*as.numeric(X1 < -3)*X3+2.5*as.numeric(X1 > -2)-
    2*as.numeric(X1 > 0)+2.5*as.numeric(X1 > 2)*X3-2.5*as.numeric(X1 > 3) +
    as.numeric(X2 > -1)-4*as.numeric(X2 > 1)*X3+2*as.numeric(X2 > 3)+rnorm(n)
  data <- data.frame(X1 = X1, X2 = X2, X3 = X3, Y = Y)
  return(data)
}

# sinusoidal
sim_tri_sin <- function(n, corr){
  mu <- c(0, 0, 0, 0)
  Sigma <- matrix(corr, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  samples <- mvrnorm(n, mu = mu, Sigma = Sigma)
  X1 <- 8 * (pnorm(samples[,1]) - 0.5)
  X2 <- 8 * (pnorm(samples[,2]) - 0.5)
  X3 <- rbinom(n, 1, plogis(0.25 * X1 - 0.25 * X2))
  Y <- X3 * (4 * sin(pi / 2 * abs(X1)) * (X2 < 0) + 4.1 * cos(pi / 2 * abs(X1)) * (X2 > 0)) + rnorm(n)
  data <- data.frame(X1 = X1, X2 = X2, X3 = X3, Y = Y)
  return(data)
}

# five-variate -----------------------------------------------------------------
# smooth
sim_five_smooth <- function(n, corr) {
  mu <- c(0, 0, 0, 0)
  Sigma <- matrix(corr, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  samples <- mvrnorm(n, mu = mu, Sigma = Sigma)
  X1 <- 8 * (pnorm(samples[,1]) - 0.5)
  X2 <- 8 * (pnorm(samples[,2]) - 0.5)
  X3 <- rbinom(n, 1, plogis(0.25 * X1 - 0.25 * X2))
  X4 <- samples[,3] + X1 * X3
  X5 <- qgamma(pnorm(samples[,4]), shape = 2, rate = 1) + 2
  Y <- X1/10 - 0.3*X1^2 + 0.25*X2 + 0.5*X2*X3 - 0.5*X4 + 0.04*X5^2 - 0.1*X5 + rnorm(n)
  data <- data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, Y = Y)
  return(data)
}

# jumps
sim_five_jump <- function(n, corr) {
  mu <- c(0, 0, 0, 0)
  Sigma <- matrix(corr, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  samples <- mvrnorm(n, mu = mu, Sigma = Sigma)
  X1 <- 8 * (pnorm(samples[,1]) - 0.5)
  X2 <- 8 * (pnorm(samples[,2]) - 0.5)
  X3 <- rbinom(n, 1, plogis(0.25 * X1 - 0.25 * X2))
  X4 <- samples[,3] + X1 * X3
  X5 <- qgamma(pnorm(samples[,4]), shape = 2, rate = 1) + 2
  Y <- -as.numeric(X1 < -3)*X3 + 0.5*as.numeric(X1 > -2) - as.numeric(X1 > 0) +
    2*as.numeric(X1 > 2)*X3 - 3*as.numeric(X1 > 3) + 1.5*as.numeric(X2 > -1) -
    5*as.numeric(X2 > 1)*X3 + 2*as.numeric(X2 > 3) + 2*as.numeric(X4 < 0) -
    as.numeric(X5 > 5) - as.numeric(X4 < 0)*as.numeric(X1 < 0) + 2*X3 + rnorm(n)
  data <- data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, Y = Y)
  return(data)
}

# sinusoidal
sim_five_sin <- function(n, corr) {
  mu <- c(0, 0, 0, 0)
  Sigma <- matrix(corr, nrow = 4, ncol = 4)
  diag(Sigma) <- 1
  samples <- mvrnorm(n, mu = mu, Sigma = Sigma)
  X1 <- 8 * (pnorm(samples[,1]) - 0.5)
  X2 <- 8 * (pnorm(samples[,2]) - 0.5)
  X3 <- rbinom(n, 1, plogis(0.25 * X1 - 0.25 * X2))
  X4 <- samples[,3] + X1 * X3
  X5 <- qgamma(pnorm(samples[,4]), shape = 2, rate = 1) + 2
  Y <- X3*(3.75*sin(pi/2*abs(X1))*(X2 < 0) + 4*cos(pi/2*abs(X1))*(X2 > 0)) +
    sin(X4*pi)*X5/10 + cos(abs(X4-X5))*X3 + rnorm(n)
  data <- data.frame(X1 = X1, X2 = X2, X3 = X3, X4 = X4, X5 = X5, Y = Y)
  return(data)
}
