rm(list=ls())
library(debiasedhmc)
library(spatstat)
library(tictoc)

# load pine saplings dataset
data(finpines)
data_x <- (finpines$x + 5) / 10 # normalize data to unit square
data_y <- (finpines$y + 8) / 10
plot(x = data_x, y = data_y, type = "p")

ngrid <- 64
grid <- seq(from = 0, to = 1, length.out = ngrid+1)
dimension <- ngrid^2
data_counts <- rep(0, dimension)
for (i in 1:ngrid){
  for (j in 1:ngrid){
    logical_y <- (data_x > grid[i]) * (data_x < grid[i+1])
    logical_x <- (data_y > grid[j]) * (data_y < grid[j+1])
    data_counts[(i-1)*ngrid + j] <- sum(logical_y * logical_x)
  }
}

# prior distribution
parameter_sigmasq <- 1.91
parameter_mu <- log(126) - 0.5 * parameter_sigmasq
parameter_beta <- 1 / 33
parameter_area <- 1 / dimension

prior_mean <- rep(parameter_mu, dimension)
prior_cov <- matrix(nrow = dimension, ncol = dimension)
for (m in 1:dimension){
  for (n in 1:dimension){
    index_m <- c( floor((m-1) / ngrid) + 1, ((m-1) %% ngrid) + 1 )
    index_n <- c( floor((n-1) / ngrid) + 1, ((n-1) %% ngrid) + 1 )
    prior_cov[m,n] <- parameter_sigmasq * exp(- sqrt(sum((index_m - index_n)^2)) / (ngrid * parameter_beta) )
  }
}
prior_precision <- solve(prior_cov)
prior_precision_chol <- t(chol(prior_precision))
prior <- list()
prior$logdensity <- function(x){
  return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), prior_mean, prior_precision_chol))
}
prior$gradlogdensity <- function(x){
  return(gradlognormal(x, prior_mean, prior_precision))
}

# likelihood function
likelihood <- list()
likelihood$log <- function(x){
  return(coxprocess_loglikelihood(matrix(x, nrow=1), data_counts, parameter_area))
}
likelihood$gradlog <- function(x){
  return( data_counts - parameter_area * exp(x) )
}

# posterior distribution
logtarget <- function(x) prior$logdensity(x) + likelihood$log(x)
gradlogtarget <- function(x) prior$gradlogdensity(x) + likelihood$gradlog(x)

# initial distribution
rinit <- function() as.numeric(fast_rmvnorm(1, prior_mean, prior_cov))

save.image(file = "inst/coxprocess/coxprocess.RData")

# compute metric tensor as in Girolami and Calderhead 2011 (Section 9)
# takes a few minutes to compute
tic()
metric <- list()
metric$tensor <- prior_precision
diag(metric$tensor) <- parameter_area * exp(prior_mean + 0.5 * diag(prior_cov)) + diag(prior_precision)
metric$inverse <- solve(metric$tensor)
cat("Inverse done")
metric$chol_inverse <- t(chol(metric$inverse))
cat("Chol done")
metric$inverse_chol_inverse <- solve(t(metric$chol_inverse))
toc()

save.image(file = "inst/coxprocess/coxprocess_with_metric.RData")

