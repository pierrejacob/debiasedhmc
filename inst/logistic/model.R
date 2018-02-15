rm(list=ls())
library(debiasedhmc)

# simulate dataset
# nsamples <- 1000
# dimension <- 2^5+2
# true_alpha <- 0
# true_beta <- rep(1, dimension-2)
# logistic_function <- function(x) 1 / (1 + exp(-x))
# design_matrix <- matrix(nrow = nsamples, ncol = dimension-2)
# response <- rep(0, nsamples)
# for (isample in 1:nsamples){
#   # sample from Rademacher distribution
#   covariate <- sample(c(-1,1), size = dimension-2, replace = TRUE, prob = rep(0.5, 2))
#   # normalize
#   covariate <- covariate / sqrt(sum(covariate^2))
#   design_matrix[isample, ] <- covariate
#   response[isample] <- ( runif(1) < logistic_function(true_alpha + sum(true_beta * covariate)) )
# }

# load german credit dataset
germancredit <- read.table(system.file("", "germancredit.txt", package = "debiasedhmc"))
design_matrix <- scale(germancredit[, 1:24])
response <- germancredit[, 25] - 1
nsamples <- nrow(design_matrix)
dimension <- ncol(design_matrix)
interaction_terms <- matrix(nrow = nsamples, ncol = dimension*(dimension-1) / 2)
index <- 1
for (j in 1:(dimension-1)){
  for (jprime in (j+1):dimension){
    interaction_terms[, index] <- design_matrix[, j] * design_matrix[, jprime]
    index <- index + 1
  }
}
design_matrix <- cbind(design_matrix, scale(interaction_terms))
colnames(design_matrix) <- NULL
dimension <- ncol(design_matrix) + 2

# posterior distribution
lambda <- 0.01 # exponential rate parameter
logtarget <- function(x) debiasedhmc:::logistic_logtarget_c(x, response, design_matrix, lambda)
gradlogtarget <- function(x) debiasedhmc:::logistic_gradlogtarget_c(x, response, design_matrix, lambda)

# initial distribution
rinit <- function() rnorm(dimension)

filename <- "inst/logistic/germancredit.RData"
save(nsamples, dimension, response, design_matrix, lambda, logtarget, gradlogtarget, rinit,
     file = filename)
