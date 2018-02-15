### This file has been created by modifying files from the tmg package on CRAN, by Ari Pakman
rm(list=ls())
library(debiasedhmc)
library(tictoc)
library(coda)

# set.seed(1)
setmytheme()
par(mfrow = c(1,1))

dimension <- 2
# Define precision matrix and linear term
M <- matrix(c(1,0,0,1), dimension, dimension)
r <- c(4,4)

# Define two linear constraints
# first is 1.1 * x - 1 * y + 0 >= 0
f <- matrix(nrow = 1, ncol = dimension)
f[1,1] <- 1.1
f[1,2] <- -1
g <- 0

# second is -1 * x + 1 * y + 0 >= 0
f <- rbind(f, c(-1, 1))
g <- c(g, 0)

# Set initial positions for the HMC chains
initial1 <- c(2,2.1)
initial2 <- c(2,2.1)

# First run HMC
n <- 11000
burnin <- 1001
tic()
  samples <- rtmgR(n, M, r, initial1, f, g, q=NULL)
toc()
samples <- samples[burnin:11000,]
plot(samples, pch=".")

# Estimate asymptotic variance
asymptotic_variance <- sum(spectrum0.ar(samples)$spec) # variance of means
asymptotic_variance <- asymptotic_variance + sum(spectrum0.ar(samples^2)$spec) # variance of second moments
acf(samples[,1])
acf(samples[,2])
acf(samples[,1]^2)
acf(samples[,2]^2)
save(samples, asymptotic_variance, file = "inst/truncatednormal/output/hmc.linear.RData")

# Simulate meeting times
nreps <- 100
max_iterations <- 100
meetingtimes <- rep(0, nreps)
tic()
for (irep in 1:nreps){
  cchains <- rtmg_coupled_R(max_iterations, M, r, initial1, initial2, f, g, q = NULL)
  distances <- sapply(1:max_iterations, function(index) sqrt(sum(((cchains$samples1[index+1,] - cchains$samples2[index,]))^2)))
  meetingtimes[irep] <- which(distances == 0)[1]
}
toc()
summary(meetingtimes)

# Use guideline for choice of k and m
k <- floor(quantile(meetingtimes, 0.9))
m <- 10 * k

# Perform repeats
nreps <- 1000
meetingtimes <- rep(0, nreps)
unbiased_estimates <- list()

tic()
for (irep in 1:nreps){
  cchains <- rtmg_coupled_R(m, M, r, initial1, initial2, f, g, q = NULL)
  distances <- sapply(1:m, function(index) sqrt(sum(((cchains$samples1[index+1,] - cchains$samples2[index,]))^2)))
  meetingtime <- which(distances == 0)[1]
  meetingtimes[irep] <- meetingtime
  cchains$meetingtime <- meetingtime
  cchains$iteration <- m

  H_bar_estimates <- matrix(nrow = 2 * dimension, ncol = m)
  for(t in 1:m){
    H_bar_estimates[, t] <- H_bar(cchains, h = function(x) c(x, x^2), k = t, m = t)
  }
  unbiased_estimates[[irep]] <- H_bar_estimates
}
toc()
save(nreps, k, m, meetingtimes, unbiased_estimates, file = "inst/truncatednormal/output/repeats.linear.RData")
