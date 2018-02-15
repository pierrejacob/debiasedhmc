### This file has been created by modifying files from the tmg package on CRAN, by Ari Pakman
rm(list=ls())
library(debiasedhmc)
library(ggplot2)
library(tictoc)
library(coda)

set.seed(18)
setmytheme()
par(mfrow = c(1,1))

dimension <- 2
# Define precision matrix and linear term
M <- matrix(c(1,0,0,1), 2, 2)
r <- c(0,0)

# Define two quadratic constraints
A1 <- matrix(c(-1/32,0,0,-1/8), 2, 2)
B1 <- c(.25,.25)
C1 <- 3/8
constr1 <- list(A1,B1,C1)

A2 <- matrix(c(4,-1,-1,8), 2, 2)
B2 <- c(0,5)
C2 <- -1
constr2 <- list(A2,B2,C2)
q <- list(constr1,constr2)

# Set initial point for the Markov chain
initial1 <- c(2,0)
initial2 <- c(2,0)

# First run HMC
n <- 11000
burnin <- 1001
tic()
  samples <- rtmgR(n, M, r, initial1, f=NULL, g=NULL, q=q);
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
save(samples, asymptotic_variance, file = "inst/truncatednormal/output/hmc.quadratic.RData")

# Simulate meeting times
nreps <- 100
max_iterations <- 100
meetingtimes <- rep(0, nreps)
tic()
for (irep in 1:nreps){
  cchains <- rtmg_coupled_R(max_iterations, M, r, initial1, initial2, f=NULL, g=NULL, q=q);
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
  cchains <- rtmg_coupled_R(m, M, r, initial1, initial2, f=NULL, g=NULL, q=q);
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
save(nreps, k, m, meetingtimes, unbiased_estimates, file = "inst/truncatednormal/output/repeats.quadratic.RData")
