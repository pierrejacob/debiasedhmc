rm(list=ls())
library(ggplot2)
setmytheme()

# process output files
dimension <- 302
ntestfunct <- 2 * dimension
njobs <- 100
nreps <- 10
nrepeats <- njobs * nreps
estimates <- matrix(nrow = nrepeats, ncol = ntestfunct)
cost <- rep(0, nrepeats)
runtime <- rep(0, nrepeats)
for (ijob in 1:njobs){
  filename <- paste("inst/logistic/output/output.hmc.germancredit.repeat", ijob, ".RData", sep = "")
  load(file = filename)
  for (irep in 1:nreps){
    estimates[(ijob-1)*nreps+irep, ] <- unbiased_estimates[irep, ]
    cat((ijob-1)*nreps+irep, "\n")
    cost[(ijob-1)*nreps+irep] <- 2 * meetingtime[irep] + max(1, m + 1 - meetingtime[irep])
    runtime[(ijob-1)*nreps+irep] <- runtimes[irep]
  }
}

# load optimal HMC asymptotic variance
load("inst/logistic/output/optimal.hmc.germancredit.RData")

# compute relative inefficiency
mean_cost <- mean(cost)
variance <- sum(apply(estimates, 2, var))
inefficiency <- mean_cost * variance
cat("cost =", mean_cost, "\n")
cat("variance =", variance, "\n")
cat("inefficiency =", inefficiency, "\n")
cat("relative inefficiency =", inefficiency / optimal_hmc$asymptotic_variance, "\n")

