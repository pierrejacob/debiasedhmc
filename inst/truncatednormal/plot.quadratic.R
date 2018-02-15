rm(list=ls())
library(ggplot2)
setmytheme()

# plot truncated normal with quadratic inequalities
load("inst/truncatednormal/output/hmc.quadratic.RData")
plot(samples, pch=".")
nsamples <- 2000
gdist <- qplot(samples[1:nsamples,1], samples[1:nsamples,2], geom = "point", alpha = 0.25)
gdist <- gdist + xlab(expression(x[1])) + ylab(expression(x[2])) + theme(legend.position = "none")
gdist <- gdist + coord_fixed()
gdist

# plot histogram of meeting times
load("inst/truncatednormal/output/repeats.quadratic.RData")
ghist <- qplot(x = meetingtimes, geom = "blank") + geom_histogram(aes(y = ..density..), binwidth = 1.2)
ghist <- ghist + xlab("Meeting times") + ylab("Density")
ghist <- ghist
ghist

# compute relative inefficiency
ntestfunct <- 4
parameter_k <- k
parameter_m <- m
estimates <- matrix(nrow = nreps, ncol = ntestfunct)
cost <- rep(0, nreps)

for (irep in 1:nreps){
  if (parameter_k == parameter_m){
    estimates[irep, ] <- unbiased_estimates[[irep]][, parameter_k]
  } else {
    estimates[irep, ] <- rowMeans(unbiased_estimates[[irep]][, parameter_k:parameter_m])
  }
  cost[irep] <- 2 * meetingtimes[irep] + max(1, parameter_m + 1 - meetingtimes[irep])
}
mean_cost <- mean(cost)
variance <- sum(apply(estimates, 2, var))
inefficiency <- mean_cost * variance
cat("cost =", mean_cost, "\n")
cat("variance =", variance, "\n")
cat("inefficiency =", inefficiency, "\n")
cat("relative inefficiency =", inefficiency / asymptotic_variance, "\n")
