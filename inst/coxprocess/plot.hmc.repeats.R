rm(list=ls())
library(spatstat)
library(ggplot2)
library(gridExtra)
library(scales)
setmytheme()

# process output files
dimension <- 64^2
ntestfunct <- 2 * dimension
njobs <- 100
nreps <- 10
nrepeats <- njobs * nreps
estimates <- matrix(nrow = nrepeats, ncol = ntestfunct)
cost <- rep(0, nrepeats)
runtime <- rep(0, nrepeats)
for (ijob in 1:njobs){
  filename <- paste("inst/coxprocess/output/output.hmc.repeat", ijob, ".RData", sep = "")
  load(file = filename)
  for (irep in 1:nreps){
    estimates[(ijob-1)*nreps+irep, ] <- unbiased_estimates[irep, ]
    cat((ijob-1)*nreps+irep, "\n")
    cost[(ijob-1)*nreps+irep] <- 2 * meetingtime[irep] + max(1, m + 1 - meetingtime[irep])
    runtime[(ijob-1)*nreps+irep] <- runtimes[irep]
  }
}

# load optimal HMC asymptotic variance
load("inst/coxprocess/output/optimal.hmc.RData")

# compute relative inefficiency
mean_cost <- mean(cost)
variance <- sum(apply(estimates, 2, var))
inefficiency <- mean_cost * variance
cat("cost =", mean_cost, "\n")
cat("variance =", variance, "\n")
cat("inefficiency =", inefficiency, "\n")
cat("relative inefficiency =", inefficiency / optimal_hmc$asymptotic_variance, "\n")

# plot data, mean of latent process and its associated variance
ngrid <- 64
grid <- seq(from = 0, to = 1, length.out = ngrid+1)
data(finpines)
data_x <- (finpines$x + 5) / 10 # normalize data to unit square
data_y <- (finpines$y + 8) / 10
finpine.df <- data.frame(x = data_x, y = data_y)
g_data <- ggplot(finpine.df, aes(x = x, y = y)) + geom_point() +
          xlab("coordinate x") + ylab("coordinate y")

unbiased_mean <- colMeans(estimates)
variance_estimate <- unbiased_mean[(dimension+1):(2*dimension)] - unbiased_mean[1:dimension]^2
plot.unbiased.df <- data.frame()
for(i in 1:ngrid){
  for(j in 1:ngrid){
    plot.unbiased.df <- rbind(plot.unbiased.df, data.frame(x = grid[i],
                                                           y = grid[j],
                                                           latent = unbiased_mean[(i - 1) * ngrid + j],
                                                           var = variance_estimate[(i - 1) * ngrid + j]))
  }
}
g_mean <- ggplot(plot.unbiased.df, aes(x = x, y = y, fill = latent)) + geom_tile() +
          scale_fill_gradient(low = "white", high = "black") +
          xlab("coordinate x") + ylab("coordinate y")
g_var <- ggplot(plot.unbiased.df, aes(x = x, y = y, fill = var)) + geom_tile() +
         scale_fill_gradient(low = "white", high = "black") +
         xlab("coordinate x") + ylab("coordinate y")

grid.arrange(arrangeGrob(g_data, g_mean, g_var, nrow = 1, ncol = 3))
