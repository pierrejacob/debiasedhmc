library(debiasedhmc)
library(coda)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
set.seed(18)
setmytheme()

dimension <- 250
mean_target <- rep(0, dimension)
Sigma_target <- diag(1, dimension, dimension)
for (i in 1:dimension){
  for (j in 1:dimension){
    Sigma_target[i,j] <- exp(-abs(i-j))
  }
}
# trajectories lengths
originaleffectiveTimes <- c(pi/4, 2*pi/4, 3*pi/4, pi, 5*pi/4, 6*pi/4)
#
load(file = "mvnorm.hmctuning.RData")

library(latex2exp)
g <- ggplot(fixednsteps.mcmc.df, aes(x = effectiveTime, y = mcmcvar)) + geom_point() + scale_y_log10()
g <- g + scale_x_continuous(breaks = effectiveTimes, labels = c(TeX("$\\pi/4$"), TeX("$2 \\pi/4$"),
                                                                TeX("$3\\pi/4$"), TeX("$4 \\pi/4$"),
                                                                TeX("$5\\pi/4$"), TeX("$6\\pi/4$")))
g <- g + ylab("HCMC variance") + xlab("trajectory length")
g
ggsave(filename = "mvnorm.hmcvariance.trajectorylength.pdf", plot = g, width = 9, height = 6)

# from this we can identify the best trajectory length, i.e. the one that leads to highest efficiency of HMC

## Then we look at contraction of two HMC trajectories
## again as a function of the trajectory length
moreeffectiveTimes <- c(pi/4, 3*pi/8, 7*pi/16, 2*pi/4, 9*pi/16, 5*pi/8, 3*pi/4, pi, 5*pi/4, 6*pi/4)
load(file = "mvnorm.contraction.RData")
distance.df %>% head

g <- ggplot(distance.df, aes(x = effectiveTime, y = distance)) + geom_point()
g <- g + scale_y_log10(breaks = c(1e-20, 1e-15, 1e-10, 1e-5, 1), limits = c(1e-20, 10))
g <- g + scale_x_continuous(breaks = originaleffectiveTimes, labels = c(TeX("$\\pi/4$"), TeX("$2 \\pi/4$"),
                                                                TeX("$3\\pi/4$"), TeX("$4 \\pi/4$"),
                                                                TeX("$5\\pi/4$"), TeX("$6\\pi/4$")))

g <- g + ylab("distance after 100 iterations") + xlab("trajectory length")
g
ggsave(filename = "mvnorm.distance100iter.trajectorylength.pdf", plot = g, width = 9, height = 6)

# We see that we obtain the contraction over a range of values that does
# not include the best values obtained for HMC.
# This brings the question: what is the efficiency loss?
#
# function to compute distance between coupled chains
distance_ <- function(cchain){
  nsteps <- nrow(cchain$samples1) - 1
  return(sapply(1:nsteps, function(index) sqrt(sum((cchain$samples1[1+index,] - cchain$samples2[index,])^2))))
}
#
effectiveTime <- pi/2
nsteps <- 20
stepsize <- effectiveTime / nsteps
load(file = "mvnorm.cchains.RData")
summary(meetingtimes)

distances_ <- foreach (irep = 1:nrep, .combine = rbind) %dorng% {
  ds <- distance_(cchains[[irep]])
  if (length(ds) < max_iterations){
    ds <- c(ds, rep(0, max_iterations - length(ds)))
  }
  ds
}
dist.df <- melt(log(t(distances_[,1:350])))
tail(dist.df)
names(dist.df) <- c("iteration", "rep", "value")
g <- ggplot(dist.df, aes(x = iteration, y = value, group = rep)) + geom_line(alpha = 0.25)
g <- g + xlab('iteration') + ylab('logarithm of distance') + ylim(-45, 5)
g
ggsave(filename = "mvnorm.distancetraces.png", plot = g, width = 9, height = 6)
##

load(file = "mvnorm.withMH.cchains.RData")
sapply(cchains, function(x) x$meetingtime) %>% summary
distances_ <- foreach (irep = 1:nrep, .combine = rbind) %dorng% {
  ds <- distance_(cchains[[irep]])
  if (length(ds) < max_iterations){
    ds <- c(ds, rep(0, max_iterations - length(ds)))
  }
  ds
}
matplot(log(t(distances_)), type = "l", xlim = c(1,300), col = "black")
dist.df <- melt(log(t(distances_[,1:350])))
tail(dist.df)
names(dist.df) <- c("iteration", "rep", "value")
g <- ggplot(dist.df, aes(x = iteration, y = value, group = rep)) + geom_line(alpha = 0.25)
g <- g + xlab('iteration') + ylab('logarithm of distance') + ylim(-45, 5)
g
ggsave(filename = "mvnorm.distancetraces.withMH.png", plot = g, width = 9, height = 6)


### now compute efficiency
k <- 50
K <- 500
mean(sapply(cchains, function(x) x$iteration)) # computing cost is equal to K
mean_estimator1 <- as.numeric(foreach(irep = 1:nrep, .combine = c) %dopar% {
  H_bar(cchains[[irep]], h = function(x) x[1], k, K)
})

summary(mean_estimator1)
mean_target[1]
# efficiency
var(mean_estimator1) * K # to be compared with best asymptotic variance for HMC
# -> 1.719483
# for HMC, we had found, averaging over 5 repeats
fixednsteps.mcmc.df %>% select(irep, effectiveTime, mcmcvar) %>% group_by(effectiveTime) %>%
  summarize(averagemcmcvar = mean(mcmcvar))
# -> 0.1630441 for the optimal choice of trajectory length equal to pi

## we can plot histograms of the target as follows
component <- 1
hist1 <- histogram_c_chains(cchains, component, k, K)
g <- plot_histogram(hist1)
g <- g + stat_function(fun = function(x) dnorm(x, mean = mean_target[component], sd = sqrt(Sigma_target[component,component])))
g


