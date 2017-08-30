library(debiasedhmc)
library(coda)
library(latex2exp)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
set.seed(18)
setmytheme()
registerDoParallel(cores = detectCores())
# data
germancredit <- read.table(system.file("", "germancredit.txt", package = "debiasedhmc"))
X <- germancredit[,1:24]
X <- scale(X)
Y <- germancredit[,25]
Y <- Y - 1
n <- nrow(X)
p <- ncol(X)
Xinter <- matrix(nrow = n, ncol = p*(p-1) / 2)
index <- 1
for (j in 1:(p-1)){
  for (jprime in (j+1):p){
    Xinter[,index] <- X[,j] * X[,jprime]
    index <- index + 1
  }
}
Xinter <- scale(Xinter)
X <- cbind(X, Xinter)
p <- ncol(X)
lambda <- 0.01
dimension <- p + 2
#
## First, investigate choice of tuning parameter for HMC
load(file = "germancredit.hmctuning.RData")
nrep <- hmctraces %>% length
mcmcvar.df <- data.frame()
for (irep in 1:nrep){
  print(effectiveTimes)
  burnin <- 5000
  mcmcvar <- c()
  par(mfrow = c(2,3))
  for (i in 1:length(effectiveTimes)){
    cat("length ", effectiveTimes[i], "\n")
    chain <- hmctraces[[irep]][[i]]
    nmcmc <- nrow(chain)
    subiterations <- floor(seq(from = 1, to = nmcmc, length.out = min(nmcmc, 1000)))
    print('accept probability:')
    print(1-mean(diff(chain[,1]) == 0)) # proportion of rejected moves
    print("ESS:")
    print(as.numeric(effectiveSize(chain[burnin:nmcmc,1])))
    print("Variance:")
    print(spectrum0(chain[burnin:nmcmc,1])$spec)
    mcmcvar <- c(mcmcvar, spectrum0(chain[burnin:nmcmc,1])$spec)
    matplot(x = subiterations, y = chain[subiterations,1], type = "l", ylim = c(-10, 1))
  }
  mcmcvar.df <- rbind(data.frame(effectiveTimes = effectiveTimes, mcmcvar = mcmcvar, irep = rep(irep, length(mcmcvar))), mcmcvar.df)
}
g <- ggplot(mcmcvar.df, aes(x = effectiveTimes, y = mcmcvar, group = irep)) + geom_point()
g <- g + ylim(0, 0.5)
g <- g + ylab("HCMC variance") + xlab("trajectory length")
g
ggsave(filename = "germancredit.hmcvariance.trajectorylength.pdf", plot = g, width = 9, height = 6)

load(file = "germancredit.contraction.RData")
g <- ggplot(distance.df %>% filter(distance > 1e-30), aes(x = effectiveTime, y = distance)) + geom_point()
g <- g + scale_y_log10(breaks = c(1e-20, 1e-15, 1e-10, 1e-5, 1), limits = c(1e-20, 10))
g <- g + ylab("distance after 1000 iterations") + xlab("trajectory length")
g
ggsave(filename = "germancredit.distance100iter.trajectorylength.pdf", plot = g, width = 9, height = 6)

distance_ <- function(cchain){
  nsteps <- nrow(cchain$samples1) - 1
  return(sapply(1:nsteps, function(index) sqrt(sum((cchain$samples1[1+index,] - cchain$samples2[index,])^2))))
}

load(file = "germancredit.withMH.cchains.RData")
cchains.poorinit <- cchains
meetingtimes <- sapply(cchains, function(x) x$meetingtime)
summary(meetingtimes)
distances_ <- foreach (irep = 1:nrep, .combine = rbind) %dorng% {
  ds <- distance_(cchains[[irep]])
  if (length(ds) < max_iterations){
    ds <- c(ds, rep(0, max_iterations - length(ds)))
  }
  ds
}
# par(mfrow = c(1,1))
# matplot(log(t(distances_)), type = "l", xlim = c(1,2500), col = "black")
dist.df <- melt(log(t(distances_[,1:1000])))
tail(dist.df)
names(dist.df) <- c("iteration", "rep", "value")
g <- ggplot(dist.df, aes(x = iteration, y = value, group = rep)) + geom_line(alpha = 0.25)
g <- g + xlab('iteration') + ylab('logarithm of distance') + ylim(-20, 5)
g
ggsave(filename = "germancredit.distancetraces.png", plot = g, width = 9, height = 6)


load(file = "germancredit.withMH.cchains.betterinit.RData")
cchains.betterinit <- cchains
meetingtimes <- sapply(cchains, function(x) x$meetingtime)
summary(meetingtimes)
distances_ <- foreach (irep = 1:nrep, .combine = rbind) %dorng% {
  ds <- distance_(cchains[[irep]])
  if (length(ds) < max_iterations){
    ds <- c(ds, rep(0, max_iterations - length(ds)))
  }
  ds
}
# par(mfrow = c(1,1))
# matplot(log(t(distances_)), type = "l", xlim = c(1,500), col = "black")
dist.df <- melt(log(t(distances_[,1:1000])))
tail(dist.df)
names(dist.df) <- c("iteration", "rep", "value")
g <- ggplot(dist.df, aes(x = iteration, y = value, group = rep)) + geom_line(alpha = 0.25)
g <- g + xlab('iteration') + ylab('logarithm of distance') + ylim(-20, 5)
g
ggsave(filename = "germancredit.distancetraces.betterinit.png", plot = g, width = 9, height = 6)

### now, compute efficiencies, using larger number of repeats
# (takes a long time, around 1 minute, to load)
load(file = "germancredit.withMH.cchains.betterinit.1000.RData")
cchains.1000 <- cchains
k <- 100
K <- 1000
nrep <- length(cchains.1000)
mean_estimator1 <- as.numeric(foreach(irep = 1:nrep, .combine = c) %dopar% {
  H_bar(cchains.1000[[irep]], h = function(x) x[1], k, K)
})
summary(mean_estimator1)

expectedcost <- mean(sapply(cchains.1000, function(x) x$iteration))
expectedcost # equal to K
# efficiency
var(mean_estimator1) * K # to be compared with best asymptotic variance for HMC
# -> 0.3981529

# to be compared with the asymptotic variance of HMC,
# obtained from long traces
load("germancredit.hmclongtraces.RData")
hmctraces[[1]] %>% length
# first corresponds to trajectory length of 0.1
burnin <- 5000
chain.0.1 <- hmctraces[[1]][[1]]
print(spectrum0(chain.0.1[burnin:nrow(chain.0.1),1])$spec)
# -> 0.3264556
# second corresponds to trajectory length of 0.3
burnin <- 5000
chain.0.3 <- hmctraces[[1]][[2]]
print(spectrum0(chain.0.3[burnin:nrow(chain.0.3),1])$spec)
# -> 0.09299325

# we can compare to average of smaller traces over 5 repeats
mcmcvar.df %>% select(irep, effectiveTimes, mcmcvar) %>% group_by(effectiveTimes) %>%
  summarize(averagemcmcvar = mean(mcmcvar))

# -> 0.36018686 for epsilon * L = 0.1
# -> 0.09132089 for epsilon * L = 0.3
#
# ## we can plot histograms of the target as follows
nclass <- 50
hist1 <- histogram_c_chains(cchains.1000, 1, k, K, nclass = nclass)
histhmc <- hist(chain.0.3[burnin:nrow(chain.0.3),1], prob = TRUE, nclass = nclass, plot = FALSE)
histhmc.df <- data.frame(x = histhmc$mids, y = histhmc$density)
g <- plot_histogram(hist1)
g <- g + geom_line(data = histhmc.df, colour = "red", aes(x = x, y = y, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL, xend = NULL, yend = NULL))
g <- g + xlab(expression(alpha))
g
ggsave(filename = "germancredit.histogram1.pdf", plot = g, width = 9, height = 6)


hist2 <- histogram_c_chains(cchains.1000, 2, k, K, nclass = nclass)
histhmc <- hist(chain.0.3[burnin:nrow(chain.0.3),2], prob = TRUE, nclass = nclass, plot = FALSE)
histhmc.df <- data.frame(x = histhmc$mids, y = histhmc$density)
g <- plot_histogram(hist2)
g <- g + geom_line(data = histhmc.df, colour = "red", aes(x = x, y = y, xmin = NULL, xmax = NULL, ymin = NULL, ymax = NULL, xend = NULL, yend = NULL))
g <- g + xlab(expression(beta[1]))
g
ggsave(filename = "germancredit.histogram2.pdf", plot = g, width = 9, height = 6)

