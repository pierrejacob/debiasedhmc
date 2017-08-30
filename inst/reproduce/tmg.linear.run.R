### This file has been created by modifying files from the tmg package on CRAN, by Ari Pakman
library(debiasedhmc)
registerDoParallel(cores = detectCores())
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
set.seed(1)
setmytheme()
par(mfrow = c(1,1))

#Define precision matrix and linear term
M = matrix(c(1,0,0,1), 2, 2)
r = c(4,4)
# Set initial point for the Markov chain
initial1 <- c(2,2.1)
initial2 <- c(2,2.1)
# Define two linear constraints
# first is 1.1 * x - 1 * y + 0 >= 0
f = matrix(nrow = 1, ncol = 2)
f[1,1] = 1.1
f[1,2] = -1
g = 0
# samples = rtmg_coupled_R(n, M, r, initial1, initial2, f, g, q=NULL);
# plot(samples$samples1, pch=".")

# second is -1 * x + 1 * y + 0 >= 0
f <- rbind(f, c(-1, 1))
g <- c(g, 0)

f %*% matrix(c(0.403504, 0.608504), ncol = 1) + g

# Set number of samples
n=10000
samples = rtmgR(n, M, r, initial1, f, g, q=NULL);
save(samples, file = "tmg.linear.samples.RData")
load("tmg.linear.samples.RData")
plot(samples, pch=".")

nrep <- 1000
K <- 100
cchains <- foreach(irep = 1:nrep) %dorng% {
  res <- rtmg_coupled_R(K, M, r, initial1, initial2, f, g, q = NULL);
  distances <- sapply(1:K, function(index) sqrt(sum(((res$samples1[index+1,] - res$samples2[index,]))^2)))
  list(samples1 = res$samples1, samples2 = res$samples2, meetingtime = which(distances == 0)[1],
       iteration = K)
}

save(nrep, K, initial1, initial2, cchains, file = "tmg.linear.cchains.RData")

