### This file has been created by modifying files from the tmg package on CRAN, by Ari Pakman
library(debiasedhmc)
registerDoParallel(cores = 6)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
set.seed(18)
setmytheme()
par(mfrow = c(1,1))

#Define precision matrix and linear term
M = matrix(c(1,0,0,1), 2, 2)
r = c(0,0)
# Set initial point for the Markov chain
initial1 <- c(2,0)
initial2 <- c(2,0)
# Define two quadratic constraints
A1 = matrix(c(-1/32,0,0,-1/8), 2, 2)
B1 = c(.25,.25)
C1 = 3/8
constr1 = list(A1,B1,C1)
#
A2 = matrix(c(4,-1,-1,8), 2, 2)
B2 = c(0,5)
C2 = -1
constr2 = list(A2,B2,C2)
q = list(constr1,constr2)

# Set number of samples
n=10000
samples = rtmgR(n, M, r, initial1, f=NULL, g=NULL, q=q);
save(samples, file = "tmg.quadratic.samples.RData")

nrep <- 1000
K <- 100
cchains <- foreach(irep = 1:nrep) %dorng% {
  res <- rtmg_coupled_R(K, M, r, initial1, initial2, f=NULL, g=NULL, q=q);
  distances <- sapply(1:K, function(index) sqrt(sum(((res$samples1[index+1,] - res$samples2[index,]))^2)))
  list(samples1 = res$samples1, samples2 = res$samples2, meetingtime = which(distances == 0)[1],
       iteration = K)
}

save(nrep, K, initial1, initial2, cchains, file = "tmg.quadratic.cchains.RData")

