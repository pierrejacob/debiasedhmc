library(debiasedhmc)
rm(list = setdiff(setdiff(ls(), "scriptfolder"), "resultsfolder"))
set.seed(18)
setmytheme()
registerDoParallel(cores = detectCores())
# define target of dimension p
dimension <- 250
mean_target <- rep(0, dimension)
Sigma_target <- diag(1, dimension, dimension)
for (i in 1:dimension){
  for (j in 1:dimension){
    Sigma_target[i,j] <- exp(-abs(i-j))
  }
}
target <- get_mvnormal(dimension, mean_target, Sigma_target)
# initial distribution of chains
rinit <- function() fast_rmvnorm(1, mean_target, Sigma_target)
# function to compute distance between coupled chains
distance_ <- function(cchain){
  nsteps <- nrow(cchain$samples1) - 1
  return(sapply(1:nsteps, function(index) sqrt(sum((cchain$samples1[1+index,] - cchain$samples2[index,])^2))))
}
#
effectiveTime <- pi/2
nsteps <- 20
stepsize <- effectiveTime / nsteps
# HMC kernels
hmc_kernel <- get_hmc_kernel(target$logtarget, target$gradlogtarget, stepsize, nsteps, dimension)
# MH kernels
omega <- 0.1
Sigma_proposal <- 1e-10 * diag(1, dimension, dimension)
mh_kernel <- get_mh_kernel(target$logtarget, Sigma_proposal, dimension)
# Mixture kernels
mixture_kernel <- function(chain_state, iter){
  if (runif(1) < omega){
    return(mh_kernel$kernel(chain_state, iter))
  } else {
    return(hmc_kernel$kernel(chain_state, iter))
  }
}
#
mixture_coupled_kernel <- function(chain_state1, chain_state2, iter){
  if (runif(1) < omega){
    return(mh_kernel$coupled_kernel(chain_state1, chain_state2, iter))
  } else {
    return(hmc_kernel$coupled_kernel(chain_state1, chain_state2, iter))
  }
}

nrep <- 100
max_iterations <- 1000
K <- 500
cchains <- foreach(irep = 1:nrep) %dorng% {
  coupled_chains(mixture_kernel, mixture_coupled_kernel, rinit, K = K, max_iterations = max_iterations)
}
save(K, max_iterations, nrep, effectiveTime, nsteps, stepsize, cchains, file = "mvnorm.withMH.cchains.RData")
# load(file = "mvnorm.withMH.cchains.RData")

