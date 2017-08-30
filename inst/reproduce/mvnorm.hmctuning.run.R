library(debiasedhmc)
library(coda)
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
#
nmcmc <- 5000
effectiveTimes <- c(pi/4, 2*pi/4, 3*pi/4, pi, 5*pi/4, 6*pi/4)
neff <- length(effectiveTimes)

nrep <- 5
nsteps <- 20
fixednsteps.mcmc.df <- data.frame()
for (i in 1:neff){
  effectiveTime <- effectiveTimes[i]
  cat("effective time:", effectiveTime, "\n")
  stepsize <- effectiveTime / nsteps
  hmc_kernel <- get_hmc_kernel(target$logtarget, target$gradlogtarget, stepsize, nsteps, dimension)
  fixednsteps.mcmc.df <- rbind(fixednsteps.mcmc.df, foreach(irep = 1:nrep, .combine = rbind) %dorng% {
    chain <- matrix(nrow = nmcmc, ncol = dimension)
    current_x <- rinit()
    for (iteration in 1:nmcmc){
      # if (iteration %% 1000 == 1) cat("iteration ", iteration, "\n")
      current_x <- hmc_kernel$kernel(current_x, iteration)$chain_state
      chain[iteration,] <- current_x
    }
    mcmcvar <- spectrum0(chain[,1])$spec
    mcmcess <- effectiveSize(chain[,1])
    data.frame(irep = irep, effectiveTime = effectiveTime, stepsize = stepsize, nsteps = nsteps, mcmcvar = mcmcvar, mcmcess = mcmcess)
  })
  save(nmcmc, nrep, effectiveTimes, nsteps,  fixednsteps.mcmc.df, file = "mvnorm.hmctuning.RData")
}
# load(file = "mvnorm.hmctuning.RData")


