rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(parallel)

# load model and dataset
load("inst/coxprocess/coxprocess_with_metric.RData")

# use the parallel RNG of L'Ecuyer et al (2002) to ensure that each processor
# has access to a stream of random numbers that do not overlap with other processors
RNGkind("L'Ecuyer-CMRG")

# index of job-array when submitting batch jobs on the compute cluster
# set igrid <- 1 for serial jobs
igrid <- as.integer(Sys.getenv('SLURM_ARRAY_TASK_ID'))
set.seed(1) # initial seed
for (i in 1:igrid){
  .Random.seed <- nextRNGStream(.Random.seed) # compute appropriate stream
}

# no. of repetitions desired for each processor
nreps <- 10

# guideline for the choice of k and m values for this application
k <- 71
m <- 710

# stepsize and no. of steps with good contraction and small average compute time
stepsize <- 0.13
nsteps <- 10

# define HMC kernel and coupled HMC kernel
hmc <- get_rm_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension, metric)

# define RWMH kernel and coupled RWMH kernel
omega <- 1 / 20 # probability of selecting coupled RWMH
Sigma_std <- 1e-3 # proposal standard deviation of RWMH
Sigma_proposal <- Sigma_std^2 * diag(1, dimension, dimension)
mh <- get_mh_kernel(logtarget, Sigma_proposal, dimension)

# define mixture kernel
mixture_kernel <- function(chain_state, current_pdf, iteration){
  if (runif(1) < omega){
    return(mh$kernel(chain_state, current_pdf, iteration))
  } else {
    return(hmc$kernel(chain_state, current_pdf, iteration))
  }
}

# define coupled mixture kernel
mixture_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration){
  if (runif(1) < omega){
    return(mh$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration))
  } else {
    return(hmc$coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration))
  }
}

# pre-allocate
runtimes <- rep(0, nreps)
meetingtime <- rep(0, nreps)
unbiased_estimates <- matrix(nrow = nreps, ncol = 2 * dimension) # first and second moment

for(irep in 1:nreps){
  # EITHER compute unbiased estimator while simulating coupled chains (lower memory requirements)
  tic()
  estimation_output <- unbiased_estimator(logtarget, mixture_kernel, mixture_coupled_kernel,
                                          rinit, h = function(x) c(x, x^2), k = k, m = m)
  timing <- toc()
  runtime <- timing$toc - timing$tic
  runtimes[irep] <- runtime
  meetingtime[irep] <- estimation_output$meetingtime
  unbiased_estimates[irep, ] <- estimation_output$uestimator

  # OR run coupled chains and compute unbiased estimator (higher memory requirements)
  # tic()
  #   cchains <- coupled_chains(logtarget, mixture_kernel, mixture_coupled_kernel, rinit, m = m)
  # timing <- toc()
  # runtime <- timing$toc - timing$tic
  # runtimes[irep] <- runtime
  # meetingtime[irep] <- cchains$meetingtime
  # unbiased_estimates[irep, ] <- H_bar(cchains, h = function(x) c(x, x^2), k = k, m = m)

  cat("Repetition:", irep, "/", nreps, "\n")
}

filename <- paste("inst/coxprocess/output/output.rm.hmc.repeat", igrid, ".RData", sep = "")
save(nreps, k, m, stepsize, nsteps, runtimes, meetingtime, unbiased_estimates,
     file = filename, safe = F)
