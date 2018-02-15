rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(coda)

# load model and dataset
load("inst/logistic/germancredit.RData")

# compute asymptotic variance
compute_variance <- function(trajectory){
  total_var <- sum(spectrum0.ar(trajectory)$spec) # variance of first moments
  total_var <- total_var + sum(spectrum0.ar(trajectory^2)$spec) # variance of second moments
  return(total_var)
}

# optimal HMC parameters
stepsize <- 0.03
nsteps <- 10

# define HMC kernel
hmc <- get_hmc_kernel(logtarget, gradlogtarget, stepsize, nsteps, dimension)

# number of HMC iterations and burn-in
nmcmc <- 11000
burnin <- 1001

# pre-allocate and initialize
chain <- matrix(nrow = nmcmc, ncol = dimension)
current_x <- rinit()
current_pdf <- logtarget(current_x)
accept <- 0
tic()
for (imcmc in 1:nmcmc){
  current_state <- hmc$kernel(current_x, current_pdf, imcmc)
  current_x <- current_state$chain_state
  current_pdf <- current_state$current_pdf
  chain[imcmc, ] <- current_x
  accept <- accept + current_state$accept
  cat("Iteration:", imcmc, "/", nmcmc, "\n")
}
timing <- toc()
runtime <- timing$toc - timing$tic

asymptotic_variance <- compute_variance(chain[burnin:nmcmc, ])
acceptprob <- accept / nmcmc
optimal_hmc <- list(asymptotic_variance = asymptotic_variance, acceptprob = acceptprob,
                    runtime = runtime)
filename <- "inst/logistic/output/optimal.hmc.germancredit.RData"
save(stepsize, nsteps, optimal_hmc, file = filename, safe = F)


