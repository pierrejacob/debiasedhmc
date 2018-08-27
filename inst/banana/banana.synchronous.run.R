rm(list = ls())
library(debiasedhmc)
library(doParallel)
library(doRNG)
library(ggplot2)
library(gridExtra)
setmytheme()
registerDoParallel(cores = detectCores()-2)
set.seed(1)

# Banana distribution
dimension <- 2

#
# target log-density evaluation
target <- function(x) -(1-x[1])^2 - 10*((x[2]-x[1]^2)^2)
# gradient of target log-density
gradtarget <- function(x) c(-2*(x[1]-1) + 40*x[1]*(x[2]-x[1]^2),
                            -20 * (x[2]-x[1]^2))

## initialize chains from the target (i.e. at equilibrium, i.e. no burn-in)
rinit <- function(dimension){
  return(c(runif(1, -5, 5), runif(1, -5, 5)))
}

# leap frog
leapfrog <- function(x, v, nsteps, stepsize){
  v <- v + stepsize * gradtarget(x) / 2
  for (step in 1:nsteps){
    x <- x + stepsize * v
    if (step != nsteps){
      v <- v + stepsize * gradtarget(x)
    }
  }
  v <- v + stepsize * gradtarget(x) / 2
  # we could negate the momentum but we don't use it here
  return(list(x = x, v = v))
}

# One step of HMC
hmc_single_kernel <- function(chain_state, current_pdf, nsteps, stepsize){
  current_v <- rnorm(length(chain_state)) # momentum/velocity (mass is one)
  leapfrog_result <- leapfrog(chain_state, current_v, nsteps, stepsize)
  proposed_v <- - leapfrog_result$v
  proposed_x <- leapfrog_result$x
  proposed_pdf <- target(proposed_x)
  accept_ratio <- proposed_pdf - current_pdf
  # the acceptance ratio also features the "kinetic energy" term of the extended target
  accept_ratio <- accept_ratio + sum(current_v^2) / 2 - sum(proposed_v^2) / 2
  accept <- FALSE
  if (is.finite(accept_ratio)){
    accept <- (log(runif(1)) < accept_ratio)
  }
  if (accept){
    chain_state <- proposed_x
    current_pdf <- proposed_pdf
    accept <- TRUE
  }
  return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
}

# One step of coupled HMC with synchronous coupling
hmc_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, nsteps, stepsize){
  current_v <- rnorm(length(chain_state1)) # velocity or momentum, shared by both chains
  leapfrog_result1 <- leapfrog(chain_state1, current_v, nsteps, stepsize)
  leapfrog_result2 <- leapfrog(chain_state2, current_v, nsteps, stepsize)
  proposed_v1 <- - leapfrog_result1$v
  proposed_x1 <- leapfrog_result1$x
  proposed_v2 <- - leapfrog_result2$v
  proposed_x2 <- leapfrog_result2$x
  # if we were a bit smarter, we would save logtarget(current_x) so as
  # to not re-evaluate it at every step
  proposed_pdf1 <- target(proposed_x1)
  proposed_pdf2 <- target(proposed_x2)
  accept_ratio1 <- proposed_pdf1 - current_pdf1
  accept_ratio2 <- proposed_pdf2 - current_pdf2
  # the acceptance ratio also features the "kinetic energy" term of the extended target
  accept_ratio1 <- accept_ratio1 + sum(current_v^2) / 2 - sum(proposed_v1^2) / 2
  accept_ratio2 <- accept_ratio2 + sum(current_v^2) / 2 - sum(proposed_v2^2) / 2
  if (is.na(accept_ratio1)){
    accept_ratio1 <- -Inf
  }
  if (is.na(accept_ratio2)){
    accept_ratio2 <- -Inf
  }
  logu <- log(runif(1)) # shared by both chains
  if (logu < accept_ratio1){
    chain_state1 <- proposed_x1
    current_pdf1 <- proposed_pdf1
  } else {
    #
  }
  if (logu < accept_ratio2){
    chain_state2 <- proposed_x2
    current_pdf2 <- proposed_pdf2
  } else {
  }
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2, current_pdf1 = current_pdf1, current_pdf2 = current_pdf2))
}

# MH kernel
mh_single_kernel <- function(chain_state, current_pdf, stepsize = 1){
  proposal_value <- chain_state + stepsize * rnorm(length(chain_state), mean = 0, sd = 1)
  proposal_pdf <- target(proposal_value)
  accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
  if (accept){
    return(list(chain_state = proposal_value, current_pdf = proposal_pdf, accept = accept))
  } else {
    return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
  }
}

# max coupling of Gaussians, centered at mu1 and mu2, and with covariance matrix stepsize^2 * Identity
gaussian_max_coupling <- function(mu1, mu2, stddev){
  d <- length(mu1)
  x <- rnorm(d, mu1, stddev)
  if (sum(dnorm(x, mu1, stddev, log = TRUE)) + log(runif(1)) < sum(dnorm(x, mu2, stddev, log = TRUE))){
    return(cbind(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rnorm(d, mu2, stddev)
      reject <- (sum(dnorm(y, mu2, stddev, log = TRUE)) + log(runif(1)) < sum(dnorm(y, mu1, stddev, log = TRUE)))
    }
    return(cbind(x,y))
  }
}

# coupled MH kernel
mh_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, stepsize){
  proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, stepsize)
  proposal1 <- proposal_value[,1]
  proposal2 <- proposal_value[,2]
  proposal_pdf1 <- target(proposal1)
  proposal_pdf2 <- target(proposal2)
  logu <- log(runif(1))
  accept1 <- FALSE
  accept2 <- FALSE
  if (is.finite(proposal_pdf1)){
    accept1 <- (logu < (proposal_pdf1 - current_pdf1))
  }
  if (is.finite(proposal_pdf2)){
    accept2 <- (logu < (proposal_pdf2 - current_pdf2))
  }
  if (accept1){
    chain_state1 <- proposal1
    current_pdf1 <- proposal_pdf1
  }
  if (accept2){
    chain_state2 <- proposal2
    current_pdf2 <- proposal_pdf2
  }
  return(list(chain_state1 = chain_state1, chain_state2 = chain_state2, current_pdf1 = current_pdf1, current_pdf2 = current_pdf2))
}

samplemeetingtime <- function(omega, mh_stepsize, hmc_stepsize, hmc_nsteps, max_iterations = Inf){
  # Mixture kernels
  mixture_single_kernel <- function(chain_state, current_pdf){
    if (runif(1) < omega){
      return(mh_single_kernel(chain_state, current_pdf, mh_stepsize))
    } else {
      return(hmc_single_kernel(chain_state, current_pdf, hmc_nsteps, hmc_stepsize))
    }
  }
  #
  mixture_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2){
    if (runif(1) < omega){
      return(mh_coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, mh_stepsize))
    } else {
      return(hmc_coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, hmc_nsteps, hmc_stepsize))
    }
  }
  #
  chain_state1 <- rinit(dimension)
  chain_state2 <- rinit(dimension)
  current_pdf1 <- target(chain_state1)
  current_pdf2 <- target(chain_state2)
  sres1 <- mixture_single_kernel(chain_state1, current_pdf1)
  chain_state1 <- sres1$chain_state
  current_pdf1 <- sres1$current_pdf
  iter <- 1
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  # iter here is 1; at this point we have X_1,Y_0 and we are going to generate successively X_t,Y_{t-1} where iter = t
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    res_coupled_kernel <- mixture_coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)
    chain_state1 <- res_coupled_kernel$chain_state1
    chain_state2 <- res_coupled_kernel$chain_state2
    current_pdf1 <- res_coupled_kernel$current_pdf1
    current_pdf2 <- res_coupled_kernel$current_pdf2
    if (all(chain_state1 == chain_state2) && !meet){
      # recording meeting time tau
      meet <- TRUE
      meetingtime <- iter
    }
    # stop after tau steps
    if (iter >= meetingtime){
      finished <- TRUE
    }
  }
  return(list(meetingtime = meetingtime, iteration = iter, finished = finished))
}

coupled_chains <- function(single_kernel, coupled_kernel, rinit, m = 1, max_iterations = Inf, preallocate = 10){
  # initialize
  chain_state1 <- rinit(dimension)
  chain_state2 <- rinit(dimension)
  current_pdf1 <- target(chain_state1)
  current_pdf2 <- target(chain_state2)
  # pre-allocate
  p <- length(chain_state1)
  samples1 <- matrix(nrow = m+preallocate+1, ncol = p)
  samples2 <- matrix(nrow = m+preallocate, ncol = p)
  nrowsamples1 <- m+preallocate+1
  samples1[1,] <- chain_state1
  samples2[1,] <- chain_state2
  current_nsamples1 <- 1
  iter <- 1

  # move first chain
  res_single_kernel <- single_kernel(chain_state1, current_pdf1)
  chain_state1 <- res_single_kernel$chain_state
  current_pdf1 <- res_single_kernel$current_pdf
  current_nsamples1 <- current_nsamples1 + 1
  samples1[current_nsamples1,] <- chain_state1

  # iterate
  meet <- FALSE
  finished <- FALSE
  meetingtime <- Inf
  while (!finished && iter < max_iterations){
    iter <- iter + 1
    if (meet){
      # only need to use single kernel after meeting
      res_single_kernel <- single_kernel(chain_state1, current_pdf1)
      chain_state1 <- res_single_kernel$chain_state
      current_pdf1 <- res_single_kernel$current_pdf
      chain_state2 <- chain_state1
      current_pdf2 <- current_pdf1
    } else {
      # use coupled kernel
      res_coupled_kernel <- coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2)
      chain_state1 <- res_coupled_kernel$chain_state1
      chain_state2 <- res_coupled_kernel$chain_state2
      current_pdf1 <- res_coupled_kernel$current_pdf1
      current_pdf2 <- res_coupled_kernel$current_pdf2
      # check if meeting happens
      if (all(chain_state1 == chain_state2) && !meet){
        # recording meeting time tau
        meet <- TRUE
        meetingtime <- iter
      }
    }
    # store coupled chains
    if ((current_nsamples1+1) > nrowsamples1){
      new_rows <- nrowsamples1 - 1
      nrowsamples1 <- nrowsamples1 + new_rows
      samples1 <- rbind(samples1, matrix(NA, nrow = new_rows, ncol = p))
      samples2 <- rbind(samples2, matrix(NA, nrow = new_rows, ncol = p))
    }
    samples1[current_nsamples1+1,] <- chain_state1
    samples2[current_nsamples1,] <- chain_state2
    current_nsamples1 <- current_nsamples1 + 1
    # stop after max(m, tau) steps
    if (iter >= max(meetingtime, m)){
      finished <- TRUE
    }
  }
  # drop redundant entries
  samples1 <- samples1[1:current_nsamples1,,drop=F]
  samples2 <- samples2[1:(current_nsamples1-1),,drop=F]
  return(list(samples1 = samples1, samples2 = samples2,
              meetingtime = meetingtime, iteration = iter, finished = finished))
}

# probability associated with MH kernel
omega <- 0.05
hmc_nsteps <- 500
hmc_stepsize <- 1/hmc_nsteps
# std deviation of proposal in MH kernel
mh_stepsize <- 1e-3

# Mixture kernels
mixture_single_kernel <- function(chain_state, current_pdf){
  if (runif(1) < omega){
    return(mh_single_kernel(chain_state, current_pdf, mh_stepsize))
  } else {
    return(hmc_single_kernel(chain_state, current_pdf, hmc_nsteps, hmc_stepsize))
  }
}
#
mixture_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2){
  if (runif(1) < omega){
    return(mh_coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, mh_stepsize))
  } else {
    return(hmc_coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, hmc_nsteps, hmc_stepsize))
  }
}


nrep <- 1000
results.df <- foreach(irep = 1:nrep) %dorng% {
coupled_chains(mixture_single_kernel, mixture_coupled_kernel, rinit, max_iterations = 1e4)
}
save(results.df, nrep, omega, hmc_nsteps, hmc_stepsize, mh_stepsize, file = "banana.synchr.meetings.RData")
load(file = "banana.synchr.meetings.RData")
summary(sapply(results.df, function(x) x$meetingtime))
hist(sapply(results.df, function(x) x$meetingtime), nclass = 100)
#

