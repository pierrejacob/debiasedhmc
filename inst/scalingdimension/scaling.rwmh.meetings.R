rm(list = ls())
library(tictoc)
library(debiasedhmc)
library(ggplot2)
library(doParallel)
library(doRNG)
registerDoParallel(cores = detectCores()-2)
setmytheme()
set.seed(1)
# scaling of RWMH on standard Normal target in arbitrary dimension

# initialize chains from the target (i.e. at equilibrium, i.e. no burn-in)
rinit <- function(dimension){
  return(rnorm(dimension))
}

# variance = 1 on each component, 0 off diagonal
# target log-density evaluation
target <- function(x) -0.5 * sum(x^2)

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

# scaling prescribed in Roberts et al.
scaling_stepsize <- function(constant, dimension){
  return(constant * (1 / dimension)^(1/2)) # MALA
}

samplemeetingtime <- function(mh_stepsize, dimension, max_iterations = Inf){
  mixture_single_kernel <- function(chain_state, current_pdf){
    return(mh_single_kernel(chain_state, current_pdf, mh_stepsize))
  }
  #
  mixture_coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2){
    return(mh_coupled_kernel(chain_state1, chain_state2, current_pdf1, current_pdf2, mh_stepsize))
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

constants <- c(1, 1.5, 2)
nrep <- 100 # replace by 1000 to get results of the paper
dimensions <- seq(from = 1, to = 10, by = 1)
results.df <- data.frame()
for (constant in constants){
  for (dimension in dimensions){
    mh_stepsize <- scaling_stepsize(constant, dimension)
    print(mh_stepsize)
    results_ <- foreach(irep = 1:nrep) %dorng% {
      samplemeetingtime(mh_stepsize, dimension, max_iterations = Inf)
    }
    meetingtimes <- sapply(results_, function(x) x$meetingtime)
    results.df <- rbind(results.df, data.frame(constant = constant,
                                               dimension = dimension,
                                               mean_meetingtime = mean(meetingtimes)))
    # save(results.df, file = "output.scaling.rwmh.meetings.RData")

  }
}
# results.df

save(results.df, file = "scaling.rwmh.meetings.v3.RData")

ggplot(results.df, aes(x = dimension, y = mean_meetingtime, group = constant, linetype = factor(constant))) + geom_line() +
  scale_linetype(name = "C =") + ylab("Average meeting time") + xlab("Dimension") + ylim(0,3000)
#
