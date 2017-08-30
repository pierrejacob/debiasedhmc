#'@rdname get_hmc_kernel
#'@title Get Hamiltonian Monte Carlo kernels
#'@description This function takes descriptions of the target
#' and tuning parameters, and returns a list containing the keys
#' \code{single_kernel}, \code{coupled_kernel} corresponding to marginal
#' and coupled HMC kernels. These kernels can then be used in the function \code{\link{coupled_chains}}.
#'
#'@param logtarget function to compute target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param gradlogtarget function to compute gradient of target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param stepsize step size in the leap-frog integrator
#'@param nsteps number of leap-frog steps
#'@param dimension dimension of the target distribution
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_hmc_kernel <- function(logtarget, gradlogtarget, stepsize, nsteps, dimension){
  # leap frog integrator
  # note that some papers use the notation U for - logtarget, so that there are minus signs everywhere
  leapfrog <- function(x, v){
    xtraj <- matrix(nrow = nsteps+1, ncol = length(x))
    xtraj[1,] <- x
    v <- v + stepsize * gradlogtarget(x) / 2
    for (step in 1:nsteps){
      x <- x + stepsize * v
      xtraj[step+1,] <- x
      if (step != nsteps){
        v <- v + stepsize * gradlogtarget(x)
      }
    }
    v <- v + stepsize * gradlogtarget(x) / 2
    # we could negate the momentum but we don't use it here
    return(list(x = x, v = v, xtraj = xtraj))
  }
  # One step of HMC
  kernel <- function(chain_state, iteration){
    current_v <- rnorm(dimension) # velocity or momentum
    leapfrog_result <- leapfrog(chain_state, current_v)
    proposed_v <- - leapfrog_result$v
    proposed_x <- leapfrog_result$x
    # if we were a bit smarter, we would save logtarget(chain_state) so as
    # to not re-evaluate it at every step
    accept_ratio <- logtarget(proposed_x) - logtarget(chain_state)
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    accept_ratio <- accept_ratio + sum(current_v^2) / 2 - sum(proposed_v^2) / 2
    # we store the entire x-trajectory here, for plotting purposes
    # the first row stores the current value (which is also the last value of the latest accepted trajectory)
    xtraj <- matrix(nrow = nsteps+1, ncol = dimension)
    accept <- FALSE
    if (log(runif(1)) < accept_ratio){
      chain_state <- proposed_x
      current_v <- proposed_v
      xtraj <- leapfrog_result$xtraj
      accept <- TRUE
    } else {
      # if we reject, the entire x-trajectory is boring
      for (istep in 1:(nsteps+1)){
        xtraj[istep,] <- chain_state
      }
    }
    return(list(chain_state = chain_state, xtraj = xtraj, accept = accept))
  }
  # One step of HMC
  coupled_kernel <- function(chain_state1, chain_state2, iteration){
    current_v <- rnorm(dimension) # velocity or momentum, shared by both chains
    leapfrog_result1 <- leapfrog(chain_state1, current_v)
    leapfrog_result2 <- leapfrog(chain_state2, current_v)
    proposed_v1 <- - leapfrog_result1$v
    proposed_x1 <- leapfrog_result1$x
    proposed_v2 <- - leapfrog_result2$v
    proposed_x2 <- leapfrog_result2$x
    # if we were a bit smarter, we would save logtarget(current_x) so as
    # to not re-evaluate it at every step
    accept_ratio1 <- logtarget(proposed_x1) - logtarget(chain_state1)
    accept_ratio2 <- logtarget(proposed_x2) - logtarget(chain_state2)
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    accept_ratio1 <- accept_ratio1 + sum(current_v^2) / 2 - sum(proposed_v1^2) / 2
    accept_ratio2 <- accept_ratio2 + sum(current_v^2) / 2 - sum(proposed_v2^2) / 2
    # we store the entire x-trajectory here, for plotting purposes
    # the first row stores the current value (which is also the last value of the latest accepted trajectory)
    xtraj1 <- matrix(nrow = nsteps+1, ncol = dimension)
    xtraj2 <- matrix(nrow = nsteps+1, ncol = dimension)
    logu <- log(runif(1)) # shared by both chains
    if (logu < accept_ratio1){
      chain_state1 <- proposed_x1
      current_v1 <- proposed_v1
      xtraj1 <- leapfrog_result1$xtraj
    } else {
      # if we reject, the entire x-trajectory is boring
      for (istep in 1:(nsteps+1)){
        xtraj1[istep,] <- chain_state1
      }
    }
    if (logu < accept_ratio2){
      chain_state2 <- proposed_x2
      current_v2 <- proposed_v2
      xtraj2 <- leapfrog_result2$xtraj
    } else {
      for (istep in 1:(nsteps+1)){
        xtraj2[istep,] <- chain_state2
      }
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2, xtraj1 = xtraj1, xtraj2 = xtraj2))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}
