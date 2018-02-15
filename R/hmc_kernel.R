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
    v <- v + stepsize * gradlogtarget(x) / 2
    for (step in 1:nsteps){
      x <- x + stepsize * v
      if (step != nsteps){
        v <- v + stepsize * gradlogtarget(x)
      }
    }
    v <- v + stepsize * gradlogtarget(x) / 2
    # we could negate the momentum but we don't use it here
    return(list(x = x, v = v))
  }

  # One step of HMC
  kernel <- function(chain_state, current_pdf, iteration){
    current_v <- rnorm(dimension) # velocity or momentum
    leapfrog_result <- leapfrog(chain_state, current_v)
    proposed_v <- - leapfrog_result$v
    proposed_x <- leapfrog_result$x
    proposed_pdf <- logtarget(proposed_x)
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

  # One step of coupled HMC
  coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration){
    current_v <- rnorm(dimension) # velocity or momentum, shared by both chains
    leapfrog_result1 <- leapfrog(chain_state1, current_v)
    leapfrog_result2 <- leapfrog(chain_state2, current_v)
    proposed_v1 <- - leapfrog_result1$v
    proposed_x1 <- leapfrog_result1$x
    proposed_v2 <- - leapfrog_result2$v
    proposed_x2 <- leapfrog_result2$x

    proposed_pdf1 <- logtarget(proposed_x1)
    proposed_pdf2 <- logtarget(proposed_x2)
    accept_ratio1 <- proposed_pdf1 - current_pdf1
    accept_ratio2 <- proposed_pdf2 - current_pdf2
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    current_kinetic_energy <- sum(current_v^2) / 2
    accept_ratio1 <- accept_ratio1 + current_kinetic_energy - sum(proposed_v1^2) / 2
    accept_ratio2 <- accept_ratio2 + current_kinetic_energy - sum(proposed_v2^2) / 2
    logu <- log(runif(1)) # shared by both chains
    accept1 <- FALSE
    accept2 <- FALSE
    if (is.finite(accept_ratio1)){
      accept1 <- (logu < accept_ratio1)
    }

    if (accept1){
      chain_state1 <- proposed_x1
      current_pdf1 <- proposed_pdf1
      accept1 <- TRUE
    }

    if (is.finite(accept_ratio2)){
      accept2 <- (logu < accept_ratio2)
    }

    if (accept2){
      chain_state2 <- proposed_x2
      current_pdf2 <- proposed_pdf2
      accept2 <- TRUE
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
                current_pdf1 = current_pdf1, current_pdf2 = current_pdf2,
                accept1 = accept1, accept2 = accept2))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}
