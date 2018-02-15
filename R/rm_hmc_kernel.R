#'@rdname get_rm_hmc_kernel
#'@title Get Riemann manifold Hamiltonian Monte Carlo kernels
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
#'@param metric metric tensor
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_rm_hmc_kernel <- function(logtarget, gradlogtarget, stepsize, nsteps, dimension, metric){
  # leap frog integrator
  # note that some papers use the notation U for - logtarget, so that there are minus signs everywhere
  leapfrog <- function(x, v){
    v <- v + stepsize * gradlogtarget(x) / 2
    for (step in 1:nsteps){
      x <- x + stepsize * as.numeric(metric$inverse %*% v)
      if (step != nsteps){
        v <- v + stepsize * gradlogtarget(x)
      }
    }
    v <- v + stepsize * gradlogtarget(x) / 2
    # we could negate the momentum but we don't use it here
    return(list(x = x, v = v))
  }

  # One step of RM-HMC
  kernel <- function(chain_state, current_pdf, iteration){
    current_v <- as.numeric(metric$inverse_chol_inverse %*% rnorm(dimension)) # velocity or momentum
    leapfrog_result <- leapfrog(chain_state, current_v)
    proposed_v <- - leapfrog_result$v
    proposed_x <- leapfrog_result$x

    proposed_pdf <- logtarget(proposed_x)
    accept_ratio <- proposed_pdf - current_pdf
    # the acceptance ratio also features the "kinetic energy" term of the extended target
    accept_ratio <- accept_ratio + 0.5 * current_v %*% metric$inverse %*% current_v -
                                   0.5 * proposed_v %*% metric$inverse %*% proposed_v
    accept <- FALSE
    if (log(runif(1)) < accept_ratio){
      chain_state <- proposed_x
      current_v <- proposed_v
      current_pdf <- proposed_pdf
      accept <- TRUE
    } else {
    }
    return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
  }

  # One step of coupled RM-HMC
  coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration){
    current_v <- as.numeric(metric$inverse_chol_inverse %*% rnorm(dimension)) # velocity or momentum, shared by both chains
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
    current_kinetic_energy <- 0.5 * current_v %*% metric$inverse %*% current_v

    accept_ratio1 <- accept_ratio1 + current_kinetic_energy - 0.5 * proposed_v1 %*% metric$inverse %*% proposed_v1
    accept_ratio2 <- accept_ratio2 + current_kinetic_energy - 0.5 * proposed_v2 %*% metric$inverse %*% proposed_v2

    logu <- log(runif(1)) # shared by both chains
    if (logu < accept_ratio1){
      chain_state1 <- proposed_x1
      current_pdf1 <- proposed_pdf1
    }
    if (logu < accept_ratio2){
      chain_state2 <- proposed_x2
      current_pdf2 <- proposed_pdf2
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
                current_pdf1 = current_pdf1, current_pdf2 = current_pdf2))
  }

  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}
