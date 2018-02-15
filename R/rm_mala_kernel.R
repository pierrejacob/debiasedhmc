#'@rdname rm_mala_kernel
#'@title Get Riemann manifold coupled Metropolis adjusted Langevin algorithm (MALA) kernels
#'that switches between synchrous coupling and maximal coupling
#'@description This function takes descriptions of the target
#' and tuning parameters, and returns a list containing the keys
#' \code{single_kernel}, \code{coupled_kernel} corresponding to marginal
#' and coupled MALA kernels. These kernels can then be used in the function \code{\link{coupled_chains}}.
#'
#'@param logtarget function to compute target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param gradlogtarget function to compute gradient of target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param stepsize standard deviation of Gaussian increment
#'@param dimension dimension of the target distribution
#'@param adaptive_tol threshold below which maximal coupling is triggered
#'@param metric metric tensor
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_rm_mala_kernel <- function(logtarget, gradlogtarget, stepsize, dimension, adaptive_tol, metric){
  stepsize_sq <- stepsize^2
  # matrices to compute maximal coupling
  Sigma_chol <- stepsize * t(metric$chol_inverse)
  inv_Sigma_chol <- metric$inverse_chol_inverse / stepsize

  # One step of RM-MALA
  kernel <- function(chain_state, iteration){
    # compute forward move
    forward_euler <- chain_state + 0.5 * stepsize_sq *
      as.numeric(metric$inverse %*% gradlogtarget(chain_state))

    # sample proposal
    proposal <- forward_euler + stepsize * as.numeric(metric$chol_inverse %*% rnorm(dimension))

    # compute backward move
    backward_euler <- proposal + 0.5 * stepsize_sq *
      as.numeric(metric$inverse %*% gradlogtarget(proposal))

    # compute acceptance probability
    accept_ratio <- logtarget(proposal) - logtarget(chain_state) -
      0.5 * (chain_state - backward_euler) %*% metric$tensor %*% (chain_state - backward_euler) / stepsize_sq +
      0.5 * (proposal - forward_euler) %*% metric$tensor %*% (proposal - forward_euler) / stepsize_sq

    # accept or reject proposal
    accept <- FALSE
    if (log(runif(1)) < accept_ratio){
      chain_state <- proposal
      accept <- TRUE
    } else {
    }
    return(list(chain_state = chain_state, accept = accept))
  }

  # one step of coupled RM-MALA
  coupled_kernel <- function(chain_state1, chain_state2, iteration){
    # compute forward move
    forward_euler1 <- chain_state1 + 0.5 * stepsize_sq *
      as.numeric(metric$inverse %*% gradlogtarget(chain_state1))
    forward_euler2 <- chain_state2 + 0.5 * stepsize_sq *
      as.numeric(metric$inverse %*% gradlogtarget(chain_state2))

    distance_between_chains <- sqrt( sum( (chain_state1 - chain_state2)^2 ) )
    if (distance_between_chains > adaptive_tol){
      # synchrous coupling
      increment1 <- as.numeric(metric$chol_inverse %*% rnorm(dimension))
      increment2 <- increment1

      # sample proposals
      proposal1 <- forward_euler1 + stepsize * increment1
      proposal2 <- forward_euler2 + stepsize * increment2
    } else {
      # sample proposals using maximal coupling
      proposal_value <- gaussian_max_coupling_cholesky_R(forward_euler1, forward_euler2,
                                                         Sigma_chol, Sigma_chol,
                                                         inv_Sigma_chol, inv_Sigma_chol)
      proposal1 <- proposal_value[,1]
      proposal2 <- proposal_value[,2]

    }

    # compute backward move
    backward_euler1 <- proposal1 + 0.5 * stepsize_sq *
      as.numeric(metric$inverse %*% gradlogtarget(proposal1))
    backward_euler2 <- proposal2 + 0.5 * stepsize_sq *
      as.numeric(metric$inverse %*% gradlogtarget(proposal2))

    # compute acceptance probability
    accept_ratio1 <- logtarget(proposal1) - logtarget(chain_state1) -
      0.5 * (chain_state1 - backward_euler1) %*% metric$tensor %*% (chain_state1 - backward_euler1) / stepsize_sq +
      0.5 * (proposal1 - forward_euler1) %*% metric$tensor %*% (proposal1 - forward_euler1) / stepsize_sq

    accept_ratio2 <- logtarget(proposal2) - logtarget(chain_state2) -
      0.5 * (chain_state2 - backward_euler2) %*% metric$tensor %*% (chain_state2 - backward_euler2) / stepsize_sq +
      0.5 * (proposal2 - forward_euler2) %*% metric$tensor %*% (proposal2 - forward_euler2) / stepsize_sq

    # accept or reject proposals
    logu <- log(runif(1)) # shared by both chains
    if (logu < accept_ratio1){
      chain_state1 <- proposal1
    } else {
    }
    if (logu < accept_ratio2){
      chain_state2 <- proposal2
    } else {
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}
