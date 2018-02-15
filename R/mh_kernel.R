#'@rdname get_mh_kernel
#'@title Get random walk Metropolis-Hastings kernels
#'@description This function takes descriptions of the target
#' and tuning parameters, and returns a list containing the keys
#' \code{single_kernel}, \code{coupled_kernel} corresponding to marginal
#' and coupled MH kernels, with Normal random walks.
#' These kernels can then be used in the function \code{\link{coupled_chains}}.
#'@param logtarget function to compute target log-density, e.g. see \code{\link{get_mvnormal}}
#'@param Sigma_proposal variance of the Normal random walk for the proposals
#'@param dimension dimension of the target distribution
#'@return A list containing the keys
#' \code{single_kernel}, \code{coupled_kernel}.
#'@export
get_mh_kernel <- function(logtarget, Sigma_proposal, dimension){
  proposal_std <- sqrt(diag(Sigma_proposal))
  Sigma1_chol <- diag(proposal_std, dimension, dimension)
  Sigma1_chol_inv <- diag(1 / proposal_std, dimension, dimension)
  zeromean <- rep(0, dimension)

  # single kernel
  kernel <- function(chain_state, current_pdf, iteration){
    proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]
    proposal_pdf <- logtarget(proposal_value)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value, current_pdf = proposal_pdf, accept = accept))
    } else {
      return(list(chain_state = chain_state, current_pdf = current_pdf, accept = accept))
    }
  }

  # coupled kernel
  coupled_kernel <- function(chain_state1, chain_state2, current_pdf1, current_pdf2, iteration){
    # sample from maximal coupling
    proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                       Sigma1_chol, Sigma1_chol,
                                                       Sigma1_chol_inv, Sigma1_chol_inv)
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]
    # overlap <- all(proposal1 == proposal2) # keep track of stuff
    proposal_pdf1 <- logtarget(proposal1)
    proposal_pdf2 <- logtarget(proposal2)

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
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2,
                current_pdf1 = current_pdf1, current_pdf2 = current_pdf2,
                accept1 = accept1, accept2 = accept2))
                # overlap = overlap))
  }

  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}
