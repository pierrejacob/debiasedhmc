#'@export
get_mh_kernel <- function(logtarget, Sigma_proposal, dimension){
  Sigma1_chol <- chol(Sigma_proposal)
  Sigma1_chol_inv <- solve(chol(Sigma_proposal))
  Sigma2_chol <- chol(Sigma_proposal)
  Sigma2_chol_inv <- solve(chol(Sigma_proposal))
  zeromean <- rep(0, dimension)
  # single kernel
  kernel <- function(chain_state, iteration){
    proposal_value <- chain_state + fast_rmvnorm_chol(1, zeromean, Sigma1_chol)[1,]
    proposal_pdf <- logtarget(proposal_value)
    current_pdf <- logtarget(chain_state)
    accept <- (log(runif(1)) < (proposal_pdf - current_pdf))
    if (accept){
      return(list(chain_state = proposal_value))
    } else {
      return(list(chain_state = chain_state))
    }
  }
  coupled_kernel <- function(chain_state1, chain_state2, iteration){
    distance_ <- mean((chain_state1 - chain_state2)^2)
    # proposal_value <- gaussian_max_coupling(chain_state1, chain_state2, Sigma_proposal, Sigma_proposal)
    proposal_value <- gaussian_max_coupling_cholesky_R(chain_state1, chain_state2,
                                                       Sigma1_chol, Sigma2_chol, Sigma1_chol_inv, Sigma2_chol_inv)
    proposal1 <- proposal_value[,1]
    proposal2 <- proposal_value[,2]
    proposal_pdf1 <- logtarget(proposal1)
    proposal_pdf2 <- logtarget(proposal2)
    current_pdf1 <- logtarget(chain_state1)
    current_pdf2 <- logtarget(chain_state2)
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
    }
    if (accept2){
      chain_state2 <- proposal2
    }
    return(list(chain_state1 = chain_state1, chain_state2 = chain_state2))
  }
  return(list(kernel = kernel, coupled_kernel = coupled_kernel))
}
