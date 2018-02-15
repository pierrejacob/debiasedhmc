#'@rdname get_mvnormal
#'@title Get Multivariate Normal target
#'@description Defines a multivariate Normal target
#' with prescribed dimension, mean, variance.
#'@param dimension dimension of the target
#'@param mu1 mean
#'@param Sigma_target variance
#'@return list with keys \code{logtarget} and \code{gradlogtarget},
#'ready to be used e.g. in \code{\link{get_hmc_kernel}}.
#'@export
get_mvnormal <- function(dimension, mean_target, Sigma_target){
  # compute precision matrix
  precision <- solve(Sigma_target)
  # compute Cholesky factor of precision matrix
  precision_chol <- t(chol(precision))
  # log-density of multivariate Normal
  logtarget <- function(x){
    return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), mean_target, precision_chol))
  }
  # gradient of log-density of multivariate Normal
  gradlogtarget <- function(x){
    return((- precision %*% matrix(x - mean_target, ncol = 1))[,1])
  }
  return(list(logtarget = logtarget, gradlogtarget = gradlogtarget))
}
