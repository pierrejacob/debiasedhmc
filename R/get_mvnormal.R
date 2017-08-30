#'@export
get_mvnormal <- function(dimension, mean, Sigma_target){
  # compute precision matrix
  precision <- solve(Sigma_target)
  # compute Cholesky factor of precision matrix
  precision_chol <- t(chol(precision))
  # log-density of multivariate Normal
  logtarget <- function(x){
    return(fast_dmvnorm_chol_inverse(matrix(x, nrow = 1), mean, precision_chol))
  }
  # gradient of log-density of multivariate Normal
  gradlogtarget <- function(x){
    return((- precision %*% matrix(x - mean_target, ncol = 1))[,1])
  }
  return(list(logtarget = logtarget, gradlogtarget = gradlogtarget))
}
