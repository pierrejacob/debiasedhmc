#'@rdname rnorm_max_coupling
#'@title Maximal coupling of two univariate Normal distributions
#'@description Sample from maximal coupling of two univariate Normal distributions,
#'specified through their means and standard deviations.
#'@param mu1 mean of first distribution
#'@param mu2 mean of second distribution
#'@param sigma1 standard deviation of first distribution
#'@param sigma2 standard deviation of second distribution
#'
#'@export
rnorm_max_coupling <- function(mu1, mu2, sigma1, sigma2){
  x <- rnorm(1, mu1, sigma1)
  if (dnorm(x, mu1, sigma1, log = TRUE) + log(runif(1)) < dnorm(x, mu2, sigma2, log = TRUE)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rnorm(1, mu2, sigma2)
      reject <- (dnorm(y, mu2, sigma2, log = TRUE) + log(runif(1)) < dnorm(y, mu1, sigma1, log = TRUE))
    }
    return(c(x,y))
  }
}



# from gaussian_max_coupling ----------------------------------------------

#'@rdname gaussian_max_coupling_cholesky_R
#'@title Maximal coupling of two multivariate Normal distributions
#'@description Sample from maximal coupling of two multivariate Normal distributions,
#'specified through their means, the Cholesky factors of their covariance matrices,
#'and the Cholesky factors of the inverse covariance matrices (i.e. the precision matrices).
#'@param mu1 mean of first distribution
#'@param mu2 mean of second distribution
#'@param Cholesky1 Cholesky factor of variance of first distribution
#'@param Cholesky2 Cholesky factor of variance of second distribution
#'@param Cholesky_inverse1 Cholesky factor of precision of first distribution
#'@param Cholesky_inverse2 Cholesky factor of precision of second distribution
#'
#'@export
gaussian_max_coupling_cholesky_R <- function(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2){
  # we need Cholesky <- chol(Sigma), not the transpose
  valid_Cholesky1 <- all(Cholesky1[lower.tri(Cholesky1)]==0)
  valid_Cholesky2 <- all(Cholesky2[lower.tri(Cholesky2)]==0)
  stopifnot(valid_Cholesky1, valid_Cholesky2)

  return(gaussian_max_coupling_cholesky(mu1, mu2, Cholesky1, Cholesky2, Cholesky_inverse1, Cholesky_inverse2))
}

#'@rdname gaussian_max_coupling
#'@title Maximal coupling of two multivariate Normal distributions
#'@description Sample from maximal coupling of two multivariate Normal distributions,
#'specified through their means and covariance matrices.
#'@param mu1 mean of first distribution
#'@param mu2 mean of second distribution
#'@param Sigma1 Variance of first distribution
#'@param Sigma2 Variance of second distribution
#'@export
gaussian_max_coupling <- function(mu1, mu2, Sigma1, Sigma2){
  return(gaussian_max_couplingC(mu1, mu2, Sigma1, Sigma2))
}
