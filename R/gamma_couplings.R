#'@rdname rgamma_coupled
#'@title Sample from maximally coupled Gamma
#'@description Sample from maximally coupled Gamma distributions with given parameters alpha1, alpha2, beta1, beta2,
#'where alpha is the shape and beta the rate, i.e. the log-density at x is
#'
#' alpha * log(beta) - lgamma(alpha) + (alpha-1) * log(x) - beta x
#'
#'@param alpha1 shape parameter of first distribution
#'@param alpha2 shape parameter of second distribution
#'@param beta1 rate parameter of first distribution
#'@param beta2 rate parameter of second distribution
#'@export
rgamma_coupled <- function(alpha1, alpha2, beta1, beta2){
  x <- rgamma(1, alpha1, rate = beta1)
  if (dgamma(x, alpha1, rate = beta1, log = TRUE) + log(runif(1)) < dgamma(x, alpha2, rate = beta2, log = TRUE)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rgamma(1, alpha2, rate = beta2)
      reject <- (dgamma(y, alpha2, rate = beta2, log = TRUE) + log(runif(1)) < dgamma(y, alpha1, rate = beta1, log = TRUE))
    }
    return(c(x,y))
  }
}
