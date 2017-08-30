#'@rdname digamma
#'@title Compute log-density of inverse gamma
#'@description Compute log-density of inverse gamma at x > 0
#' and with given parameters alpha, beta, given by
#'
#'  alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
#'
#'@export
digamma <- function(x, alpha, beta){
  return(alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x)
}

#'@rdname rigamma
#'@title Sample from inverse gamma
#'@description Sample from inverse gamma with given parameters alpha, beta,
#' with log-density at x given by
#'
#' alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
#'
#'@export
rigamma <- function(n, alpha, beta){
  return(1/rgamma(n = n, shape = alpha, rate = beta))
}

#'@rdname rigamma_coupled
#'@title Sample from maximally coupled inverse gamma
#'@description Sample from maximally coupled inverse gamma with given parameters alpha1, alpha2, beta1, beta2,
#'where the parametrization is that the log-density of IG(alpha, beta) at x is
#'
#' alpha * log(beta) - lgamma(alpha) - (alpha+1) * log(x) - beta / x
#'
#'@export
rigamma_coupled <- function(alpha1, alpha2, beta1, beta2){
  x <- rigamma(1, alpha1, beta1)
  if (digamma(x, alpha1, beta1) + log(runif(1)) < digamma(x, alpha2, beta2)){
    return(c(x,x))
  } else {
    reject <- TRUE
    y <- NA
    while (reject){
      y <- rigamma(1, alpha2, beta2)
      reject <- (digamma(y, alpha2, beta2) + log(runif(1)) < digamma(y, alpha1, beta1))
    }
    return(c(x,y))
  }
}

