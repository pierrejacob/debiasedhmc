#include <Rcpp.h>
#include "mvnorm.h"

using namespace Rcpp;


// [[Rcpp::export]]
double logcosh(double x){
  double result = 0.;
  if (x > 0.){
    result = x + log(1.0 + exp(-2.0*x)) - 0.6931472;
  } else {
    result = -x + log(1.0 + exp(2.0*x)) - 0.6931472;
  }
  return result;
}


// [[Rcpp::export]]
NumericVector gaussian_max_couplingC(const NumericVector & mu1,
                                     const NumericVector & mu2,
                                     const NumericMatrix & Sigma1,
                                     const NumericMatrix & Sigma2){
    RNGScope scope;
    int p = Sigma1.cols();

    NumericMatrix x(p,2);

    NumericMatrix x1 = rmvnorm(1, mu1, Sigma1);
    x.column(0) = x1(_,0);

    double d1 = dmvnorm(x1, mu1, Sigma1)(0);
    double u1 = log(runif(1,0,1)(0));
    double d2 = dmvnorm(x1, mu2, Sigma2)(0);
    if ((u1 + d1) <= d2){
      x.column(1) = x1(_,0);
    } else {
      bool accept = FALSE;
      while (accept==FALSE){
        NumericMatrix x2 = rmvnorm(1, mu2, Sigma2);
        double d2 = dmvnorm(x2, mu2, Sigma2)(0);
        double u2 = log(runif(1,0,1)(0));
        double d1 = dmvnorm(x2, mu1, Sigma1)(0);
        if ((u2 + d2) > d1){
          accept = TRUE;
          x.column(1) = x2(_,0);
        }
      }
    }
    return x;
}

// [[Rcpp::export]]
NumericVector gaussian_max_coupling_cholesky(const NumericVector & mu1,
                                     const NumericVector & mu2,
                                     const Eigen::MatrixXd & Cholesky1,
                                     const Eigen::MatrixXd & Cholesky2,
                                     const Eigen::MatrixXd & Cholesky_inverse1,
                                     const Eigen::MatrixXd & Cholesky_inverse2){
  RNGScope scope;
  int p = Cholesky1.cols();

  NumericMatrix x(p,2);

  NumericMatrix x1 = rmvnorm_cholesky(1, mu1, Cholesky1);
  x.column(0) = x1(_,0);

  double d1 = dmvnorm_cholesky_inverse(x1, mu1, Cholesky_inverse1)(0);
  double u1 = log(runif(1,0,1)(0));
  double d2 = dmvnorm_cholesky_inverse(x1, mu2, Cholesky_inverse2)(0);
  if ((d1 + u1) <= d2){
    x.column(1) = x1(_,0);
  } else {
    bool accept = FALSE;
    while (accept==FALSE){
      NumericMatrix x2 = rmvnorm_cholesky(1, mu2, Cholesky2);
      double d2 = dmvnorm_cholesky_inverse(x2, mu2, Cholesky_inverse2)(0);
      double u2 = log(runif(1,0,1)(0));
      double d1 = dmvnorm_cholesky_inverse(x2, mu1, Cholesky_inverse1)(0);
      if ((u2 + d2) > d1){
        accept = TRUE;
        x.column(1) = x2(_,0);
      }
    }
  }
  return x;
}
