#include <RcppEigen.h>
using namespace Rcpp;

// [[Rcpp::export]]
double logistic_logtarget_c(NumericVector chain_state, const NumericVector & Y, const NumericMatrix & X, double lambda){
  int p = X.cols();
  int n = X.rows();
  double alpha = chain_state(0);
  NumericVector beta(p);
  for (int j = 0; j < p; j++){
    beta(j) = chain_state(j+1);
  }
  double logsigma2 = chain_state(p+1);
  double sigma2 = exp(logsigma2);
  NumericVector xbeta(n);
  double eval = 0;
  NumericVector logvalues(n);
  for (int i = 0; i < n; i ++){
    xbeta(i) = alpha;
    for (int j = 0; j < p; j ++){
      xbeta(i) += X(i,j) * beta(j);
    }
    eval -= xbeta(i) * (1. - Y(i));
    logvalues(i) = log(1 + exp(- xbeta(i)));
  }
  // double eval = - sum(xbeta * (1-Y));
  // NumericVector logvalues = log(1 + exp(- xbeta));
  double maxlogvalues = max(logvalues);
  eval -= (sum(logvalues - maxlogvalues) + n * maxlogvalues);
  eval += - 0.5 * (alpha*alpha + sum(beta*beta)) / sigma2 - 0.5 * ((double) p) * logsigma2 - lambda * sigma2;
  eval += logsigma2; // Jacobian term, from reparametrization
  return eval;
}


// [[Rcpp::export]]
NumericVector logistic_gradlogtarget_c(NumericVector chain_state, NumericVector Y, NumericMatrix X, double lambda){
  int p = X.cols();
  int n = X.rows();
  double alpha = chain_state(0);
  NumericVector beta(p);
  for (int j = 0; j < p; j++){
    beta(j) = chain_state(j+1);
  }
  double logsigma2 = chain_state(p+1);
  double sigma2 = exp(logsigma2);
  NumericVector xbeta(n);
  // NumericVector precomputed_vec(n);
  double precomp = 0.;
  NumericVector gradient(p+2);
  // std::fill(gradient.begin(), gradient.end(), 0);
  gradient(0) = - alpha / sigma2;
  for (int j = 0; j < p; j++){
    gradient(j+1) = - beta(j) / sigma2;
  }
  for (int i = 0; i < n; i ++){
    xbeta(i) = alpha;
    for (int j = 0; j < p; j ++){
      xbeta(i) += X(i,j) * beta(j);
    }
    // precomputed_vec(i) = (1. + exp(xbeta(i)));
    precomp = (1. + exp(xbeta(i)));
    gradient(0) += (Y(i)-1) + 1. / precomp;
    for (int j = 0; j < p; j++){
      gradient(j+1) += (Y(i)-1) * X(i,j) + X(i,j) / precomp;
    }
  }
  gradient(p+1) = 0.5 * (alpha*alpha + sum(beta*beta)) / sigma2 - 0.5 * ((double) p) - lambda * sigma2 + 1.;
  return gradient;
}
