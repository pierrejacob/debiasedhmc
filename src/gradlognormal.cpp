#include <Rcpp.h>
using namespace Rcpp;

//' @export
// [[Rcpp::export]]
NumericVector gradlognormal(const NumericVector & x, const NumericVector & mean, const NumericMatrix & precision) {
  int d = x.size();
  NumericVector xc = mean - x;
  NumericVector output(d);

  for(int i = 0; i < d; ++i){
    double total = 0;
    for(int j = 0; j < d; ++j){
      total += xc(j) * precision(j,i);
    }
    output[i] = total;
  }
  return(output);
}
