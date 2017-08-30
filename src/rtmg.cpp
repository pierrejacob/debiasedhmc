// This file is taken from the tmg package on CRAN, by Ari Pakman
//#include "rtmg.h"
#include "HmcSampler.h"

#include <RcppEigen.h>
#include <cstdlib>
#include <iostream>
#include <Eigen/Dense>


using namespace Eigen;
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix rtmg(int n_, int seed_, NumericVector initial_, int numlin_, NumericMatrix F_, NumericVector g_,
                   int numquad_, NumericMatrix quadratics_ ){


  const int n = n_;
  const Map<VectorXd> initial_value = as<Map<VectorXd> > (initial_);

  int dim = initial_value.rows();
  int seed = seed_;

  HmcSampler hmc1(dim, seed);


  const int numlin = numlin_;
  const int numquad = numquad_;



  if (numlin >0){
	  const Map<MatrixXd> F = as<Map<MatrixXd> > (F_);
	  const Map<VectorXd> g = as<Map<VectorXd> > (g_);

	  for(int i=0; i<numlin; i++){
		hmc1.addLinearConstraint(F.row(i),g(i));
  	}
   }


  if (numquad >0){
	  const Map<MatrixXd> Q = as<Map<MatrixXd> > (quadratics_);

	  for(int i=0; i<numquad; i++){
	 	MatrixXd A = Q.block( i*(dim+2), 0, dim, dim);  //block(firstRow, firstCol, rows, cols)
	        VectorXd B = Q.row(i*(dim+2) + dim);
	        double C = Q(i*(dim+2) + dim+1,0);
		hmc1.addQuadraticConstraint(A,B,C);
  		}

  }



  hmc1.setInitialValue(initial_value);

  MatrixXd samples(n,dim);

  for (int i=0; i<n; i++){
      samples.row(i) = hmc1.sampleNext();
 	 }


 return wrap(samples);

}

