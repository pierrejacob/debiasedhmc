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
List rtmg_coupled(int n_, int seed_, NumericVector initial_1, NumericVector initial_2, int numlin_, NumericMatrix F_, NumericVector g_,
                   int numquad_, NumericMatrix quadratics_ ){


  const int n = n_;
  const Map<VectorXd> initial_value1 = as<Map<VectorXd> > (initial_1);
  const Map<VectorXd> initial_value2 = as<Map<VectorXd> > (initial_2);
  int dim = initial_value1.rows();
  int seed = seed_;

  HmcSampler hmc1(dim, seed);
  HmcSampler hmc2(dim, seed);


  const int numlin = numlin_;
  const int numquad = numquad_;



  if (numlin >0){
    const Map<MatrixXd> F = as<Map<MatrixXd> > (F_);
    const Map<VectorXd> g = as<Map<VectorXd> > (g_);

    for(int i=0; i<numlin; i++){
      hmc1.addLinearConstraint(F.row(i),g(i));
      hmc2.addLinearConstraint(F.row(i),g(i));
    }
  }


  if (numquad >0){
    const Map<MatrixXd> Q = as<Map<MatrixXd> > (quadratics_);

    for(int i=0; i<numquad; i++){
      MatrixXd A = Q.block( i*(dim+2), 0, dim, dim);  //block(firstRow, firstCol, rows, cols)
      VectorXd B = Q.row(i*(dim+2) + dim);
      double C = Q(i*(dim+2) + dim+1,0);
      hmc1.addQuadraticConstraint(A,B,C);
      hmc2.addQuadraticConstraint(A,B,C);
    }

  }



  hmc1.setInitialValue(initial_value1);
  hmc2.setInitialValue(initial_value2);

  MatrixXd samples1(n+1,dim);
  MatrixXd samples2(n,dim);
  Map<VectorXd> velocity = as<Map<VectorXd> > (rnorm(dim));
  samples1.row(0) = hmc1.sampleNext_given_velocity(velocity);

  for (int i=0; i<n; i++){
    // std::cerr << "iteration" << i << std::endl;
    velocity = as<Map<VectorXd> > (rnorm(dim));
    samples1.row(i+1) = hmc1.sampleNext_given_velocity(velocity);
    // std::cerr << samples1.row(i+1) << std::endl;
    samples2.row(i) = hmc2.sampleNext_given_velocity(velocity);
    // std::cerr << samples2.row(i) << std::endl;
    // std::cerr << hmc1._verifyConstraints(samples1.row(i+1)) << std::endl;
    // std::cerr << hmc1._verifyConstraints(samples2.row(i)) << std::endl;
  }

    return List::create(Named("samples1") = wrap(samples1),
                        Named("samples2") = wrap(samples2));

}

