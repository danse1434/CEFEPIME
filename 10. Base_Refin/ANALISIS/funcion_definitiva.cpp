#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]
arma::vec slop_signC(arma::vec A, int j, int l){
  arma::vec w(l+1);
  std::iota(w.begin(), w.end(), 1);
  
  arma::mat x(l+1,2);
  x.col(0).ones();
  
  for(int i = 0; i < (l+1); i++){
    x(i,1) = w[i];
  }
  
  arma::vec y(l+1);
  
  for(int i = 0; i < (l+1); i++){
    y[i] = A[j - l - 1 + i];
  }
  
  
  arma::vec theta = inv(x.t()*x)*x.t()*y;
  
  int n = x.n_rows;
  int p = theta.n_rows;
  
  arma::vec D = y - (x*theta);
  arma::vec M = (D.t() * D)/(n-p);
  
  double MSE = M(0,0);
  
  arma::vec S = arma::sqrt(MSE*arma::diagvec(inv(x.t()*x)) );
  
  return abs(theta)/S;
}



