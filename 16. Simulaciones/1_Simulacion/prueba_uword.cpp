#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export()]]

double func1(int int1, int int2){
  double ratio = ((double) int1)/int2;
  
  arma::uword word1 = 666;
  double word2 = ratio * word1;
  
  return word2;
  
}


/***R
func1(1,2)
*/