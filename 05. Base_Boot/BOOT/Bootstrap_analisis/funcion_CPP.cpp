#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export()]]
NumericVector func1(int l){
  NumericVector y(abs(l - 1) + 1);
  std::iota(y.begin(), y.end(), 1);
  
  return y;
  
}


// [[Rcpp::export()]]
arma::vec func2(int l){
  NumericVector w(abs(l - 0) + 1);
  std::iota(w.begin(), w.end(), 1);
  arma::vec x = as<arma::vec>(w);
  
  return x;
  
}


// [[Rcpp::export()]]
arma::mat func3(arma::vec A, int j, int l){
  arma::mat y(l+1, 2);
  
  y.col(0).ones();
  
  for(int i = 0; i < (l+1); i++){
    y(i,1) = A[j - l - 1 + i];
  }

  return y;
}



// [[Rcpp::export()]]
arma::mat func4(arma::vec A, int j, int l){
  arma::vec w(l+1);
  std::iota(w.begin(), w.end(), 1);
  
  arma::mat y(l+1, 2);
  
  y.col(0).ones();
  
  for(int i = 0; i < (l+1); i++){
    y(i,1) = w[i];
  }
  
  return y;
}

// [[Rcpp::export()]]
arma::vec func5(arma::vec A, int j, int l){
  arma::vec y(l+1);
  
  for(int i = 0; i < (l+1); i++){
    y[i] = A[j - l - 1 + i];
  }
  
  return y;
}




// [[Rcpp::export()]]
arma::vec func6(arma::mat x, arma::vec y){
  arma::vec theta = inv(x.t()*x)*x.t()*y;
  return theta;
}

  
// [[Rcpp::export()]]
arma::vec func7(arma::mat x, arma::vec y){
  arma::vec theta = inv(x.t()*x)*x.t()*y;
  
  int n = x.n_rows;
  int p = theta.n_rows;
  
  arma::vec D = y - (x*theta);
  arma::vec M = (D.t() * D)/(n-p);
  
  double MSE = M(0,0);
  
  arma::vec S = arma::sqrt(MSE*arma::diagvec(inv(x.t()*x)) );
  
  return S;
}


// [[Rcpp::export()]]
arma::vec func8(arma::mat x, arma::vec y){
  arma::vec theta = inv(x.t()*x)*x.t()*y;
  
  int n = x.n_rows;
  int p = theta.n_rows;
  
  arma::vec D = y - (x*theta);
  arma::vec M = (D.t() * D)/(n-p);
  
  double MSE = M(0,0);
  
  arma::vec S = arma::sqrt(MSE*arma::diagvec(inv(x.t()*x)) );
  
  return abs(theta)/S;
}


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



// [[Rcpp::export()]]
double likelihood(arma::vec yobs, arma::vec ypred) {
  arma::vec e = yobs - ypred;
  arma::vec r = pow(e,2);
  
  int n = r.size();
  double S;
  for(int i = 0; i < n; i++){
    S += r[i];
  }
  
  double var;
  for(int i = 0; i < n; i++){
    var += r[i]/n;
  }
  
  double L = -(log(2*M_PI)*n/2) - (log(var)*n/2) - (S/(2*var));
  
  return(L);
}


