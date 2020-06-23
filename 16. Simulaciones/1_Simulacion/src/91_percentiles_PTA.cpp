#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//-------------------------------------------------------------------------------
//' Función de obtención de percentiles
//' 
//' @param data (data.frame) la primera columna corresponde a valores de MIC, el 
//' resto corresponde al índice PK-PD para varios individuos.
//' @param prob_vec (vector) es un vector que indica las probabilidades para la 
//' obtención de cuantiles para cada MIC (fila).
//' 
//' @return (lista) contiene un data.frame, la primera columna son valores de MIC, 
//' y el resto son los percentiles correspondientes.
//' @export
//' 
//' @example
//' perc_fT(RESPTA1$fTmasMIC[[1]], c(0.1,0.2,0.3,0.4,0.5)) 
//' 
// [[Rcpp::export()]]
List perc_fT(DataFrame data, arma::frowvec prob_vec){
  arma::mat A = as<arma::mat>(internal::convert_using_rfunction(data, "as.matrix"));
  
  // "fT_mat" matrix con valores de índice PK-PD por MIC (filas), y percentiles 
  //(columnas)
  arma::mat fT_mat = A.cols(1, A.n_cols-1);
  // 
  // Creación de data.frame
  DataFrame x = DataFrame::create(
    _["MIC"] = A.col(0),
    _["p"]   = arma::quantile(fT_mat, prob_vec, 1)
  );
  // Creación de lista resultado
  List Y = List::create(
    _["x"] = x
  );
  
  return Y;
}

/*
perc_fT(RESPTA1$fTmasMIC[[1]], c(0.1,0.2,0.3,0.4,0.5))
*/
