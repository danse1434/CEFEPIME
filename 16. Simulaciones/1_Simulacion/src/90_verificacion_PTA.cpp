#include <math.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]

//------------------------------------------------------------------------------#
//' Verificación de PTA
//' 
//' @param PRED objeto de tipo matriz con predicciones (columnas: ID, filas: 
//' observaciones). Este objeto está compuesto por concentraciones predichas. 
//' @param MIC vector con valores de MIC a probar
//' @param crit: criterio con el cual se declara que se ha alcanzado el objetivo 
//' farmacoterapéutico
//' 
//' @export
//' @return
//' 
// [[Rcpp::export()]]
List pta_verificador(arma::mat PRED, arma::vec MIC, double crit){
  // Tamaño del vector de MIC
  arma::uword n_mic = MIC.size();
  // Tamaño de la matriz de predicciones
  arma::uword n_pred_1 = PRED.n_rows;
  arma::uword n_pred_2 = PRED.n_cols;
  // Vector de PTA
  arma::fcolvec PTA_vec = arma::zeros<arma::fcolvec>(n_mic); //
  // Matriz con valores de fT>MIC vs MIC para cada individuo
  arma::fmat Pmat = arma::fmat(n_mic,n_pred_2,arma::fill::zeros);
  
  // Vector temporal de tamaño con número de columnas
  arma::fcolvec TmasMIC = arma::zeros<arma::fcolvec>(n_pred_2); //f
  arma::ucolvec TmasMIC_ind = arma::zeros<arma::ucolvec>(n_pred_2);
  
  /* Matriz temporal con el tamaño de "PRED"
   Se crea una matriz con enteros 1 o 0 */
  arma::umat M_ind = arma::umat(n_pred_1, n_pred_2, arma::fill::zeros);
  
  // Objetos de temporales de conteo
  float M_tot = 0;
  float TmasMIC_tot = 0;
  
  // Nivel de MIC
  for(arma::uword k = 0; k < n_mic; k++)
  {
    // Nivel individuos (ID_j)
    for(arma::uword j = 0; j < n_pred_2; j++)
      {
      // Nivel predicciones (t_i)
      for(arma::uword i = 0; i < n_pred_1; i++)
        {
        // Llena matriz con 1 si es mayor que MIC //
        if(PRED.at(i,j) >= MIC.at(k)){
          M_ind.at(i,j) = 1;
        } 
        } //END
      // Prop. instancias con valor a cero por individuo
      /* La suma o acumulación no sirven bien si la columna es entera, 
         se crea una variable interna con la suma. */
      M_tot = arma::accu( M_ind.col(j) );
      
      TmasMIC.at(j) = M_tot/n_pred_1;
      
      M_ind.zeros();
      
      if(TmasMIC.at(j) >= crit){
        TmasMIC_ind.at(j) = 1;
      }
      
      // Llenar matriz con valores de fT>MIC
      Pmat.at(k,j) = TmasMIC.at(j);
      
      } //END
    
    TmasMIC_tot = arma::accu(TmasMIC_ind);
    PTA_vec.at(k) =  TmasMIC_tot/ n_pred_2;
    TmasMIC_ind.zeros();
  } //END
  
  /* Valores de PTA */  
  DataFrame Y = DataFrame::create(
    _["MIC"] = MIC, 
    _["PTA"] = PTA_vec
  );
  
  /* Objetivo Farmacodinámico */  
  DataFrame P_df = DataFrame::create(
    _["MIC"] = MIC,
    _["ID"] = Pmat
  );
  
  /* Lista de Data.Frames Resultado */
  List Z = List::create(
    _["PTA"] = Y,
    _["fTmasMIC"] = P_df
  );
  
  return Z;
}

//------------------------------------------------------------------------------
//' Cálculo de Cmax, tmax y AUC (exposición)
//' 
//' @param data: es un data.frame de R (acepta tibble) que contiene datos de 
//' concentraciones predichas por individuo (cada columna es un ID), y en la 
//' primera columna tiene una referencia del tiempo al que corresponde cada 
//' observación (Cpred). 
//' @export
//' @return lista con dos tablas "Pmaximos" que contiene un resumen por 
//' individuo de la localización de Cmax y tmax, y "AUC_estm" que contiene 
//' un resumen por individuo del valor de AUC calculado sin tener en cuenta 
//' la fase terminal de cada curva.
//' 
//' @example
//' df1 <- res$ufCc %>% 
//' filter(group == 1) %>% 
//' pivot_wider( id_cols = time,
//'         names_from = id,
//'         names_prefix = 'ID', 
//'         values_from = ufCc )
//' UDF_exposure(df1)
//' 
// [[Rcpp::export()]]
List UDF_exposure(DataFrame data){
  
  arma::mat A = as<arma::mat>(internal::convert_using_rfunction(data, "as.matrix")); 
  
  arma::uword A_rows = A.n_rows; // No. observaciones
  arma::uword A_cols = A.n_cols; // No. individuos
  
  // Vector de tiempos
  arma::vec t_vec = A.col(0);
  // Matriz de predicciones
  arma::mat pred_vec = A.cols(1, A.n_cols-1);
  
  // Vector de concentraciones y tiempos máximos por individuo
  arma::vec tmax = arma::vec(pred_vec.n_cols);
  arma::vec cmax = arma::vec(pred_vec.n_cols);
  
  /* Localización de Cmax y tmax */
  arma::uword t; // Indice reemplazable
  for(arma::uword i = 0; i < pred_vec.n_cols; i++)
  {
    //std::cout << i << std::endl;
    t = pred_vec.col(i).index_max();
    tmax(i) = t_vec(t);
    cmax(i) = pred_vec(t,i);
  }
  
  /* Vectores de locación de resultados AUC */
  arma::fcolvec auc_vec = arma::zeros<arma::fcolvec>(A_rows);
  arma::fcolvec auc_id = arma::zeros<arma::fcolvec>(pred_vec.n_cols);
  
  // Cálculo de AUC
  for(arma::uword k = 0; k < pred_vec.n_cols; k++)
  {
    for(arma::uword j = 0; j < A_rows; j++)
    {
      if(j == 0){
        auc_vec(j) = 0;
      } else {
        // Regla del trapecio;
        auc_vec(j) = (t_vec(j) - t_vec(j-1)) * (pred_vec(j-1, k) + pred_vec(j, k))/2;
      }
    }
    auc_id(k) = sum(auc_vec);
  }
  
  // Vector de ID
  StringVector id_vec(pred_vec.n_cols);
  for(arma::uword k = 0; k < pred_vec.n_cols; k++){
    id_vec(k) = "ID" + std::to_string(k+1);
  }
  
  // Tabla Puntos máximos y mínimos en perfil
  DataFrame max_points =  DataFrame::create(
    _["ID"]   = id_vec,
    _["tmax"] = tmax,
    _["Cmax"] = cmax
  );
  
  // Tabla Resumen AUC
  DataFrame AUCvsID = DataFrame::create(
    _["ID"] = id_vec,
    _["AUC"] = auc_id
  );
  
  // Lista con tablas
  List Y = List::create(
    _["Pmaximos"] = max_points,
    _["AUC_estm"] = AUCvsID
  );
  
  return Y;
}

