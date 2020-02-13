#' Retropredicción de concentraciones desde modelo lineal
#'
#' @param model # Un modelo de regresión en formato de objeto clase 'lm'
#' @param data # Una columna con datos de la variable dependiente
#' @param ... 
#'
#' @return # lista con resultados de retro-predicción
#' @export
#'
#' @examples
retro.predic = function(model, data, weights,...){
    if(!("chemCal" %in% (.packages()))){
      library("chemCal", lib.loc="C:/Program Files/R/R-3.6.0/library")}
    
  model = model; 
  D = data;   
  N = dim(D)[1]
  p1 = c(rep(NA, length = N));   p2 = c(rep(NA, length = N)); p3 = c(rep(NA, length = N))
  p4 = c(rep(NA, length = N));   p5 = c(rep(NA, length = N))
  
  if (missing(weights)) {
    W = model$weights
    if (is.double(W) == TRUE) {
      for (i in 1:N) {
        p1[i] = inverse.predict(object=model, newdata=D[[i,'m_Diam']], alpha=0.05, ws = W[i])$Prediction
        p2[i] = inverse.predict(object=model, newdata=D[[i,'m_Diam']], alpha=0.05, ws = W[i])$`Standard Error`
        p3[i] = inverse.predict(object=model, newdata=D[[i,'m_Diam']], alpha=0.05, ws = W[i])$Confidence
        p4[i] = p1[i] - p3[i]
        p5[i] = p1[i] + p3[i]}
    }
    
    else {
      for (i in 1:N) {
        p1[i] = inverse.predict(object=model, newdata=D[[i,'m_Diam']], alpha=0.05)$Prediction
        p2[i] = inverse.predict(object=model, newdata=D[[i,'m_Diam']], alpha=0.05)$`Standard Error`
        p3[i] = inverse.predict(object=model, newdata=D[[i,'m_Diam']], alpha=0.05)$Confidence
        p4[i] = p1[i] - p3[i]
        p5[i] = p1[i] + p3[i]}
    }
    
    
  } else {
    W = weights
    for (i in 1:N) {
        p1[i] = inverse.predict(object=model, newdata=D[[i,'m_Diam']], alpha=0.05, ws = W[i])$Prediction
        p2[i] = inverse.predict(object=model, newdata=D[[i,'m_Diam']], alpha=0.05, ws = W[i])$`Standard Error`
        p3[i] = inverse.predict(object=model, newdata=D[[i,'m_Diam']], alpha=0.05, ws = W[i])$Confidence
        p4[i] = p1[i] - p3[i]
        p5[i] = p1[i] + p3[i]} 
  }
    
    #####################################################################################################-
    # Predicción inversa de resultados
    pG1 = data.frame(p1,p2,p3,p4,p5)
    colnames(pG1) = c("Pred.","SE","Conf.","LI","LS")
    # Predicción en dominio normal
    pG2 = exp(pG1)
    # colnames(pG2) = c("Pred.","SE","Conf.","LI","LS")
    # Predicción en dominio normal
    # pG3 = data.frame(s=seq(1,N,1),(pG2-(D$C_pCre))/D$C_pCre)[,c(1,2,5,6)]
    # pG3
return(list(pG1,pG2))}

# retro.predic = Vectorize(FUN = retro.predic, vectorize.args = c('weights'))

