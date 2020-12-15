#-------------------------------------------------------------------------------#
#' Mean Absolute Error (MAE)
#' Son los valores absolutos de cada residual y luego se calcula la media de 
#' los mismos.
#' 
#' @param predictions numeric vector with model predictions
#' @param true_values numeric vector with original values
#'
#' @return (double) summary statistic
#' @export
#'
#' @examples
#' set.seed(3)
#' x = rnorm(1000, 100, 10)
#' f = (x * 3.14) + 1.56
#' y = sapply(f, function(x) rnorm(1,x,3))
#' mod <- lm(y ~ x)
#' ypred <- predict(mod)
#' print( MAE(ypred, y) )
#'  
MAE <- function(predictions, true_values) {
  stopifnot(length(predictions) == length(true_values))
  
  x = mapply(function(x,y){abs(x-y)}, predictions, true_values)
  return(mean(x, na.rm = TRUE))
}

#-------------------------------------------------------------------------------#
#' Root Mean Squared Error (RMSE)
#' Son los cuadrados de cada residual que se promedian y luego se obtiene la 
#' raíz cuadrada.
#'
#' @param predictions vector numérico con predicciones de un modelo
#' @param true_values vector numérico con valores originales 
#'
#' @return (double) Estadístico de resumen.
#' @export
#'
#' @examples
#' set.seed(3)
#' x = rnorm(1000, 100, 10)
#' f = (x * 3.14) + 1.56
#' y = sapply(f, function(x) rnorm(1,x,3))
#' mod <- lm(y ~ x)
#' ypred <- predict(mod)
#' print( RMSE(ypred, y) )
#' 
RMSE <- function(predictions, true_values) {
  stopifnot(length(predictions) == length(true_values))
  
  x = mapply(function(x,y){ (x-y)^2 }, predictions, true_values)
  return(sqrt(mean(x, na.rm = TRUE)))
}

#-------------------------------------------------------------------------------#
#' Mean Absolute Percentage Error (MAPE)
#' Es el valor absoluto de cada residual dividido por el valor actual de 
#' ese punto para obtener un porcentaje, luego promediar.
#'
#' @param predictions vector numérico con predicciones de un modelo
#' @param true_values vector numérico con valores originales 
#'
#' @return (double) Estadístico de resumen.
#' @export
#'
#' @examples
#' set.seed(3)
#' x = rnorm(1000, 100, 10)
#' f = (x * 3.14) + 1.56
#' y = sapply(f, function(x) rnorm(1,x,3))
#' mod <- lm(y ~ x)
#' ypred <- predict(mod)
#' print( MAPE(ypred, y) )
#' 
MAPE <- function(predictions, true_values) {
  stopifnot(length(predictions) == length(true_values))
  
  x <- mapply(function(x,y) {abs((x-y)/y)}, predictions, true_values)
  return(mean(x, na.rm = TRUE))
}

#-------------------------------------------------------------------------------#
#' Convertir funciones de rendimiento predictivo en funciones para Boot
#' Esta función permite obtener 
#' 
#'
#' @param data datos
#' @param predictions vector numérico con predicciones de un modelo
#' @param true_values vector numérico con valores originales 
#' @param fun función de rendimiento predictivo (MAE, RMSE o MAPE)
#' @param indices indices
#'
#' @return
#' @export
#'
#' @examples
#' set.seed(3)
#' x = rnorm(1000, 100, 10)
#' f = (x * 3.14) + 1.56
#' y = sapply(f, function(x) rnorm(1,x,3))
#' mod <- lm(y ~ x)
#' ypred <- predict(mod)
#' 
#' df=data.frame(y,ypred)
#' 
#' print( bootPredPerf(df, ypred, y, MAPE, 1:30) )
#' 
bootPredPerf <- function(data, predictions, true_values, fun, indices) {
  ypred = data[indices,predictions]
  yobs  = data[indices,true_values]
  return( fun(ypred, yobs) )
}


#-------------------------------------------------------------------------------#
#' Cálculo de intervalos de confianza por bootstrap no paramétrico 
#'
#' @param data datos
#' @param predictions vector numérico con predicciones de un modelo
#' @param true_values vector numérico con valores originales  
#' @param fun función de rendimiento predictivo (MAE, RMSE, MAPE)
#' @param R número de replicados
#' @param stats estadístico para conversión a bootstrap
#' @param type tipo de intervalo percentil (perc), normal (norm), básico (basic), 
#'  estudentizado (stud), bca (bca).
#'
#' @return lista con valor original (t0), límite inferior (li) y superior (ls) 
#'  de intervalo de confianza del 95%.
#' @export
#'
#' @examples
#' set.seed(3)
#' x = rnorm(1000, 100, 10)
#' f = (x * 3.14) + 1.56
#' y = sapply(f, function(x) rnorm(1,x,3))
#' mod <- lm(y ~ x)
#' ypred <- predict(mod)
#' 
#' df=data.frame(y,ypred)
#' 
#' # Retorna la lista con valores del intervalo
#' bootICPredPerf(df, ypred, y, MAPE) 
#' 
#' 
bootICPredPerf <- function(data, predictions, true_values, fun, 
                           R=1e3, stats=bootPredPerf, type='perc') {
  b = boot::boot(data=data, statistic=stats, R = R, 
           predictions=predictions, true_values=true_values, fun=fun)
  b1 = boot::boot.ci(b, type = type)
  
  return(list('t0'=b1$t0, 'li'=b1$percent[4], 'ls'=b1$percent[5]))
}

#-------------------------------------------------------------------------------#
#' Etiqueta con valores de rendimiento predictivo
#'
#' @param data dataFrame que contiene las predicciones y valores verdaderos
#' @param predictions vector numérico con predicciones de un modelo 
#' @param true_values vector numérico con valores originales  
#' @param x posición relativa en x de caja (empieza aquí a la derecha)
#' @param y posición relativa en y de caja (empieza aquí hacía arriba)
#' @param xlim límite verdadero del eje x en el gráfico
#' @param ylim límite verdadero del eje y en el gráfico
#' @param round para redondear las cifras
#' @param size tamaño de la letra
#' @param boot indica si se deben calcular los intervalos
#' @param R número de replicaciones
#' @param stats convertir a función sensible a índice
#' @param type tipo de percentil (default='perc')
#'
#' @return objeto de tipo ggplot con la etiqueta
#' @export
#'
#' @examples
#' require(ggplot2)
#' set.seed(3)
#' x = rnorm(1000, 100, 10)
#' f = (x * 3.14) + 1.56
#' y = sapply(f, function(x) rnorm(1,x,3))
#' mod <- lm(y ~ x)
#' ypred <- predict(mod)
#' 
#' df=data.frame(x,y)
#' 
#' ggplot(df, aes(x=x, y=y)) +
#' predictivePerformaceLabel(df, 'x','y')
#' 
#' 
predictivePerformaceLabel <- function(data, predictions, true_values, 
                                      x = 0.80, y = 0.2, xlim = c(0, 0), ylim = c(0, 0), 
                                      round = 3, size=3.5, 
                                      boot=FALSE, R=1e1, 
                                      stats=bootPredPerf, type='perc') {
  # Datos de predicción
  ypred = data[,predictions][predictions]
  yobs  = data[true_values][true_values]
  
  # Matriz con dimensiones en X y Y
  mat   = matrix(c(min(ypred), max(ypred), min(yobs), max(yobs)), 2, byrow=T)
  
  # Posicionamiento de etiqueta
  if (xlim[1]!=0 | xlim[2]!=0) {
    xpos = (xlim[2] - xlim[1]) * x + xlim[1]
  } else {
    xpos = (mat[1, 2] - mat[1, 1]) * x + mat[1, 1] 
  }
  
  if (ylim[1]!=0 | ylim[2]!=0) {
    ypos = (ylim[2] - ylim[1]) * y + ylim[1]
  } else {
    ypos = max((mat[2, 2] - mat[2, 1]) * y + mat[2, 1], y*ylim)
  }
  
  # Etiqueta con estadísticos de rendimiento predictivo
  if (boot) {
    MAE_estim = bootICPredPerf(data, predictions, true_values, MAE, R=R, stats=stats, type=type)
    MAE_estim = lapply(MAE_estim, round, c('digits' = round))
    
    RMSE_estim = bootICPredPerf(data, predictions, true_values, RMSE, R=R, stats=stats, type=type)
    RMSE_estim = lapply(RMSE_estim, round, c('digits' = round))
    
    MAPE_estim = bootICPredPerf(data, predictions, true_values, MAPE, R=R, stats=stats, type=type)
    MAPE_estim = lapply(MAPE_estim, round, c('digits' = round))
    
    m <- glue::glue(
      'MAE = {MAE_estim$t0} ({MAE_estim$li}-{MAE_estim$ls})\n',
      'RMSE = {RMSE_estim$t0} ({RMSE_estim$li}-{RMSE_estim$ls})\n',
      'MAPE = {MAPE_estim$t0*100} ({MAPE_estim$li*100}-{MAPE_estim$ls*100})%'
      )
    
  } else {
    m <- paste0(
      'MAE   = ',  round(MAE(ypred, yobs), round), '\n',
      'RMSE = ', round(RMSE(ypred, yobs), round),'\n',
      'MAPE = ', round(MAPE(ypred, yobs), round)*100, '%'
    )
  }
  # Objeto de tipo
  object = list(
    ggplot2::geom_label(mapping=aes(x=xpos, y=ypos, label=m), hjust=0, vjust=0, size=size)
  )
  
  return(object)
}