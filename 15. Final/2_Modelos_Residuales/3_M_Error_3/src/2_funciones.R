#-------------------------------------------------------------------------------#
#' Función de gráficos diagnósticos de bondad de ajuste
#' Esta función asume que existe una variable de tipo 
#' 
#' @param x variable en eje X
#' @param y variable en eje Y
#' @param xspline variable de eje X en tabla de spline
#' @param yspline variable de eje Y en tabla de spline 
#' @param xconfint variable de eje X en tabla de IC
#' @param yconfint_lo variable de eje Y inferior en tabla de IC
#' @param yconfint_up variable de eje Y superior en tabla de IC
#' @param colourp color de gráfico diagnóstico
#' @param xlab etiqueta de eje X
#' @param ylab etiqueta de eje Y
#'
#' @return
#' @export este se basa en un objeto llamado "y_1_obsVsPred", y 
#' otro llamado "y_1_visualGuides"
#'
#' @examples
#'   GOF_PRED(x = popPred, y = y_1, xspline = popPred_spline_abscissa, 
#'   yspline = popPred_spline, xconfint = popPred_ci_abscissa, 
#'   yconfint_lo = popPred_piLower, yconfint_up = popPred_piUpper, 
#'   colourp = 'blue1', xlab = 'PRED', ylab = 'OBS')
#' 
GOF_PRED <- function(x, y, xspline, yspline, xconfint, yconfint_lo, 
                     yconfint_up, colourp, xlab, ylab) {
  x           = rlang::ensym(x)
  y           = rlang::ensym(y)
  xspline     = rlang::ensym(xspline)
  yspline     = rlang::ensym(yspline)
  xconfint    = rlang::ensym(xconfint)
  yconfint_lo = rlang::ensym(yconfint_lo)
  yconfint_up = rlang::ensym(yconfint_up)
  
  y_1_obsVsPred %>% 	
    ggplot(mapping = aes(x = !!x, y = !!y, group = ID)) +	
    geom_point(shape = 1) + 	
    xlab(xlab) + ylab(ylab) +	
    geom_abline(slope = 1, intercept = 0, lty = 'dotted') + 	
    geom_line(y_1_visualGuides, 	
              mapping =  aes(x = !!xspline, y = !!yspline), 	
              inherit.aes = F, colour = colourp) +	
    geom_ribbon(y_1_visualGuides, 	
                mapping =  aes(x = !!xconfint, 	
                               ymin = !!yconfint_lo,	
                               ymax = !!yconfint_up), 	
                inherit.aes = F, fill = colourp, alpha = 0.1) +	
    coord_cartesian(xlim = c(0,90), ylim = c(0,90))	%>% 
    return(.)
}

#-------------------------------------------------------------------------------#
#' Función de gráficos de residuales vs concentración
#'
#' @param x variable en eje X
#' @param y variable en eje Y
#' @param xspline variable de eje X en tabla de spline
#' @param yspline variable de eje Y en tabla de spline
#' @param perc_data tabla de percentiles
#' @param xlab etiqueta de eje X
#' @param ylab etiqueta de eje Y
#'
#' @return se obtiene un gráfico de residuales con datos de *y_1_residuals* y 
#' *y_1_spline*
#' @export
#'
#' @examples 
#' RES_TSFD(x = time, y = pwRes, xspline = time_pwRes, yspline = time_pwRes_spline,
#' perc_data = y_1_time_percentiles_pwRes, xlab = 'TSFD', ylab = 'PWRES')
#' 
RES_TSFD <- function(x, y, xspline, yspline, perc_data, xlab, ylab) {
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  xspline <- rlang::ensym(xspline)
  yspline <- rlang::ensym(yspline)
  stopifnot(is.data.frame(perc_data) == TRUE)
  
  y_1_residuals %>% 
    ggplot(mapping = aes(x = !!x, y = !!y)) +
    geom_hline(yintercept = 0) +
    geom_point(col = '#4682B4') +
    geom_line(y_1_spline, 
              mapping = aes(x = !!xspline, y = !!yspline), 
              col = '#EBA213', lty = 'solid', size = 1) +
    geom_line(perc_data, 
              mapping = aes(x = !!x, y = empirical_median), 
              col = '#EBA213', lty = 'dashed', size = 1.2) + 
    geom_line(perc_data, 
              mapping = aes(x = !!x, y = empirical_lower), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(perc_data, 
              mapping = aes(x = !!x, y = empirical_upper), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    xlab(xlab) + ylab(ylab) %>% return(.)
}

#-------------------------------------------------------------------------------#
#' Función de gráficos de residuales vs tiempo
#'
#' @param x variable en eje X
#' @param y variable en eje Y
#' @param xlab etiqueta de eje X
#' @param ylab etiqueta de eje Y
#'
#' @return
#' @export
#'
#' @examples
#' RES_TAD(x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')
#' 
RES_TAD <- function(x, y, xlab, ylab) {
  #-------------------------------------------------------------------------------#
  # Volver las variables x y y en expresiones para evaluacián tardáa
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  
  #................................................................................
  #  Calculo de percentiles empíricos
  #  1 Tomar la variable de datos a explorar.
  #  2 Crear una variable que tome 6 grupos a partir de los datos ordenados 
  #  con la función dplyr::ntile().
  #  3 Agrupar por la variable recién creado.
  #  4 Resumir por la media de TAD, y los percentiles P5%, P50%, y P95% de 
  #  los datos de residuales.
  #................................................................................

  #-------------------------------------------------------------------------------#
  A <- y_1_residuals %>%
    mutate(gr = ntile(TAD, 6)) %>% 
    group_by(gr) %>% 
    summarise(TIME = mean(TAD), 
              ME = quantile(x = !!y, probs = 0.50), 
              LI = quantile(x = !!y, probs = 0.05), 
              LS = quantile(x = !!y, probs = 0.95))
  
  #-------------------------------------------------------------------------------#
  # Crear el gráfico de residuales con las especificaciones mostradas
  y_1_residuals %>% 
    ggplot(mapping = aes(x = !!x, y = !!y)) +
    geom_hline(yintercept = 0) +
    geom_point(col = '#4682B4') +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    xlab(xlab) + ylab(ylab) +
    stat_smooth(method = 'loess', formula = y ~ x, se = FALSE, 
                col = '#EBA213', lty = 'solid', size = 1) +
    geom_line(data = A, aes(TIME, ME), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(data = A, aes(TIME, LI), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(data = A, aes(TIME, LS), 
              col = '#EBA213', lty = 'dashed', size = 1.2) + 
    geom_text_repel(data = filter(y_1_residuals, abs(!!y) > 2), 
                    aes(label = ID))  %>% 
    return(.)
}


#-------------------------------------------------------------------------------#
#' Creación de gráfico de residuales
#'
#' @param x Variable X (en *y_1_residuals*) para puntos
#' @param y Variable Y (en *y_1_residuals*) para puntos
#' @param xspline Variable X (en *y_1_spline*) para spline
#' @param yspline Variable Y (en *y_1_spline*) para spline
#' @param perc_data Tabla con datos de percentiles empíricos
#' @param xlab Etiqueta de eje X personalizada
#' @param ylab Etiqueta de eje Y personalizada
#'
#' @return Gráfico de residuales vs predicciones
#' @export
#'
#' @examples
#' RES_PRE(x = prediction_pwRes, y = pwRes, xspline = prediction_pwRes, 
#' yspline = prediction_pwRes_spline, perc_data = y_1_prediction_percentiles_pwRes,
#' xlab = 'PRED', ylab = 'WRES')
#' 
RES_PRE <- function(x, y, xspline, yspline, perc_data, xlab, ylab) {
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  xspline <- rlang::ensym(xspline)
  yspline <- rlang::ensym(yspline)
  stopifnot(is.data.frame(perc_data) == TRUE)
  
  y_1_residuals %>% 
    ggplot(mapping = aes(x = !!x, y = !!y)) +
    geom_hline(yintercept = 0) +
    geom_point(col = '#8721B8') +
    geom_line(y_1_spline, 
              mapping = aes(x = !!xspline, y = !!yspline), 
              col = '#EBA213', lty = 'solid', size = 1) +
    geom_line(data = perc_data,
              mapping = aes(x = prediction, y = empirical_median),
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(data = perc_data,
              mapping = aes(x = prediction, y = empirical_lower),
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(data = perc_data,
              mapping = aes(x = prediction, y = empirical_upper),
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    # Texto de etiquetas
    geom_text_repel(data = filter(y_1_residuals, abs(!!y) > 2), 
                    aes(label = ID)) + 
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    xlab(xlab) + ylab(ylab) %>% return(.)
}
