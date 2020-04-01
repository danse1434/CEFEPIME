##------------------------------------------------------------------------#
## Nombre del Script: Estudio de tamizaje de covariables ------------------
##  
## Propósito del Script: realizar un estudio de la relación entre covariables 
##  con desviaciones eta con medios gráficos y regresión.
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion:  19-03-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
# Carga de paquetes
require(rlang)
require(tidyverse)
require(grid)

# Función de diseño layout
vplayout <- function(x, y) {
    viewport(layout.pos.row = x, layout.pos.col = y)
  }

# Selección de directorio principal
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '11. Screening_Covariables'))

##########################################################################-
# Introducción ------------------------------------------------------------
##########################################################################-
# Carga de archivo de datos original
data_ori <- read_delim("../10. Base_Refin/BASE/interv_censored.csv", 
             ";", escape_double = FALSE, trim_ws = TRUE, na = ".", 
             col_types = cols() )
# data_ori

##########################################################################-
# Apertura de archivos con parámetros individuales 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Apertura de archivo con parámetros estimados *indiv_param*
##  2 Apertura de archivo con parámetros estimados *indiv_eta*
##  3 Selección de las variables id y cualquiera que contenga "SAEM"
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

indiv_param <- read_csv(file.path('..', '10. Base_Refin', 'BASE', 'BASE_MODEL', 
                     'IndividualParameters', 'estimatedIndividualParameters.txt'), 
                     col_types = cols())

indiv_eta <- read_csv(file.path('..', '10. Base_Refin', 'BASE', 'BASE_MODEL', 
                                'IndividualParameters', 'estimatedRandomEffects.txt'),
                      col_types = cols())

indiv_param <- indiv_param %>% 
  select(id, matches("SAEM"))

indiv_eta <- indiv_eta %>% 
  select(id, matches("SAEM"))

##########################################################################-
# Procesamiento de tabla de datos -----------------------------------------
##########################################################################-
# Unión de tablas
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Unir tabla indiv_param a la tabla data_ori
##  2 Unir tabla indiv_eta a la tabla data_ori
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data_ori <- data_ori %>% 
  left_join(indiv_param, by = c('ID' = 'id')) %>% 
  left_join(indiv_eta, by = c('ID' = 'id'))

##########################################################################-
# Procesamiento de tabla
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar de tabla completa
##  2 Seleccionar de eventos de administración y que tenga en tiempo en cero, 
##  esto deja una tabla con las covariables y el valor de los parámetros 
##  (etas) individuales estimados.
##  3 Eliminar variables innecesarias
##  4 Eliminar variables que empiezen por valor de parámetro (par. individ)
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data <- data_ori %>%
  filter(EVID == '1' & TIME == 0) %>%
  select(-DV, -MDV, -EVID, -AMT, -RATE, -ADDL, -II, -SS, -TIME, -YTYPE,
         -LIMIT, -ADM, -CENS) %>% 
  select(-matches("^Cl\\_|^V1|^Q|^V2"))


##########################################################################-
# Creación de gráficos de manera programática -----------------------------
##########################################################################-
#' Función de correlación
#' 
#' Crear quosures con las variables \code{x} y \code{y}
#' Extraer las variables seleccionadas por x y y como vectores atómicos
#' Crear un carácter con la correlación a 3 dígitos
#'
#' @param df Tabla de datos con variables
#' @param x Variable dependiente (covariable)
#' @param y Variable independiente (desviación eta)
#'
#' @return Correlación entre variables \code{x} y \code{y} dentro de la 
#' tabla especificada
#' @export
#' @examples
#' corr_eq(data, AGEA, eta_Cl_SAEM)
corr_eq <- function(df, x, y) {
  qx <- rlang::enquo(x); qy <- rlang::enquo(y)

  x <- pull(df, !!qx); y <- pull(df, !!qy)

  paste0("r = ", round(cor(x, y), digits = 3)) %>% return(.)
}

##########################################################################-
#' Función de creación de gráficos *et_pf*
#' 
#' Se crea un gráfico con inclusión de valor de correlación entre desviaciones 
#' eta, y covariables continuas.
#' 
#' @param x Variable dependiente (covariable)
#' @param y Variable independiente (desviación eta)
#' @param xlab Etiqueta de eje X
#' @param ylab Etiqueta de eje Y 
#' @param col Color de línea y relleno; azul (predeterminado).
#' @param type Tipo de gráfico: 1 dispersión (predeterminado), 2 diagrama 
#' de cajas.
#' @param ymin Valor mínimo del eje Y
#' @param ymax Valor máximo del eje Y
#'
#' @return Gráfico con correlación entre covariable (\code{x}) y desviación 
#' eta (\code{y}), este puede ser dispersión o boxplot. 
#' @export
#' @examples
#' et_pf(AGEA, eta_Cl_SAEM, 'Edad', bquote(eta[Cl]), type = 1)
#' 
et_pf <- function(x, y, xlab, ylab, col = "blue", type = 1, 
                         ymin = -1.5, ymax = 1.5) {
  
  x <- rlang::ensym(x); y <- rlang::ensym(y)
  
  if (type == 1) {
    c <- pull(data, !!x)
    xpos <- min(c) + (0.8 * (max(c) - min(c)))
    lab_df <- tribble(~ x, ~ y, xpos, 1)
    
    G1 <- data %>% ggplot(aes(x = !!x, y = !!y)) +
      geom_hline(yintercept = 0, lty = 'dashed') +
      geom_point(shape = 22, fill = col, colour = 'black', alpha = 0.5) + 
      stat_smooth(formula = 'y ~ x', method = 'lm', 
                  se = TRUE, colour = col, fill = col, alpha = 0.1) +
      coord_cartesian(ylim = c(ymin, ymax)) +
      theme_bw() + theme(panel.grid = element_line(colour = NA)) +
      ylab(ylab) + xlab(xlab) + 
      geom_text(data = lab_df, aes(x = x, y = y, label = corr_eq(data, !!x, !!y))
                )
    return(G1)
  } else if (type == 2) {
    G2 <- data %>% ggplot(aes(x = factor(!!x), y = !!y)) +
      geom_hline(yintercept = 0, lty = 'dashed') +
      geom_boxplot(fill = col, colour = 'black', alpha = 0.5) +
      coord_cartesian(ylim = c(ymin, ymax)) +
      theme_bw() + theme(panel.grid = element_line(colour = NA)) +
      ylab(ylab) + xlab(xlab)
    return(G2)
  }
}

##########################################################################-
#' Creación de lista con gráficos para cada variable en el set de datos
#'
#' @param ETA Desviación eta (parámetros individuales) de interés
#' @param x Covariable (variable dependiente)
#' @param col Tema de color de los gráficos
#'
#' @return Objeto tipo lista con gráficos para cada covariable de interés
#' @export
#' @examples list_object(eta_Cl_SAEM, 'Cl', col = 'blue')
#' 
list_object <- function(ETA, x, col) {
  # Función de parada
  stopifnot(is.character(x), is.character(col))
  
  ETA = rlang::ensym(ETA)
  ls <- list()
  yexpr = bquote(eta[.(x)])
  ls[[1]] <- et_pf(SEXF, !!ETA, 'Sexo', yexpr, type = 2, col = col)
  ls[[2]] <- et_pf(AGEA, !!ETA, 'Edad (años)', yexpr, col = col)
  ls[[3]] <- et_pf(WTKG, !!ETA, 'Peso (kg)', yexpr, col = col)
  ls[[4]] <- et_pf(HCM, !!ETA, 'Talla (cm)', yexpr, col = col)
  ls[[5]] <-
    et_pf(IMCKGM2, !!ETA, expression(IMC ~ (kg / cm ^ 2)), yexpr, col = col)
  ls[[6]] <-
    et_pf(SCRMGDL, !!ETA, expression(S[CR] ~ (mg / dL)), yexpr, col = col)
  ls[[7]] <-
    et_pf(CLCRMLMIN, !!ETA, expression(CrCl ~ (mL / min / 1.73 ~ m ^ 2)), 
          yexpr, col = col)
  ls[[8]] <- et_pf(PROGDL, !!ETA, 'Proteínas (g/dL)', yexpr, col = col)
  ls[[9]] <- et_pf(ALBGDL, !!ETA, 'Albúmina (g/dL)', yexpr, col = col)
  ls[[10]] <- et_pf(DIND, !!ETA, 'DIND', yexpr, col = col)
  ls[[11]] <-
    et_pf(RAN, !!ETA, expression(RAN ~ (mm ^ -3)), yexpr, col = col)
  ls[[12]] <- et_pf(RAL, !!ETA, 'RAL (/mm^3)', yexpr, col = col)
  ls[[13]] <- 
    et_pf(ANTU, !!ETA, 'Uso de antibiótico', yexpr, type = 2, col = col)
  ls[[14]] <- et_pf(Dx, !!ETA, 'Tipo Leucemia', yexpr, type = 2, col = col)
  
  return(ls)
}

##########################################################################-
# Creación de listas para cada tipo de desviación eta con colores difentes
Cl_ls <- list_object(eta_Cl_SAEM, 'Cl', col = 'blue')
V1_ls <- list_object(eta_V1_SAEM, 'V1', col = 'red')
Q_ls <- list_object(eta_V1_SAEM, 'Q', col = 'green')
V2_ls <- list_object(eta_V2_SAEM, 'V2', col = 'purple')

##########################################################################-
# Creación de documentos con gráficas en forma de layout

pdf('Figuras/Cl_correlation.pdf', width = 6*1.7, height = 4*1.7);{
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 4)))
  print(Cl_ls[[1]], vp = vplayout(1, 1))
  print(Cl_ls[[2]], vp = vplayout(1, 2))
  print(Cl_ls[[3]], vp = vplayout(1, 3))
  print(Cl_ls[[4]], vp = vplayout(1, 4))
  
  print(Cl_ls[[5]], vp = vplayout(2, 1))
  print(Cl_ls[[6]], vp = vplayout(2, 2))
  print(Cl_ls[[7]], vp = vplayout(2, 3))
  print(Cl_ls[[8]], vp = vplayout(2, 4))
  
  print(Cl_ls[[9]], vp = vplayout(3, 1))
  print(Cl_ls[[10]], vp = vplayout(3, 2))
  print(Cl_ls[[13]], vp = vplayout(3, 3))
  print(Cl_ls[[14]], vp = vplayout(3, 4))
};dev.off()

pdf('Figuras/V1_correlation.pdf', width = 6*1.7, height = 4*1.7);{
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 4)))
  print(V1_ls[[1]], vp = vplayout(1, 1))
  print(V1_ls[[2]], vp = vplayout(1, 2))
  print(V1_ls[[3]], vp = vplayout(1, 3))
  print(V1_ls[[4]], vp = vplayout(1, 4))
  
  print(V1_ls[[5]], vp = vplayout(2, 1))
  print(V1_ls[[6]], vp = vplayout(2, 2))
  print(V1_ls[[7]], vp = vplayout(2, 3))
  print(V1_ls[[8]], vp = vplayout(2, 4))
  
  print(V1_ls[[9]], vp = vplayout(3, 1))
  print(V1_ls[[10]], vp = vplayout(3, 2))
  print(V1_ls[[13]], vp = vplayout(3, 3))
  print(V1_ls[[14]], vp = vplayout(3, 4))
};dev.off()

pdf('Figuras/Q_correlation.pdf', width = 6*1.7, height = 4*1.7);{
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 4)))
  print(Q_ls[[1]], vp = vplayout(1, 1))
  print(Q_ls[[2]], vp = vplayout(1, 2))
  print(Q_ls[[3]], vp = vplayout(1, 3))
  print(Q_ls[[4]], vp = vplayout(1, 4))
  
  print(Q_ls[[5]], vp = vplayout(2, 1))
  print(Q_ls[[6]], vp = vplayout(2, 2))
  print(Q_ls[[7]], vp = vplayout(2, 3))
  print(Q_ls[[8]], vp = vplayout(2, 4))
  
  print(Q_ls[[9]], vp = vplayout(3, 1))
  print(Q_ls[[10]], vp = vplayout(3, 2))
  print(Q_ls[[13]], vp = vplayout(3, 3))
  print(Q_ls[[14]], vp = vplayout(3, 4))
};dev.off()

pdf('Figuras/V2_correlation.pdf', width = 6*1.7, height = 4*1.7);{
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 4)))
  print(V2_ls[[1]], vp = vplayout(1, 1))
  print(V2_ls[[2]], vp = vplayout(1, 2))
  print(V2_ls[[3]], vp = vplayout(1, 3))
  print(V2_ls[[4]], vp = vplayout(1, 4))
  
  print(V2_ls[[5]], vp = vplayout(2, 1))
  print(V2_ls[[6]], vp = vplayout(2, 2))
  print(V2_ls[[7]], vp = vplayout(2, 3))
  print(V2_ls[[8]], vp = vplayout(2, 4))
  
  print(V2_ls[[9]], vp = vplayout(3, 1))
  print(V2_ls[[10]], vp = vplayout(3, 2))
  print(V2_ls[[13]], vp = vplayout(3, 3))
  print(V2_ls[[14]], vp = vplayout(3, 4))
};dev.off()

