##=========================================================================#
## Nombre del Script: Obtencion de gráficos a partir de datos de figuras---
## generados por Monolix GUI. 
##  
## Proposito del Script: crear y almacenar gráficos generados a partir de 
## los datos generados por la suite de Monolix, se debe colocar en la misma 
## carpeta en la que se encuentra el proyecto. Este script lee en la carpeta 
## ChartsData, que tiene como subdirectorios a cada gráfico generado por la 
## suite de Monolix.    
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion: 04-feb-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##=========================================================================#
##########################################################################-
# Introducción -----------------------------------------------------
##########################################################################-
# Carga de paquetes
require(tidyverse)
require(rlang)

##########################################################################-
# Selección de directorio principal
setwd(file.path('F:','Documentos','(Proyecto)_Estudio_PKPD','CEFEPIME',
                '02. Base','M2CPTM','Project_M2CPTM'))
##########################################################################-
# Carga de ajustes individuales
y_1_fits <- read_csv("ChartsData/IndividualFits/y_1_fits.txt")
y_1_observations <- read_csv("ChartsData/IndividualFits/y_1_observations.txt")

y_1_fits %>% 
ggplot(mapping = aes(x = time, y = pop)) +
  geom_line() +
  facet_wrap(~ID, ncol = 4)


y_1_observations %>% 
  ggplot(mapping = aes(x = time, y = y_1)) +
  geom_point() +
  facet_wrap(~ID, ncol = 4)






##########################################################################-
# Bondad de ajuste -------------------------------------------------
##########################################################################-
# Ajuste de una variable para guardar la subcarpeta que contiene datos de 
# gráfico de bondad de ajuste
gof_dir <- file.path('ChartsData','ObservationsVsPredictions')

# Selección de tema
theme_set(theme_classic() +
          theme(panel.border = element_rect(fill = NA, colour = 'black')))

# Apertura de archivo de datos
y_1_obsVsPred <- # Observaciones vs predicciones
  read_csv(file.path(gof_dir, 'y_1_obsVsPred.txt'))

y_1_visualGuides <- # Ayudas visuales
  read_csv(file.path(gof_dir, 'y_1_visualGuides.txt'))

y_1_obsVsSimulatedPred <- # Observaciones vs predicciones simuladas
  read_csv(file.path(gof_dir, 'y_1_obsVsSimulatedPred.txt'))

diagnostic_PRED <- function(data, x, y) {
  # Vector X
  u = data[[x]]
  # Vector Y
  v = data[[y]]
  # Número de datos
  n = length(u)
  # Matriz de diseño
  A = matrix(data = c(rep(1, times = n), u),
             nrow = n,
             ncol = 2)
  # Modelo LM
  W = lm.fit(y = v, x = A)
  # Varianza de error residual
  MSE = sum(W$residuals ^ 2) / (n - W$rank)
  # Suma cuadrados de X
  Sxx = sum((u - mean(u)) ^ 2)
  # Desviación estándar de predicciones
  SEy = sqrt(MSE * ((1 / n) + ((u - mean(u)) ^ 2 / Sxx)))
  # Error de predicción ponderado
  WPE = W$residuals / SEy
  # Error al cuadrado de predicción ponderado
  WSPE = WPE ^ 2
  # Sesgo
  S = sum(WPE / n)
  # Imprecisión
  I = sum(WSPE) / (n - S ^ 2)
  # Correlación
  R2 = cor(x = u, y = v)^2
  
  # SD pendiente
  Var_theta_0 = MSE * ((1 / n) + (mean(u) ^ 2 / Sxx))
  # SD intercepto
  Var_theta_1 = (MSE/Sxx)
  # Error
  Err = qt(p = 1-0.05/2, df = n-W$rank)
  
  # IC pendiente
  LI_pend = W$coefficients[2] - (Err * Var_theta_1)
  LS_pend = W$coefficients[2] + (Err * Var_theta_1)
  
  # IC intercepto
  LI_inter = W$coefficients[1] - (Err * Var_theta_0)
  LS_inter = W$coefficients[1] + (Err * Var_theta_0)
  
  # Efectos 
  Inter = paste0(
    "Inter = ",
    round(W$coefficients[1], 4),
    ", IC95%[",
    round(LI_inter, 4),
    "-",
    round(LS_inter, 4),
    "]"
  )
  Pend = paste0(
    "Pendiente = ",
    round(W$coefficients[2], 4),
    ", IC95%[",
    round(LI_pend, 4),
    "-",
    round(LS_pend, 4),
    "]"
  )
  
  return(list(Intercepto = Inter,
              Pendiente = Pend,
              Sesgo = round(S,4),
              Imprecision = round(I,4),
              R2 = R2))
}

##########################################################################-
# Bondad de ajuste OBS vs PRED
G_PRED_OBS_PRED <-
  y_1_obsVsPred %>% 
  ggplot(mapping = aes(x = popPred, y = y_1, group = ID)) +
  geom_point(shape = 1) + 
  xlab('PRED') + ylab('OBS') +
  geom_abline(slope = 1, intercept = 0, lty = 'dotted') + 
  # Modelo LOESS
  geom_line(y_1_visualGuides, 
            mapping =  aes(x = popPred_spline_abscissa, y = popPred_spline), 
            inherit.aes = F, colour = 'blue4') +
  geom_ribbon(y_1_visualGuides, 
              mapping =  aes(x = popPred_ci_abscissa, 
                             ymin = popPred_piLower,
                             ymax = popPred_piUpper), 
              inherit.aes = F, fill = 'blue1', alpha = 0.1) +
  coord_cartesian(xlim = c(0,90), ylim = c(0,90))

diagnostic_PRED(y_1_obsVsPred, "popPred", "y_1")

##########################################################################-
# Bondad de ajuste OBS vs IPRED
G_PRED_OBS_IPRED <-
  y_1_obsVsPred %>% 
  ggplot(mapping = aes(x = indivPredMean, y = y_1, group = ID)) +
  geom_point(shape = 1) + 
  xlab('IPRED') + ylab('OBS') +
  geom_abline(slope = 1, intercept = 0, lty = 'dotted') + 
  # Modelo LOESS
  geom_line(y_1_visualGuides, 
            mapping =  aes(x = indivPred_spline_abscissa, y = indivPred_spline), 
            inherit.aes = F, colour = 'red4') +
  geom_ribbon(y_1_visualGuides, 
              mapping =  aes(x = indivPred_ci_abscissa, 
                             ymin = indivPred_piLower,
                             ymax = indivPred_piUpper), 
              inherit.aes = F, fill = 'red1', alpha = 0.1) +
  coord_cartesian(xlim = c(0,90), ylim = c(0,90))
  
diagnostic_PRED(y_1_obsVsPred, "indivPredMean", "y_1")
##########################################################################-
# Bondad de ajuste OBS vs IPRED
G_PRED_OBS_PPRED <- 
  y_1_obsVsSimulatedPred    %>%  
    ggplot(mapping = aes(x = indivPredSimulated, y = y_1)) +
    geom_point(shape = 1) + 
    xlab('PPRED') + ylab('OBS') +
    stat_smooth(method = 'loess') +
    geom_abline(slope = 1, intercept = 0, lty = 'dotted') +
    coord_cartesian(xlim = c(0, 90), ylim = c(0, 90))
    
##########################################################################-
# Transformación de logarítmos
abreaks <- c(1, seq(2,10,2), seq(20,100,20))

G_PRED_OBS_IPREDLOG <-
  G_PRED_OBS_IPRED + 
  scale_y_continuous(trans = 'pseudo_log', breaks = abreaks) +
  scale_x_continuous(trans = 'pseudo_log', breaks = abreaks)

G_PRED_OBS_PREDLOG <-
  G_PRED_OBS_PRED + 
  scale_y_continuous(trans = 'pseudo_log', breaks = abreaks) +
  scale_x_continuous(trans = 'pseudo_log', breaks = abreaks)

# Almacenamiento en pdf de los gráficos
ggsave('FIGURAS/G_PRED_OBS_PRED.pdf', G_PRED_OBS_PRED, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')  

ggsave('FIGURAS/G_PRED_OBS_IPRED.pdf', G_PRED_OBS_IPRED, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')

ggsave('FIGURAS/G_PRED_OBS_PREDLOG.pdf', G_PRED_OBS_PREDLOG, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')  

ggsave('FIGURAS/G_PRED_OBS_IPREDLOG.pdf', G_PRED_OBS_IPREDLOG, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')

ggsave('FIGURAS/G_PRED_OBS_PPRED.pdf', G_PRED_OBS_PPRED, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')  

##########################################################################-
# Residuales --------------------------------------------------------------
##########################################################################-
# Ajuste de una variable para guardar la subcarpeta que contiene datos de 
# gráfico de residuales
res_dir <- file.path('ChartsData','ScatterPlotOfTheResiduals')
# Todos los Residuales
y_1_residuals <-
  read_csv(file.path(res_dir, 'y_1_residuals.txt'))
# Residuales simulados
y_1_simulatedResiduals <-
  read_csv(file.path(res_dir, 'y_1_simulatedResiduals.txt'))
# Percentiles de residuales vs prediccion
y_1_prediction_percentiles_iwRes <-
  read_csv(file.path(res_dir, 'y_1_prediction_percentiles_iwRes.txt'))
y_1_prediction_percentiles_pwRes <-
  read_csv(file.path(res_dir, 'y_1_prediction_percentiles_pwRes.txt'))
y_1_prediction_percentiles_npde <-
  read_csv(file.path(res_dir, 'y_1_prediction_percentiles_npde.txt'))
# Percentiles de residuales vs tiempo
y_1_time_percentiles_iwRes <-
  read_csv(file.path(res_dir, 'y_1_time_percentiles_iwRes.txt'))
y_1_time_percentiles_pwRes <-
  read_csv(file.path(res_dir, 'y_1_time_percentiles_pwRes.txt'))
y_1_time_percentiles_npde <-
  read_csv(file.path(res_dir, 'y_1_time_percentiles_npde.txt'))
# Especificicación de contenedores (bins) en los datos
y_1_individualBins <-
  read_csv(file.path(res_dir, 'y_1_individualBins.txt'))
y_1_populationBins <-
  read_csv(file.path(res_dir, 'y_1_populationBins.txt'))
y_1_timeBins <-
  read_csv(file.path(res_dir, 'y_1_timeBins.txt'))
# Especificación de línea de tendencia
y_1_spline <-
  read_csv(file.path(res_dir, 'y_1_spline.txt'))

##########################################################################-
  # Gráficos de residuales vs Tiempo --------------------------------------
##########################################################################-

G_RES_T_PWRES <-
  y_1_residuals %>% 
    ggplot(mapping = aes(x = time, y = pwRes)) +
    geom_hline(yintercept = 0) +
    geom_point(col = '#4682B4') +
    geom_line(y_1_spline, 
              mapping = aes(x = time_pwRes, y = time_pwRes_spline), 
              col = '#EBA213', lty = 'solid', size = 1) +
    geom_line(y_1_time_percentiles_pwRes, 
              mapping = aes(x = time, y = empirical_median), 
              col = '#EBA213', lty = 'dashed', size = 1.2) + 
    geom_line(y_1_time_percentiles_pwRes, 
            mapping = aes(x = time, y = empirical_lower), 
            col = '#EBA213', lty = 'dotted', size = 1.2) +
    geom_line(y_1_time_percentiles_pwRes, 
              mapping = aes(x = time, y = empirical_upper), 
              col = '#EBA213', lty = 'dotted', size = 1.2) +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    xlab('TSFD') + ylab('WRES')
  
G_RES_T_IWRES <-
  y_1_residuals %>% 
    ggplot(mapping = aes(x = time, y = iwRes_mean)) +
    geom_hline(yintercept = 0) +
    geom_line(y_1_spline, 
              mapping = aes(x = time_pwRes, y = time_iwRes_spline), 
              col = '#EBA213', lty = 'solid', size = 1) +
    
    geom_line(y_1_time_percentiles_iwRes, 
              mapping = aes(x = time, y = empirical_median), 
              col = '#EBA213', lty = 'dashed', size = 1.2) + 
    geom_line(y_1_time_percentiles_iwRes, 
              mapping = aes(x = time, y = empirical_lower), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(y_1_time_percentiles_iwRes, 
              mapping = aes(x = time, y = empirical_upper), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_point(col = '#4682B4') +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    xlab('TSFD') + ylab('IWRES')
  
G_RES_T_NPDE <-
  y_1_residuals %>% 
    ggplot(mapping = aes(x = time, y = npde)) +
    geom_hline(yintercept = 0) +
    geom_line(y_1_spline, 
              mapping = aes(x = time_pwRes, y = time_npde_spline), 
              col = '#EBA213', lty = 'solid', size = 1) +
    geom_line(y_1_time_percentiles_npde, 
              mapping = aes(x = time, y = empirical_median), 
              col = '#EBA213', lty = 'dashed', size = 1.2) + 
    geom_line(y_1_time_percentiles_npde, 
              mapping = aes(x = time, y = empirical_lower), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(y_1_time_percentiles_npde, 
              mapping = aes(x = time, y = empirical_upper), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_point(col = '#4682B4') +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    xlab('TSFD') + ylab('NPDE')
  

ggsave(filename = './FIGURAS/G_RES_T_PWRES.pdf', plot = G_RES_T_PWRES, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')

ggsave(filename = './FIGURAS/G_RES_T_IWRES.pdf', plot = G_RES_T_IWRES, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')  

ggsave(filename = './FIGURAS/G_RES_T_NPDE.pdf', plot = G_RES_T_NPDE, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in') 

##########################################################################-
# Gráfico de residuales vs TAD --------------------------------------------
##########################################################################-
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Apertura de archivo de datos con los datos originales, y que tiene 
##  delimitador ";", este fue usado para realizar el modelamiento en sí mismo
##  2 Crear un archivo *last_dose* que contiene la última hora de adminis-
##  tración del antibiótico en el archivo *data* para cada individuo. 
##    2a Agrupar por ID
##    2b Filtrar los datos con EVID igual a 1 (sólo los eventos de 
##    administración)
##    2c Seleccionar los últimos en cada grupo (los datos ya están ordenados 
##    de menor a mayor)
##    2d Seleccionar sólo las columnas ID y TIME
##  3 Modificar el archivo de residuales *y_1_residuals* al unir los datos 
##  de últimnas observaciones calculadas.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data <- read_delim("../1_DATA/Monolix_data.csv", delim = ';')

last_dose <- data %>% 
  group_by(ID) %>% 
  filter(EVID == 1) %>% 
  slice(n()) %>% 
  select(ID, TIME)

y_1_residuals <-
  y_1_residuals %>% 
  left_join(., last_dose, by = 'ID') %>% 
  mutate(TAD = time - TIME)

RES_TAD <- function(x, y, xlab, ylab) {
  ##########################################################################-
  # Volver las variables x y y en expresiones para evaluación tardía
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ##  Calculo de percentiles empiricos
  ##  1 Tomar la variable de datos a explorar.
  ##  2 Crear una variable que tome 6 grupos a partir de los datos ordenados 
  ##  con la función dplyr::ntile().
  ##  3 Agrupar por la variable recién creado.
  ##  4 Resumir por la media de TAD, y los percentiles P5%, P50%, y P95% de 
  ##  los datos de residuales.
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  A <- y_1_residuals %>%
    mutate(gr = ntile(TAD, 6)) %>% 
    group_by(gr) %>% 
    summarise(TIME = mean(TAD), 
              ME = quantile(x = !!y, probs = 0.50), 
              LI = quantile(x = !!y, probs = 0.05), 
              LS = quantile(x = !!y, probs = 0.95))
  ##########################################################################-
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
              col = '#EBA213', lty = 'dashed', size = 1.2) %>% 
    return(.)
  }
 
G_RES_TAD_PWRES <- 
  RES_TAD(x = TAD, y = pwRes, xlab = 'TAD', ylab = 'WRES')

G_RES_TAD_IWRES <- 
  RES_TAD(x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')

G_RES_TAD_NPDE <- 
  RES_TAD(x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')

##########################################################################-
# Almacenamiento de los datos en formato PDF
ggsave(filename = './FIGURAS/G_RES_TAD_PWRES.pdf', plot = G_RES_TAD_PWRES, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')

ggsave(filename = './FIGURAS/G_RES_TAD_IWRES.pdf', plot = G_RES_TAD_IWRES, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')

ggsave(filename = './FIGURAS/G_RES_TAD_NPDE.pdf', plot = G_RES_TAD_NPDE, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')

##########################################################################-
# Gráficos de residuales vs Concentraciones -------------------------------
##########################################################################-
# Residuales vs concentraciones

#' Creación de gráfico de residuales
#'
#' @param x Variable X (en y_1_residuals) para puntos
#' @param y Variable Y (en y_1_residuals) para puntos
#' @param xspline Variable X (en y_1_spline) para spline
#' @param yspline Variable Y (en y_1_spline) para spline
#' @param perc_data Tabla con datos de percentiles empíricos
#' @param xlab Etiqueta de eje X personalizada
#' @param ylab Etiqueta de eje Y personalizada
#'
#' @return Gráfico de residuales vs predicciones
#' @export
#'
#' @examples
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
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    xlab(xlab) + ylab(ylab) %>% return(.)
}

G_RES_C_PWRES <- 
  RES_PRE(x = prediction_pwRes, y = pwRes, 
          xspline = prediction_pwRes, yspline = prediction_pwRes_spline,
          perc_data = y_1_prediction_percentiles_pwRes,
          xlab = 'PRED', ylab = 'WRES')

G_RES_C_IWRES <- 
  RES_PRE(x = prediction_iwRes_mean, y = iwRes_mean, 
          xspline = prediction_iwRes, yspline = prediction_iwRes_spline,
          perc_data = y_1_prediction_percentiles_iwRes,
          xlab = 'IPRED', ylab = 'IWRES')

G_RES_C_NPDE <- 
  RES_PRE(x = prediction_npde, y = npde, 
          xspline = prediction_npde, yspline = prediction_npde_spline,
          perc_data = y_1_prediction_percentiles_npde,
          xlab = 'PRED', ylab = 'NPDE')
  

# Almacenamiento en pdf de los gráficos
ggsave('./FIGURAS/G_RES_C_PWRES.pdf', G_RES_C_PWRES, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')

ggsave('./FIGURAS/G_RES_C_IWRES.pdf', G_RES_C_IWRES, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')  

ggsave('./FIGURAS/G_RES_C_NPDE.pdf', G_RES_C_NPDE, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')  


##########################################################################-
# Cálculo manual de encogimiento eta y epsilon ----------------------------
##########################################################################-
eta <- 
  read_csv("ChartsData/CorrelationBetweenRandomEffects/eta.txt")


eta_shrink_V1 <- 1 - (var(eta$eta_V1_mean) / (0.297 ^ 2))
eta_shrink_V2 <- 1 - (var(eta$eta_V2_mean) / (0.822 ^ 2))
eta_shrink_Cl <- 1 - (var(eta$eta_Cl_mean) / (0.262 ^ 2))
eta_shrink_Q <- 1 - (var(eta$eta_Q_mean) / (0.836 ^ 2))

epsilon_shrink <- 1-sd(y_1_residuals$iwRes_mean)




##########################################################################-
# Gráfico de densidad distribución parámetros -----------------------------
##########################################################################-
# 
data1 <-
    read_csv("ChartsData/DistributionOfTheIndividualParameters/cdf.txt")


  
  
  