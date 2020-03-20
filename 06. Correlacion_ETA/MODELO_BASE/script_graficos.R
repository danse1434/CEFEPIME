##=========================================================================#
## Nombre del Script: Obtencion de gr?ficos a partir de datos de figuras---
## generados por Monolix GUI. 
##  
## Proposito del Script: crear y almacenar gr?ficos generados a partir de 
## los datos generados por la suite de Monolix, se debe colocar en la misma 
## carpeta en la que se encuentra el proyecto. Este script lee en la carpeta 
## ChartsData, que tiene como subdirectorios a cada gr?fico generado por la 
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
# Introducci?n -----------------------------------------------------
##########################################################################-
# Carga de paquetes
require(tidyverse)
require(rlang)

##########################################################################-
# Selecci?n de directorio principal
setwd(file.path('C:','Users','Daniel','OneDrive','Documents','(Proyecto)_Estudio_PKPD','CEFEPIME',
                '04. Residual', 'Modelo_8', 'RES_M1'))

##########################################################################-
# Bondad de ajuste -------------------------------------------------
##########################################################################-
# Ajuste de una variable para guardar la subcarpeta que contiene datos de 
# gr?fico de bondad de ajuste
gof_dir <- file.path('ChartsData','ObservationsVsPredictions')

# Selecci?n de tema
theme_set(theme_classic() +
          theme(panel.border = element_rect(fill = NA, colour = 'black')))

# Apertura de archivo de datos
y_1_obsVsPred <- # Observaciones vs predicciones
  read_csv(file.path(gof_dir, 'y_1_obsVsPred.txt'))

y_1_visualGuides <- # Ayudas visuales
  read_csv(file.path(gof_dir, 'y_1_visualGuides.txt'))

y_1_obsVsSimulatedPred <- # Observaciones vs predicciones simuladas
  read_csv(file.path(gof_dir, 'y_1_obsVsSimulatedPred.txt'))

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

##########################################################################-
GOF_PRED <- function(x, y, xspline, yspline, xconfint, yconfint_lo, 
                     yconfint_up, colourp, xlab, ylab) {
  x = rlang::ensym(x)
  y = rlang::ensym(y)
  xspline = rlang::ensym(xspline)
  yspline = rlang::ensym(yspline)
  xconfint = rlang::ensym(xconfint)
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

G_PRED_OBS_PRED <- 
  GOF_PRED(x = popPred, y = y_1, 
         xspline = popPred_spline_abscissa, 
         yspline = popPred_spline, 
         xconfint = popPred_ci_abscissa, 
         yconfint_lo = popPred_piLower, 
         yconfint_up = popPred_piUpper, 
         colourp = 'blue1', 
         xlab = 'PRED', ylab = 'OBS')

G_PRED_OBS_IPRED <- 
  GOF_PRED(x = indivPredMean, y = y_1, 
           xspline = indivPred_spline_abscissa, 
           yspline = indivPred_spline, 
           xconfint = indivPred_ci_abscissa, 
           yconfint_lo = indivPred_piLower, 
           yconfint_up = indivPred_piUpper, 
           colourp = 'red1', 
           xlab = 'IPRED', ylab = 'OBS')


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
# Transformaci?n de logar?tmos
abreaks <- c(1, seq(2,10,2), seq(20,100,20))

G_PRED_OBS_IPREDLOG <-
  G_PRED_OBS_IPRED + 
  scale_y_continuous(trans = 'pseudo_log', breaks = abreaks) +
  scale_x_continuous(trans = 'pseudo_log', breaks = abreaks)

G_PRED_OBS_PREDLOG <-
  G_PRED_OBS_PRED + 
  scale_y_continuous(trans = 'pseudo_log', breaks = abreaks) +
  scale_x_continuous(trans = 'pseudo_log', breaks = abreaks)

# Almacenamiento en pdf de los gr?ficos
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
# gr?fico de residuales
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
# Especificicaci?n de contenedores (bins) en los datos
y_1_individualBins <-
  read_csv(file.path(res_dir, 'y_1_individualBins.txt'))
y_1_populationBins <-
  read_csv(file.path(res_dir, 'y_1_populationBins.txt'))
y_1_timeBins <-
  read_csv(file.path(res_dir, 'y_1_timeBins.txt'))
# Especificaci?n de l?nea de tendencia
y_1_spline <-
  read_csv(file.path(res_dir, 'y_1_spline.txt'))

##########################################################################-
  # Gr?ficos de residuales vs Tiempo --------------------------------------
##########################################################################-
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

G_RES_T_PWRES <- RES_TSFD(
  x = time,
  y = pwRes,
  xspline = time_pwRes,
  yspline = time_pwRes_spline,
  perc_data = y_1_time_percentiles_pwRes,
  xlab = 'TSFD',
  ylab = 'PWRES'
)

G_RES_T_IWRES <- RES_TSFD(
  x = time,
  y = iwRes_mean,
  xspline = time_iwRes,
  yspline = time_iwRes_spline,
  perc_data = y_1_time_percentiles_iwRes,
  xlab = 'TSFD',
  ylab = 'IWRES'
)

G_RES_T_NPDE <- RES_TSFD(
  x = time,
  y = npde,
  xspline = time_npde,
  yspline = time_npde_spline,
  perc_data = y_1_time_percentiles_npde,
  xlab = 'TSFD',
  ylab = 'NPDE'
)

##########################################################################-
# Almacenar los archivos

ggsave(filename = './FIGURAS/G_RES_T_PWRES.pdf', plot = G_RES_T_PWRES, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')

ggsave(filename = './FIGURAS/G_RES_T_IWRES.pdf', plot = G_RES_T_IWRES, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')  

ggsave(filename = './FIGURAS/G_RES_T_NPDE.pdf', plot = G_RES_T_NPDE, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in') 

##########################################################################-
# Gr?fico de residuales vs TAD --------------------------------------------
##########################################################################-
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Apertura de archivo de datos con los datos originales, y que tiene 
##  delimitador ";", este fue usado para realizar el modelamiento en s? mismo
##  2 Crear un archivo *last_dose* que contiene la ?ltima hora de adminis-
##  traci?n del antibi?tico en el archivo *data* para cada individuo. 
##    2a Agrupar por ID
##    2b Filtrar los datos con EVID igual a 1 (s?lo los eventos de 
##    administraci?n)
##    2c Seleccionar los ?ltimos en cada grupo (los datos ya est?n ordenados 
##    de menor a mayor)
##    2d Seleccionar s?lo las columnas ID y TIME
##  3 Modificar el archivo de residuales *y_1_residuals* al unir los datos 
##  de ?ltimnas observaciones calculadas.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data <- read_delim("Monolix_data.csv", delim = ';')

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
  # Volver las variables x y y en expresiones para evaluaci?n tard?a
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ##  Calculo de percentiles empiricos
  ##  1 Tomar la variable de datos a explorar.
  ##  2 Crear una variable que tome 6 grupos a partir de los datos ordenados 
  ##  con la funci?n dplyr::ntile().
  ##  3 Agrupar por la variable reci?n creado.
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
  # Crear el gr?fico de residuales con las especificaciones mostradas
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
  RES_TAD(x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')

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
# Gr?ficos de residuales vs Concentraciones -------------------------------
##########################################################################-
# Residuales vs concentraciones

#' Creaci?n de gr?fico de residuales
#'
#' @param x Variable X (en y_1_residuals) para puntos
#' @param y Variable Y (en y_1_residuals) para puntos
#' @param xspline Variable X (en y_1_spline) para spline
#' @param yspline Variable Y (en y_1_spline) para spline
#' @param perc_data Tabla con datos de percentiles emp?ricos
#' @param xlab Etiqueta de eje X personalizada
#' @param ylab Etiqueta de eje Y personalizada
#'
#' @return Gr?fico de residuales vs predicciones
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
  

# Almacenamiento en pdf de los gr?ficos
ggsave('./FIGURAS/G_RES_C_PWRES.pdf', G_RES_C_PWRES, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')

ggsave('./FIGURAS/G_RES_C_IWRES.pdf', G_RES_C_IWRES, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')  

ggsave('./FIGURAS/G_RES_C_NPDE.pdf', G_RES_C_NPDE, 
       device = 'pdf', width = 5/1.5, height = 4/1.5, units = 'in')  

##########################################################################-
# Gr?ficos individuales-----------------------------------------------------
##########################################################################-
# Carga de ajustes individuales
y_1_fits <- read_csv("ChartsData/IndividualFits/y_1_fits.txt")
y_1_observations <- read_csv("ChartsData/IndividualFits/y_1_observations.txt") %>% 
  left_join(., last_dose, by = 'ID') %>% 
  mutate(TAD = time - TIME)

G_CP_TAD <-
  y_1_fits %>%
  left_join(., last_dose, by = 'ID') %>%
  mutate(TAD = time - TIME) %>%
  
  ggplot(mapping = aes(x = TAD)) +
  geom_line(aes(y = pop), lty = 'dashed') +
  geom_line(aes(y = indivPredMean), lty = 'solid', col = 'green3') +
  facet_wrap( ~ ID, ncol = 4, labeller = labeller(.cols = label_both)) +
  geom_point(data = y_1_observations,
             mapping = aes(x = TAD, y = y_1),
             col = 'green4') +
  xlim(0, NA) +
  xlab('TAD') + ylab(expression(C[P]))
  
ggsave(filename = 'FIGURAS/G_CP_TAD.pdf', plot = G_CP_TAD,
       device = 'pdf', width = 5, height = 6)

##########################################################################-
# C?lculo manual de encogimiento eta y epsilon ----------------------------
##########################################################################-
eta <- 
  read_csv("ChartsData/CorrelationBetweenRandomEffects/eta.txt")

pop_eta <-
  read_csv("populationParameters.txt")

pop_eta1 <- pop_eta %>% 
  rownames_to_column %>% 
  gather(var, value, -rowname) %>% 
  spread(rowname, value) 

colnames(pop_eta1) <- pop_eta1[1,]

pop_eta1 <- pop_eta1[-1, ] %>%
  filter(parameter == 'value') %>%
  select_at(vars(matches("omega\\_"))) %>% 
  mutate_all(~as.double(.))


var_eta <- eta %>% 
  summarise_at(., vars(matches("\\_mean")), .funs = list(var = ~ var(.))) %>%
  rename_all(.funs = list(~ str_replace(., "\\_mean\\_var", ''))) %>%
  rename_all(.funs = list(~ str_replace(., "eta\\_", ''))) %>%
  rename_all(.funs = list(~ paste0(., '_pop'))) 
  
# Resultados de Eta
eta_res <- 1 - (var_eta / (pop_eta1^ 2))

# Resultados de IWRES
epsilon_shrink <- 1-sd(y_1_residuals$iwRes_mean)




##########################################################################-
# Gr?fico de densidad distribuci?n par?metros -----------------------------
##########################################################################-
# 
data1 <-
    read_csv("ChartsData/DistributionOfTheIndividualParameters/cdf.txt")


  
  
  