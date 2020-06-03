##------------------------------------------------------------------------------#
## Nombre del Script: Análisis de Malla de Lambda en TBS ------------------------
##  
## Propósito del Script: Análisis de grid de lambda con TBS introducido de manera 
## artificial en variable dependiente e independiente. 
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 03-06-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
# Carga de paquetes
require(ggrepel)
require(animation)
require(patchwork)
require(rlang)
require(tidyverse)
require(gganimate)

#-------------------------------------------------------------------------------#
# Introducción ------------------------------------------------------------------
#-------------------------------------------------------------------------------#
# Vector de valores de lambda introducidos
lamb_vec <- c(seq(-3, -0.1, 0.1), seq(+0.1, +3, 0.1)) 

#-------------------------------------------------------------------------------#
# Análisis de Log-verosimilitud -------------------------------------------------
#-------------------------------------------------------------------------------#
# Pre-alocación de **log_ls**
log_ls <- list()

# Búsqueda de registros de verosimilitud
for (i in 1:60) {
  dis <- format(round(lamb_vec[[i]], 2), nsmall = 2)
  dir <- file.path('tbs_folder_1',paste0('lambda_',dis),
                   'model_TBS','LogLikelihood','logLikelihood.txt')
  log_ls[[i]] <- read_csv(dir)
}

# Transformación de LL como data.frame
#................................................................................
#   1 Conversión en data.frame
#   2 Adicionar una columna con valores de lambda
#   3 Expandir tabla con columnas con los criterios: BIC, AIC, LL, SE
#................................................................................

log_ls1 <- log_ls %>% 
  map_dfr(~.x, .id = 'ID') %>%
  add_column(lambda = rep(lamb_vec, each = 5), .before = 'criteria') %>% 
  pivot_wider(names_from = criteria, values_from = importanceSampling)

log_ls1 %>% 
  mutate(dif2LL = (`-2LL` - lag(`-2LL`))/(lambda-lag(lambda))) %>% 
  ggplot(aes(x = lambda, y = dif2LL)) +
  geom_line()

# Ajuste de tema de gráficos
theme_set(theme_bw())

# Gráfico de seguimiento de OFV y error estándar del modelo
g1 <- ggplot(log_ls1, aes(x = lambda, y = `-2LL`)) +
  geom_line(col='blue1') + geom_point(col='blue4') + 
  xlab(expression(lambda)) 

g2 <- ggplot(log_ls1, aes(x = lambda, y = stdError)) +
  geom_line(col='red1') + geom_point(col='red4') + 
  xlab(expression(lambda)) + ylab('Error Estándar')

ggsave('figures/1_ML_grid.pdf', g1 + g2, 'pdf', width = 6, height = 3)

#-------------------------------------------------------------------------------#
# Pruebas de normalidad ----------------------------------------------------------
#-------------------------------------------------------------------------------#
# Pre-alocación de **nor_ls**
nor_ls <- list()

# Búsqueda de registros de pruebas de normalidad
for (i in 1:60) {
  dis <- format(round(lamb_vec[[i]], 2), nsmall = 2)
  dir <- file.path('tbs_folder_1',paste0('lambda_',dis),
                   'model_TBS','Tests','normalityResiduals.txt')
  nor_ls[[i]] <- read_csv(dir)
}

# 
#................................................................................
# 1 Convertir a la columna p en carácter, debido a que en algunos df hay valores 
#   p indicados por intervalo.
# 2 Convertir la lista en un data.frame
# 3 Adicionar columna con los valores de lambda correspondientes.
# 4 Convertir el valor p en variable numérica
#................................................................................

nor_ls1 <- nor_ls %>%
  map(~ mutate(.x, `p-value` = as.character(`p-value`))) %>% 
  map_dfr(~.x, .id = 'ID') %>%
  add_column(lambda = rep(lamb_vec, each = 3), .before = 'residuals') %>% 
  mutate(`p-value`=as.numeric(str_replace_all(`p-value`, "\\<", '')))


#-------------------------------------------------------------------------------#
#' Crear gráfico de resultados de pruebas de normalidad 
#'
#' @param data  datos
#' @param resid tipo de residual iwRes_mean, pwRes, y npde
#' @param yvar  variable *y* que puede ser estadístico t o valor p
#' @param col   color de tema gráfico
#' @param ylim  (vector) límites del eje y
#'
#' @return gráfico ggplot2 con resultados de residuales por lambda 
#' @export
#'
#' @examples
#' graf.res(nor_ls1, PWRES_y_1, statistics, ylim = c(0,1))
#' 
graf.res <- function(data, resid, yvar, col='blue4', ylim, 
                     sub_label) {
  
  if (missing(sub_label)) {sub_label = 'TBS - Normalidad: '}
  
  resid_quo <- rlang::ensym(resid)
  yvar_quo  <- rlang::ensym(yvar)
  
  switch (expr_text(resid_quo),
    "IWRES_y_1" = {sublab = 'IWRES'},
    "PWRES_y_1" = {sublab = 'PWRES'},
    "NPDE_y_1"  = {sublab = 'NPDE'},
  )
  
  switch (expr_text(yvar_quo),
    "`p-value`"  = {ylab = 'Valor p'},
    "statistics" = {ylab = 'Estadístico t'}
  )
  
  data %>% 
    filter(residuals == expr_text(resid_quo)) %>% 
  ggplot(aes(lambda, !!yvar_quo)) +
  geom_line(col=col) + geom_point(col=col) +
  labs(subtitle = paste0(sub_label, sublab)) +
  xlab(expression(lambda)) + ylab(ylab) +
  coord_cartesian(ylim = ylim)
}

# Lista de gráficos de RES vs lambda
p_ls <- list()

p_ls[[1]] <- graf.res(nor_ls1, PWRES_y_1, statistics, ylim = c(0,1))
p_ls[[2]] <- graf.res(nor_ls1, IWRES_y_1, statistics, 'red2', ylim = c(0,1))
p_ls[[3]] <- graf.res(nor_ls1, NPDE_y_1, statistics, 'green3', ylim = c(0,1)) 

p_ls[[4]] <- graf.res(nor_ls1, PWRES_y_1, "p-value", ylim = c(0,1)) +
  geom_hline(yintercept = 0.05, lty = 'dotted')
p_ls[[5]] <- graf.res(nor_ls1, IWRES_y_1, "p-value", 'red2', ylim = c(0,1)) +
  geom_hline(yintercept = 0.05, lty = 'dotted')
p_ls[[6]] <- graf.res(nor_ls1, NPDE_y_1, "p-value", 'green3',ylim = c(0,1)) +
  geom_hline(yintercept = 0.05, lty = 'dotted')

# Creación de gráfico conjugado
p_ls_conj <- 
(p_ls[[1]] + p_ls[[4]]) /
(p_ls[[2]] + p_ls[[5]]) /
(p_ls[[3]] + p_ls[[6]])

ggsave('figures/2_seguimiento_normalidad.pdf', p_ls_conj, 'pdf', 
       width = 6, height = 6, units = 'in')

#-------------------------------------------------------------------------------#
# Pruebas de simetría -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Prea-locación de lista con resultados de simetría
sim_ls <- list()

# Búsqueda de registros de pruebas de simetría
for (i in 1:60) {
  dis <- format(round(lamb_vec[[i]], 2), nsmall = 2)
  dir <- file.path('tbs_folder_1',paste0('lambda_',dis),
                   'model_TBS','Tests','symmetryResiduals.txt')
  sim_ls[[i]] <- read_csv(dir)
}

# 
#................................................................................
# 1 Convertir a la columna p en carácter, debido a que en algunos df hay valores 
#   p indicados por intervalo.
# 2 Convertir la lista en un data.frame
# 3 Adicionar columna con los valores de lambda correspondientes.
# 4 Convertir el valor p en variable numérica
#................................................................................

sim_ls1 <- sim_ls %>%
  map(~ mutate(.x, `p-value` = as.character(`p-value`))) %>% 
  map_dfr(~.x, .id = 'ID') %>%
  add_column(lambda = rep(lamb_vec, each = 3), .before = 'residuals') %>% 
  mutate(`p-value`=as.numeric(str_replace_all(`p-value`, "\\<", '')))

sim_ls1

# Lista de gráficos de RES vs lambda
s_ls <- list()

s_ls[[1]] <- graf.res(sim_ls1, PWRES_y_1, statistics, ylim = c(-5,10), 
                      sub_label = 'TBS - Simetría: ')
s_ls[[2]] <- graf.res(sim_ls1, IWRES_y_1, statistics, 'red2', ylim = c(-5,10), 
                      sub_label = 'TBS - Simetría: ')
s_ls[[3]] <- graf.res(sim_ls1, NPDE_y_1, statistics, 'green3', ylim = c(-5,10), 
                      sub_label = 'TBS - Simetría: ') 
s_ls[[4]] <- graf.res(sim_ls1, PWRES_y_1, "p-value", ylim = c(0,1), 
                      sub_label = 'TBS - Simetría: ') +
  geom_hline(yintercept = 0.05, lty = 'dotted')
s_ls[[5]] <- graf.res(sim_ls1, IWRES_y_1, "p-value", 'red2', ylim = c(0,1), 
                      sub_label = 'TBS - Simetría: ') +
  geom_hline(yintercept = 0.05, lty = 'dotted')
s_ls[[6]] <- graf.res(sim_ls1, NPDE_y_1, "p-value", 'green3',ylim = c(0,1),
                      sub_label = 'TBS - Simetría: ') +
  geom_hline(yintercept = 0.05, lty = 'dotted')

# Creación de gráfico conjugado
s_ls_conj <- 
  (s_ls[[1]] + s_ls[[4]]) /
  (s_ls[[2]] + s_ls[[5]]) /
  (s_ls[[3]] + s_ls[[6]])

s_ls_conj

ggsave('figures/3_seguimiento_simetría.pdf', s_ls_conj, 'pdf', 
       width = 6, height = 6, units = 'in')

#-------------------------------------------------------------------------------#
# QQplot de Residuales -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Prealocación de lista
res_ls <- list()

# Búsqueda de registros de residuales
for (i in 1:60) {
  dis <- format(round(lamb_vec[[i]], 2), nsmall = 2)
  dir <- file.path('tbs_folder_1',paste0('lambda_',dis),
                   'model_TBS', 'ChartsData','ScatterPlotOfTheResiduals',
                   'y_1_residuals.txt')
  res_ls[[i]] <- read_csv(dir)
}

# Conversión de lista a data.frame
res_ls1 <- res_ls %>%
  map_dfr(~.x, .id = 'ID') %>%
  add_column(lambda = rep(lamb_vec, each = 90), .before = 'time') 

res_ls1

#-------------------------------------------------------------------------------#
#' Función de gráficos QQ
#'
#' @param data datos
#' @param param parámetro a graficar
#' @param colour color de temático de gráfico
#'
#' @return gráfico de tipo QQ en ggplot2 de dispersión (sin identificación de ID), 
#' con línea unitaria, línea de tendencia lineal y bandas de confianza.
#' @export 
#' @examples
#' normal_dist_plot(1, 'iwRes_mean','red3')
#' normal_dist_plot(2, 'npde','red3')
#' 
normal_dist_plot <- function(data, param, colour, pnorm, psimt) {
  
  param_quo <- rlang::ensym(param)
  if(is_missing(colour)){colour = 'blue4'}
  
  g1 <- data %>%
    ggplot(mapping = aes(sample = !!param_quo)) +
    stat_qq() +
    qqplotr::stat_qq_line(distribution = 'norm', col = colour) +
    qqplotr::stat_qq_band(distribution = 'norm', fill = colour, alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = -3:+3, minor_breaks = seq(-3,3,0.5)) + 
    scale_y_continuous(breaks = -3:+3, minor_breaks = seq(-3,3,0.5)) + 
    coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) + 
    xlab('Cuantiles teóricos') + ylab('Cuantiles de muestreo') + 
    theme(panel.grid.major = element_line(colour = 'gray90'),
          panel.grid.minor = element_line(colour = 'gray97')) +
    labs(title = glue::glue('QQplot - Efecto de TBS en {quo_name(param_quo)}'), 
         subtitle = glue::glue('Valor de lambda = {lamb_vec[[i]]}; \\
                               Valor-p (norm.): {pnorm}; \\
                               Valor-p (simt.): {psimt}'))
  
  return(g1)
}

# Gráfico animado de PWRES
saveGIF({
  for (i in 1:60) {
    res_ls1_frame = res_ls1 %>%
      filter(lambda == lamb_vec[[i]])
    
    pnorm = nor_ls1 %>% 
      filter(lambda==lamb_vec[[i]] & residuals == 'NPDE_y_1') %>% 
      magrittr::use_series(`p-value`) %>%  round(3)
    
    psimt = nor_ls1 %>% 
      filter(lambda==lamb_vec[[i]] & residuals == 'NPDE_y_1') %>% 
      magrittr::use_series(`p-value`) %>% round(3)
    
    p = normal_dist_plot(res_ls1_frame, pwRes, pnorm=pnorm, psimt=psimt)
    print(p)
  }
}, movie.name = 'qplot_PWRES_TB.gif', interval = 0.8, 
ani.width = 720*0.85, ani.height = 480*0.85)

#-------------------------------------------------------------------------------#
# Gráfico animado IWRES
saveGIF({
  for (i in 1:60) {
    res_ls1_frame = res_ls1 %>%
      filter(lambda == lamb_vec[[i]])
    
    pnorm = nor_ls1 %>% 
      filter(lambda==lamb_vec[[i]] & residuals == 'NPDE_y_1') %>% 
      magrittr::use_series(`p-value`) %>%  round(3)
    
    psimt = nor_ls1 %>% 
      filter(lambda==lamb_vec[[i]] & residuals == 'NPDE_y_1') %>% 
      magrittr::use_series(`p-value`) %>% round(3)
    
    p = normal_dist_plot(res_ls1_frame, iwRes_mean, pnorm=pnorm, psimt=psimt, 
                         colour = 'red3')
    print(p)
  }
}, movie.name = 'qplot_IWRES_TB.gif', interval = 0.8, 
ani.width = 720*0.85, ani.height = 480*0.85)

#-------------------------------------------------------------------------------#
# Gráfico animado NPDE
saveGIF({
  for (i in 1:60) {
    res_ls1_frame = res_ls1 %>%
      filter(lambda == lamb_vec[[i]])
    
    pnorm = nor_ls1 %>% 
      filter(lambda==lamb_vec[[i]] & residuals == 'NPDE_y_1') %>% 
      magrittr::use_series(`p-value`) %>%  round(3)
    
    psimt = nor_ls1 %>% 
      filter(lambda==lamb_vec[[i]] & residuals == 'NPDE_y_1') %>% 
      magrittr::use_series(`p-value`) %>% round(3)
    
    p = normal_dist_plot(res_ls1_frame, npde, pnorm=pnorm, psimt=psimt, 
                         colour = 'green3')
    print(p)
  }
}, movie.name = 'qplot_NPDE_TB.gif', interval = 0.8, 
ani.width = 720*0.85, ani.height = 480*0.85)

#-------------------------------------------------------------------------------#
# Distribución de residuales -----------------------------------------------------
#-------------------------------------------------------------------------------#
#' Gráficar histograma de residuales vs lambda con tiempo
#'
#' @param data archivo de datos de residuales
#' @param xvar tipo de residual a graficar
#' @param col color de tema de gráfico
#'
#' @return
#' @export
#'
#' @examples
#' 
hist.residual <- function(data, xvar, col) {
  xvar_quo <- rlang::enquo(xvar)
  
  data %>% 
    ggplot(aes(x = !!xvar_quo)) +
    geom_vline(xintercept = 0) +
    geom_histogram(aes(y = ..density..), fill = col, col='black', bins = 20) +
    geom_density(aes(y = ..density..), col=col, fill=alpha(col, 0.2)) +
    coord_cartesian(xlim = c(-4,+4), ylim = c(0,1)) +
    ylab('Densidad') + 
    xlab(toupper(quo_name(xvar_quo))) +
    labs(title = glue::glue('Histograma - Efecto de TBS en {quo_name(xvar_quo)}'), 
         subtitle = glue::glue('Valor de lambda = {lamb_vec[[i]]}'))
}


# Histograma + Densidad Animado de PWRES
saveGIF({
  for (i in 1:60) {
    res_ls1_frame = res_ls1 %>%
      filter(lambda == lamb_vec[[i]])
    
    p = res_ls1_frame %>% 
      hist.residual(pwRes, 'blue4') 
    print(p)
  }
}, movie.name = 'histograma_PWRES_TBS.gif', interval = 0.8, 
ani.width = 720*0.85, ani.height = 480*0.85)

# Histograma + Densidad Animado de IWRES
saveGIF({
  for (i in 1:60) {
    res_ls1_frame = res_ls1 %>%
      filter(lambda == lamb_vec[[i]])
    
    p = res_ls1_frame %>% 
      hist.residual(iwRes_mean, 'red3') 
    print(p)
  }
}, movie.name = 'histograma_IWRES_TBS.gif', interval = 0.8, 
ani.width = 720*0.85, ani.height = 480*0.85)

# Histograma + Densidad Animado de NPDE
saveGIF({
  for (i in 1:60) {
    res_ls1_frame = res_ls1 %>%
      filter(lambda == lamb_vec[[i]])
    
    p = res_ls1_frame %>% 
      hist.residual(npde, 'green2') 
    print(p)
  }
}, movie.name = 'histograma_NPDE_TBS.gif', interval = 0.8, 
ani.width = 720*0.85, ani.height = 480*0.85)


# Gráfico compuesto QQplot + Histograma
saveGIF({
  for (i in 1:60) {
    res_ls1_frame = res_ls1 %>%
      filter(lambda == lamb_vec[[i]])
    
    pnorm = nor_ls1 %>% 
      filter(lambda==lamb_vec[[i]] & residuals == 'PWRES_y_1') %>% 
      magrittr::use_series(`p-value`) %>%  round(3)
    
    psimt = nor_ls1 %>% 
      filter(lambda==lamb_vec[[i]] & residuals == 'PWRES_y_1') %>% 
      magrittr::use_series(`p-value`) %>% round(3)
    
    p1 = res_ls1_frame %>% 
      normal_dist_plot(pwRes, pnorm=pnorm, psimt=psimt)
    
    p2 = res_ls1_frame %>% 
      hist.residual(pwRes, 'blue4') 
    
    pt = p1 + p2 
    
    print(pt)
  }
}, movie.name = 'hist_QQ_PWRES.gif', interval = 0.8, 
ani.width = 720, ani.height = 480/2)

#-------------------------------------------------------------------------------#
# Efecto en los parámetros estimados --------------------------------------------
#-------------------------------------------------------------------------------#
# Pre-alocación de lista de parámetros
param_ls <- list()

# Búsqueda en directorios de archivos de parámetros poblacionales estimadso
for (i in 1:60) {
  dis <- format(round(lamb_vec[[i]], 2), nsmall = 2)
  dir <- file.path('tbs_folder_1', paste0('lambda_', dis), 
                   'model_TBS', 'populationParameters.txt')
  param_ls[[i]] <- read_csv(dir)
}

# Gráfico de facet con valor de parámetro estimado por lambda
g1 <- param_ls %>% 
  map_dfr(~.x, .id = 'ID') %>%
  add_column(lambda = rep(lamb_vec, each = 11), 
             .before = 'value') %>% 
  mutate(parameter = 
           factor(parameter, 
                  levels = c('Cl_pop', 'V1_pop', 'Q_pop', 'V2_pop', 
                             
                             'omega_Cl','omega_V1','omega_Q','omega_V2',
                             'a', 'b', 'beta_Cl_tSCRMGDL'))) %>% 
  ggplot(aes(x=lambda, y=value, col=parameter)) +
  geom_line() +
  facet_wrap(~ parameter, ncol = 4, scales = 'free_y') +
  theme(legend.position = 'none')

# Almacenamiento de gráfico
ggsave('figures/4_parametros_TBS_1.pdf', plot = g1, 'pdf', 
       width = 7, height = 5)

ggsave('figures/5_parametros_TBS_2.pdf', 
       plot = g1 + geom_ribbon(aes(ymin = value - se_sa, ymax = value + se_sa,
                                   fill = parameter), alpha=0.2), 
       'pdf', width = 7, height = 5)


#-------------------------------------------------------------------------------#
# Ajuste de modelo/observaciones  --------------------------------------
#-------------------------------------------------------------------------------#
# Pre-alocación de **perfil_ls**
perfil_ls_1 <- list()
perfil_ls_2 <- list()

# Búsqueda de registros de verosimilitud
for (i in 1:60) {
  dis <- format(round(lamb_vec[[i]], 2), nsmall = 2)
  dir <- file.path('tbs_folder_1',paste0('lambda_',dis),
                   'model_TBS','ChartsData','IndividualFits')
  filename1 <- file.path(dir, 'y_1_fits.txt')
  filename2 <- file.path(dir, 'y_1_observations.txt')
  
  perfil_ls_1[[i]] <- read_csv(filename1)
  perfil_ls_2[[i]] <- read_csv(filename2)
}

perfil_ls_1 %<>% 
  map_dfr(~.x, .id = 'M') %>%
  mutate(M = as.numeric(M),
         lambda = magrittr::extract(lamb_vec, M)) 

perfil_ls_2 %<>%  
  map_dfr(~.x, .id = 'M') %>%
  mutate(M = as.numeric(M),
         lambda = magrittr::extract(lamb_vec, M))


#-------------------------------------------------------------------------------#
# Cálculo de TAD a partir de TSFD
last_dose <- read_delim('generador/data/1_data_TSFD.csv', ';',na = '.') %>% 
  group_by(ID) %>% 
  filter(EVID == 1) %>% 
  slice(n()) %>% 
  select(ID, TIME)

perfil_ls_1 %<>%
  left_join(last_dose, by = 'ID') %>%
  mutate(TAD = time - TIME) %>% 
  mutate(across(c(pop, popPred, indivPredMean, indivPredMode), 
                ~((.x*lambda) + 1)^(1/lambda)))

perfil_ls_2 %<>%
  left_join(last_dose, by = 'ID') %>%
  mutate(TAD = time - TIME) %>% 
  mutate(across(c(y_1, median, piLower, piUpper), 
                ~((.x*lambda) + 1)^(1/lambda)))


# Bondad de ajuste del modelo
saveGIF({
  for (i in 1:60) {
    obs1 = filter(perfil_ls_1, M == i)
    obs2 = filter(perfil_ls_2, M == i)
    
    p = ggplot(obs1, mapping = aes(x = TAD)) +
      geom_line(aes(y = pop), lty = 'dashed') +
      geom_line(aes(y = indivPredMean), lty = 'solid', col = 'green3') +
      facet_wrap( ~ ID, ncol = 4, labeller = labeller(.cols = label_both)) +
      geom_point(data = obs2,
                 mapping = aes(x = TAD, y = y_1), col = 'green4') +
      coord_cartesian(xlim = c(0, NA), ylim = c(0, 100)) + 
      labs(title = glue::glue('Efecto TBS en predicciones'),
           subtitle = glue::glue('Valor de lambda: {lamb_vec[[i]]}')) +
      xlab('TAD') + ylab(expression(C[P]))
    
    print(p)
  }
}, movie.name = 'perfil_OBS_PRED.gif', interval = 0.4, 
ani.width = 720*0.85, ani.height = 480*0.85)



