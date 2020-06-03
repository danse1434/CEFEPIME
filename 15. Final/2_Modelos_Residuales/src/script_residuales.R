##------------------------------------------------------------------------------#
## Nombre del Script: Análisis de distribución de errores residuales ------------
##  
## Propósito del Script: Análisis de distribución de error residual con 
## diferentes modelos
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  06-02-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#

# Carga de paquetes
require(ggrepel)
require(patchwork)
require(readxl)
require(tidyverse)
require(rlang)

#-------------------------------------------------------------------------------#
# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

#-------------------------------------------------------------------------------#
# Gráfico de residuales vs TAD --------------------------------------------
#-------------------------------------------------------------------------------#
#................................................................................
#  1 Apertura de archivo de datos con los datos originales, y que tiene 
#  delimitador ";", este fue usado para realizar el modelamiento en sí mismo
#  2 Crear un archivo *last_dose* que contiene la última hora de adminis-
#  tración del antibiótico en el archivo *data* para cada individuo. 
#    2a Agrupar por ID
#    2b Filtrar los datos con EVID igual a 1 (sólo los eventos de 
#    administración)
#    2c Seleccionar los últimos en cada grupo (los datos ya están ordenados 
#    de menor a mayor)
#    2d Seleccionar s?lo las columnas ID y TIME
#  3 Modificar el archivo de residuales *y_1_residuals* al unir los datos 
#  de las últimas observaciones calculadas.
#................................................................................

data <- read_delim('data/1_data_TSFD.csv', ';', na = '.')


last_dose <-group_by(data, ID) %>% 
  filter(EVID == 1) %>% 
  slice(n()) %>% 
  select(ID, TIME)

#-------------------------------------------------------------------------------#
# Apertura de datos de residuales
pathres <- 'ChartsData/ScatterPlotOfTheResiduals/y_1_residuals.txt'

# Modelo de error aditivo
res_1 <- read_csv(file.path('1_M_Error_1', 'BASE_MODEL', pathres))
# Modelo de error proporcional
res_2 <- read_csv(file.path('2_M_Error_2', 'BASE_MODEL', pathres))
# Modelo de error combinado 1
res_3 <- read_csv(file.path('3_M_Error_3', 'BASE_MODEL', pathres))
# Modelo de error combinado 2 (original)
res_4 <- read_csv(file.path('..', '1_Modelo_Inicial', 'BASE_MODEL', pathres))

# Adicionar columna con TAD
TAD_add <- function(data) { 
  left_join(data, last_dose, by = 'ID') %>%
    mutate(TAD = time - TIME)
}

res_1 <- TAD_add(res_1)
res_2 <- TAD_add(res_2)
res_3 <- TAD_add(res_3)
res_4 <- TAD_add(res_4)


#-------------------------------------------------------------------------------#

RES_C <- function(data, x, y, xlab, ylab, 
                  colpoint = '#335E82', col_line = '#4682B4', 
                  colloess = '#252E36', ...) {
  # Convertir *y* en expresión para evaluación tardía
  y_quo <- rlang::ensym(y)
  #................................................................................
  #  Cálculo de percentiles empíricos
  #  1 Tomar la variable de datos a explorar.
  #  2 Crear una variable que tome 6 grupos a partir de los datos ordenados 
  #  con la función dplyr::ntile().
  #  3 Agrupar por la variable recién creado.
  #  4 Resumir por la media de TAD, y los percentiles P5%, P50%, y P95% de 
  #  los datos de residuales.
  #................................................................................
  A <- data %>%
    mutate(gr = ntile({{ x }}, 6)) %>% 
    group_by(gr) %>% 
    summarise(TIME = mean({{ x }}), 
              ME = quantile(x = !!y_quo, probs = 0.50), 
              LI = quantile(x = !!y_quo, probs = 0.05), 
              LS = quantile(x = !!y_quo, probs = 0.95))
  #-------------------------------------------------------------------------------#
  # Crear el gráfico de residuales con las especificaciones mostradas
  data %>% 
    ggplot(mapping = aes(x = {{ x }}, y = !!y_quo)) +
    geom_hline(yintercept = c(0,-1.25,1.25), lty = 'dashed', col='gray50') +
    geom_point(aes(col = !!y_quo < 2)) +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    xlab(xlab) + ylab(ylab) +
    scale_color_manual(values = 
                         setNames(c(colpoint, 'red'), c(T, F))) +
    stat_smooth(method = 'loess', formula = y ~ x, se = FALSE, 
                col = colloess, lty = 'solid', size = 0.9) +
    geom_line(data = A, aes(TIME, ME), 
              col = col_line, lty = 'solid', size = 0.9) +
    geom_line(data = A, aes(TIME, LI), 
              col = col_line, lty = 'solid', size = 0.9) +
    geom_line(data = A, aes(TIME, LS), 
              col = col_line, lty = 'solid', size = 0.9) + 
    geom_text_repel(aes(label = ifelse(!!y_quo <2, NA, ID))) +
    theme(legend.position = 'none') %>% 
    return(.)
}

G_RES_T_PWRES <- list()
G_RES_T_PWRES[[1]] <- RES_C(res_1, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')
G_RES_T_PWRES[[2]] <- RES_C(res_2, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')
G_RES_T_PWRES[[3]] <- RES_C(res_3, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')
G_RES_T_PWRES[[4]] <- RES_C(res_4, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')

G_RES_T_IWRES <- list()
G_RES_T_IWRES[[1]] <- RES_C(res_1, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')
G_RES_T_IWRES[[2]] <- RES_C(res_2, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')
G_RES_T_IWRES[[3]] <- RES_C(res_3, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')
G_RES_T_IWRES[[4]] <- RES_C(res_4, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')

G_RES_T_NPDE <- list()
G_RES_T_NPDE[[1]] <- RES_C(res_1, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')
G_RES_T_NPDE[[2]] <- RES_C(res_2, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')
G_RES_T_NPDE[[3]] <- RES_C(res_3, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')
G_RES_T_NPDE[[4]] <- RES_C(res_4, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')

G_RES_T_PWRES_TONJ <-
((G_RES_T_PWRES[[1]] + G_RES_T_PWRES[[2]]) /
 (G_RES_T_PWRES[[3]] + G_RES_T_PWRES[[4]])) + 
  plot_annotation(tag_levels = 'A')


G_RES_T_IWRES_TONJ <-
((G_RES_T_IWRES[[1]] + G_RES_T_IWRES[[2]]) /
 (G_RES_T_IWRES[[3]] + G_RES_T_IWRES[[4]])) + 
  plot_annotation(tag_levels = 'A')


G_RES_T_NPDE_CONJ <-
((G_RES_T_NPDE[[1]] + G_RES_T_NPDE[[2]]) /
(G_RES_T_NPDE[[3]] + G_RES_T_NPDE[[4]])) + 
plot_annotation(tag_levels = 'A')

ggsave('figures/T_PWRES.pdf', G_RES_T_PWRES_TONJ, 'pdf', 
       width = 7, height = 6)

ggsave('figures/T_IWRES.pdf', G_RES_T_IWRES_TONJ, 'pdf', 
       width = 7, height = 6)

ggsave('figures/T_NPDE.pdf', G_RES_T_NPDE_CONJ, 'pdf', 
       width = 7, height = 6)



#-------------------------------------------------------------------------------#
palls <- list(colpoint = '#CB92DD', col_line = '#B64DD6', colloess = '#4A1F57')

G_RES_C_PWRES <- list()
G_RES_C_PWRES[[1]] <- expr(RES_C(res_1, prediction_pwRes, pwRes, 'PRED', 'PWRES', 
                            !!!palls)) %>% eval()
G_RES_C_PWRES[[2]] <- expr(RES_C(res_2, prediction_pwRes, pwRes, 'PRED', 'PWRES', 
                                 !!!palls)) %>% eval()
G_RES_C_PWRES[[3]] <- expr(RES_C(res_3, prediction_pwRes, pwRes, 'PRED', 'PWRES', 
                                 !!!palls)) %>% eval()
G_RES_C_PWRES[[4]] <- expr(RES_C(res_4, prediction_pwRes, pwRes, 'PRED', 'PWRES', 
                                 !!!palls)) %>% eval()

G_RES_C_IWRES <- list()
G_RES_C_IWRES[[1]] <- expr(RES_C(res_1, prediction_iwRes_mean, iwRes_mean, 
                                 'PRED', 'IWRES', !!!palls)) %>% eval()
G_RES_C_IWRES[[2]] <- expr(RES_C(res_2, prediction_iwRes_mean, iwRes_mean, 
                                 'PRED', 'IWRES', !!!palls)) %>% eval()
G_RES_C_IWRES[[3]] <- expr(RES_C(res_3, prediction_iwRes_mean, iwRes_mean, 
                                 'PRED', 'IWRES', !!!palls)) %>% eval()
G_RES_C_IWRES[[4]] <- expr(RES_C(res_4, prediction_iwRes_mean, iwRes_mean, 
                                 'PRED', 'IWRES', !!!palls)) %>% eval()

G_RES_C_NPDE <- list()
G_RES_C_NPDE[[1]] <- expr(RES_C(res_1, prediction_npde, npde, 'PRED', 'NPDE',
                                  !!!palls)) %>% eval()
G_RES_C_NPDE[[2]] <- expr(RES_C(res_2, prediction_npde, npde, 'PRED', 'NPDE',
                                !!!palls)) %>% eval()
G_RES_C_NPDE[[3]] <- expr(RES_C(res_3, prediction_npde, npde, 'PRED', 'NPDE',
                                !!!palls)) %>% eval()
G_RES_C_NPDE[[4]] <- expr(RES_C(res_4, prediction_npde, npde, 'PRED', 'NPDE',
                                !!!palls)) %>% eval()

G_RES_C_PWRES_CONJ <-
  ((G_RES_C_PWRES[[1]] + G_RES_C_PWRES[[2]]) /
     (G_RES_C_PWRES[[3]] + G_RES_C_PWRES[[4]])) + 
  plot_annotation(tag_levels = 'A')


G_RES_C_IWRES_CONJ <-
  ((G_RES_C_IWRES[[1]] + G_RES_C_IWRES[[2]]) /
     (G_RES_C_IWRES[[3]] + G_RES_C_IWRES[[4]])) + 
  plot_annotation(tag_levels = 'A')


G_RES_C_NPDE_CONJ <-
  ((G_RES_C_NPDE[[1]] + G_RES_C_NPDE[[2]]) /
     (G_RES_C_NPDE[[3]] + G_RES_C_NPDE[[4]])) + 
  plot_annotation(tag_levels = 'A')

ggsave('figures/C_PWRES.pdf', G_RES_C_PWRES_CONJ, 'pdf', 
       width = 7, height = 6)
ggsave('figures/C_IWRES.pdf', G_RES_C_IWRES_CONJ, 'pdf', 
       width = 7, height = 6)
ggsave('figures/C_NPDE.pdf', G_RES_C_NPDE_CONJ, 'pdf', 
       width = 7, height = 6)


#-------------------------------------------------------------------------------#
# Pruebas de normalidad en los residuales ---------------------------------
#-------------------------------------------------------------------------------#
# A continuación, se realizan pruebas de normalidad univariadas para los 
# residuales en cada uno de los modelos.

#-------------------------------------------------------------------------------#
# Función que aplica una batería de pruebas de normalidad a los datos
## *data* ingrese el objeto de datos a analizar.
## *vector* ingrese la posición de las columnas donde se encuentran los 
## residuales de interés.
## *alpha* ingrese el valor de probabilidad para los tests. 

normtest_batery = function(data, vector){
  
  df = matrix(nrow = length(unique(vector)), ncol = 7)
  
  for (j in vector) {
    X = dplyr::pull(data,j) # Selecciona como un vector atómico a una columna 
    i = match(j,vector) # Encuentra la posición en el vector
    df[i,1] = colnames(data[,j])
    df[i,2] = shapiro.test(X)[["p.value"]]            # Shapiro-Wilks
    df[i,3] = nortest::ad.test(X)[["p.value"]]        # Anderson-Darling
    df[i,4] = nortest::cvm.test(X)[["p.value"]]       # Cramer von Mises
    df[i,5] = nortest::lillie.test(X)[["p.value"]]    # Liliefors
    df[i,6] = nortest::pearson.test(X)[["p.value"]]   # Pearson
    df[i,7] = nortest::sf.test(X)[["p.value"]]        # Shapiro Francia
  }
  colnames(df) = c('Variable','Shapiro', 'Anderson_Darling', 'Cramer_von_Mises',
                   'Liliefors','Pearson','Shapiro_Francia')
  df %>% 
    as_tibble(.) %>% 
    mutate(across(!matches('Variable'), ~as.numeric(.x))) %>% 
    return(.)
}

#-------------------------------------------------------------------------------#
# Aplicación de batería de pruebas de normalidad a la tabla
#................................................................................
#  1 Crear un objeto de tipo lista para almacenar los resultados de pruebas 
#  de normalidad.
#  2 Para todos los elementos de la lista se asigna uno de los data frames 
#  y se aplica la batería de pruebas de normalidad a los residuales que se 
#  encuentran en las posiciones 4 (pwRes), 7 (iwRes), y 13 (Npde) del 
#  dataframe.
#................................................................................

norm_resid_list = list()

for (i in 1:4) {
  norm_resid_list[[i]] = 
    normtest_batery(data = get(paste0("res_", i)),
                    vector = c(4, 7, 13))
}

norm_resid_list %>% 
  map_dfr( ~ as_tibble(.x), .id = 'ID') %>% 
  write_csv(., path = 'figures/resid_list.txt')

#-------------------------------------------------------------------------------#
#' Función para marcar valores que cumplen condición
#'
#' @param data objeto gt precursor 
#' @param var columna de tabla 
#' @param alpha nivel de significancia
#'
#' @return objeto GT con tabla formateada
#' @export
#'
#' @examples
#' marcar_pval_inf(., Shapiro)
#' 
marcar_pval_inf <- function(data, var, alpha = 0.05) {
  var_quo <- rlang::ensym(var)
  gt::text_transform(
    data = data,
    locations = gt::cells_body(
      columns = quo_name(var_quo), 
      rows = (!!var_quo >= alpha)), 
    fn = function(x) {paste0('<b><span style=', "\"color:#0003D1;\"",
                             'opacity: 0.95>', x, '</span></b>')}
    )}

# Especificación de valores de tabla

tabla_normalidad <- norm_resid_list %>% 
  map_dfr( ~ as_tibble(.x), .id = 'ID') %>%
  mutate(ID = glue::glue('Modelo {ID}'), 
    Variable = fct_recode(Variable,
                          PWRES = 'pwRes',
                          IWRES = 'iwRes_mean', 
                          NPDE  = 'npde')) %>% 
  gt::gt(groupname_col = 'ID') %>% 
  gt::fmt_scientific(
    columns = vars("Shapiro", "Anderson_Darling", "Cramer_von_Mises", 
                   "Liliefors", "Pearson", "Shapiro_Francia"), 
    decimals = 2) %>% 
  gt::cols_label(Variable         = 'Tipo de Residual',
                 Shapiro          = gt::md('_Shapiro_'),
                 Anderson_Darling = gt::md('_Anderson Darling_'),
                 Cramer_von_Mises = gt::md('_Cramer von Mises_'),
                 Liliefors        = gt::md('_Liliefors_'),
                 Pearson          = gt::md('_Pearson_'),
                 Shapiro_Francia  = gt::md('_Shapiro Francia_')) %>% 
  gt::tab_spanner(label = 'Prueba de Normalidad', 
                  columns = vars("Shapiro", "Anderson_Darling", "Cramer_von_Mises", 
                                 "Liliefors", "Pearson", "Shapiro_Francia")) %>% 
  gt::tab_header(title = gt::md('**Verificación de batería de pruebas 
                                de normalidad**')) %>% 
  marcar_pval_inf(., Shapiro) %>% 
  marcar_pval_inf(., Anderson_Darling) %>% 
  marcar_pval_inf(., Cramer_von_Mises) %>% 
  marcar_pval_inf(., Liliefors) %>% 
  marcar_pval_inf(., Pearson) %>% 
  marcar_pval_inf(., Shapiro_Francia) %>%  
  gt::tab_footnote(footnote = gt::html('Se reportan los valores <i>p</i> para cada 
                                     prueba de normalidad. Se marcan en azul los valores 
                                     en los que no se puede rechazar la H<sub>0</sub> de 
                                       normalidad.'),
                   locations = gt::cells_column_spanners('Prueba de Normalidad'))  %>%
  gt::tab_options(
    column_labels.font.size = "medium",
    table.font.size = "medium",
    data_row.padding = gt::px(3)
  )

#-------------------------------------------------------------------------------#
# Almacenamiento de tablas

tabla_normalidad %>% 
  gt::gtsave('tab_1.html',
             path = file.path(getwd(), 'figures', 'dist_residuales'))

tabla_normalidad %>% 
  gt::gtsave('tab_1.tex',
             path = file.path(getwd(), 'figures', 'dist_residuales'))

#-------------------------------------------------------------------------------#
# Gráficos QQ -------------------------------------------------------------------
#-------------------------------------------------------------------------------#
#' Función de gráficos QQ
#'
#' @param model número de modelo para formar la ref. a objeto res_{i}
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
normal_dist_plot <- function(model, param, colour) {
  param_quo <- rlang::ensym(param)
  if(is_missing(colour)){colour = 'blue4'}
  
  g1 <- get(glue::glue('res_{i}')) %>%
    ggplot(mapping = aes(sample = !!param_quo)) +
    stat_qq() +
    qqplotr::stat_qq_line(distribution = 'norm', col = colour) +
    qqplotr::stat_qq_band(distribution = 'norm', fill = colour, alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0) +
    scale_x_continuous(breaks = -3:+3, minor_breaks = seq(-3,3,0.5)) + 
    scale_y_continuous(breaks = -3:+3, minor_breaks = seq(-3,3,0.5)) + 
    coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) + 
    annotate(geom = 'label', x = -1.5, y = 2, 
             label = paste0('Modelo: (',letters[model], ')')) +
    xlab('Cuantiles teóricos') + ylab('Cuantiles de muestreo') + 
    theme(panel.grid.major = element_line(colour = 'gray90'),
          panel.grid.minor = element_line(colour = 'gray97'))
  
  return(g1)
}

#-------------------------------------------------------------------------------#
# Creación de gráficos QQ conjugados 

G_IWRES_CONJ <- (normal_dist_plot(1, 'iwRes_mean','red3') + 
                   normal_dist_plot(2, 'iwRes_mean', 'red3')) / 
  (normal_dist_plot(3, 'iwRes_mean', 'red3') + 
     normal_dist_plot(4, 'iwRes_mean', 'red3'))
  
G_PWRES_CONJ <- (normal_dist_plot(1, 'pwRes','blue3') + 
                   normal_dist_plot(2, 'pwRes', 'blue3')) / 
  (normal_dist_plot(3, 'pwRes', 'blue3') + 
     normal_dist_plot(4, 'pwRes', 'blue3'))

G_NPDE_CONJ <- (normal_dist_plot(1, 'npde','green3') + 
                   normal_dist_plot(2, 'npde', 'green3')) / 
  (normal_dist_plot(3, 'npde', 'green3') + 
     normal_dist_plot(4, 'npde', 'green3'))

# Almacenamiento
ggsave('figures/qqplot/IWRES_conj.pdf', G_IWRES_CONJ, 'pdf', 
       width = 8, height = 6)
ggsave('figures/qqplot/PWRES_conj.pdf', G_PWRES_CONJ, 'pdf', 
       width = 8, height = 6)
ggsave('figures/qqplot/NPDE_conj.pdf', G_NPDE_CONJ, 'pdf', 
       width = 8, height = 6)

#-------------------------------------------------------------------------------#
# Apertura de archivos con forma de la distribución -----------------------
#-------------------------------------------------------------------------------#
pathdis <- 'ChartsData/DistributionOfTheResiduals/theoreticalGuides.txt'

resid_guides = list()

# Modelo de error aditivo
resid_guides[[1]] <- read_csv(file.path('1_M_Error_1', 'BASE_MODEL', pathdis))
# Modelo de error proporcional
resid_guides[[2]] <- read_csv(file.path('2_M_Error_2', 'BASE_MODEL', pathdis))
# Modelo de error combinado 1
resid_guides[[3]] <- read_csv(file.path('3_M_Error_3', 'BASE_MODEL', pathdis))
# Modelo de error combinado 2 (original)
resid_guides[[4]] <- read_csv(file.path('..', '1_Modelo_Inicial', 'BASE_MODEL', 
                                        pathdis))
# 
resid_guides1 <- resid_guides %>% 
  map_dfr(~.x, .id = 'Modelo') %>% 
  mutate(Modelo = glue::glue('Modelo {Modelo}'))

resid_guides1 %>% 
ggplot(aes(x = abscissa, y = pdf, col = Modelo)) +
  geom_line()

#-------------------------------------------------------------------------------#
# 
pathdis1 <- 'ChartsData/DistributionOfTheResiduals/y_1_pdf.txt'

resid_pdf = list()

# Modelo de error aditivo
resid_pdf[[1]] <- read_csv(file.path('1_M_Error_1', 'BASE_MODEL', pathdis1))
# Modelo de error proporcional
resid_pdf[[2]] <- read_csv(file.path('2_M_Error_2', 'BASE_MODEL', pathdis1))
# Modelo de error combinado 1
resid_pdf[[3]] <- read_csv(file.path('3_M_Error_3', 'BASE_MODEL', pathdis1))
# Modelo de error combinado 2 (original)
resid_pdf[[4]] <- read_csv(file.path('..', '1_Modelo_Inicial', 'BASE_MODEL', 
                                        pathdis1))

resid_pdf1 <- resid_pdf %>% 
  map_dfr(~.x, .id = 'Modelo') %>% 
  mutate(Modelo = glue::glue('Modelo {Modelo}'))

resid_pdf1 %>% 
  ggplot(aes(color =Modelo)) +
  geom_bar(aes(x = pwRes_abscissa, y = pwRes_pdf, fill = Modelo), stat = 'identity') +
  geom_line(data = resid_guides1, aes(x = abscissa, y = pdf)) + 
  facet_wrap(.~Modelo, ncol=2)
  


#-------------------------------------------------------------------------------#
# 
pathdis2 <- 'ChartsData/DistributionOfTheResiduals/y_1_cdf.txt'

resid_cdf = list()

# Modelo de error aditivo
resid_cdf[[1]] <- read_csv(file.path('1_M_Error_1', 'BASE_MODEL', pathdis2))
# Modelo de error proporcional
resid_cdf[[2]] <- read_csv(file.path('2_M_Error_2', 'BASE_MODEL', pathdis2))
# Modelo de error combinado 1
resid_cdf[[3]] <- read_csv(file.path('3_M_Error_3', 'BASE_MODEL', pathdis2))
# Modelo de error combinado 2 (original)
resid_cdf[[4]] <- read_csv(file.path('..', '1_Modelo_Inicial', 'BASE_MODEL', 
                                     pathdis2))

resid_cdf1 <- resid_cdf %>%
  map_dfr(~.x, .id = 'Modelo') %>%
  mutate(Modelo = glue::glue('Modelo {Modelo}')) 

resid_cdf1 %>% 
  ggplot(aes(col = Modelo)) +
  geom_line(aes(x = pwRes_abscissa, y = pwRes_cdf)) +
  geom_line(data = resid_guides1, aes(x = abscissa, y = cdf)) + 
  facet_wrap(.~Modelo, ncol=2) + 
  scale_color_viridis_d() + 
  theme(panel.grid.major = element_line(colour = 'gray90'))


#-------------------------------------------------------------------------------#
# Adición de datos de residuales
pathdis3 <- 'ChartsData/ScatterPlotOfTheResiduals/y_1_residuals.txt'

resid_ls = list()

# Modelo de error aditivo
resid_ls[[1]] <- read_csv(file.path('1_M_Error_1', 'BASE_MODEL', pathdis3))
# Modelo de error proporcional
resid_ls[[2]] <- read_csv(file.path('2_M_Error_2', 'BASE_MODEL', pathdis3))
# Modelo de error combinado 1
resid_ls[[3]] <- read_csv(file.path('3_M_Error_3', 'BASE_MODEL', pathdis3))
# Modelo de error combinado 2 (original)
resid_ls[[4]] <- read_csv(file.path('..', '1_Modelo_Inicial', 'BASE_MODEL', 
                                     pathdis3))

resid_ls1 <- resid_ls %>%
  map_dfr(~.x, .id = 'Modelo') %>%
  mutate(Modelo = glue::glue('Modelo {Modelo}'))

#-------------------------------------------------------------------------------#
# Función para generar gráfico de distribución de residuales --------------
#-------------------------------------------------------------------------------#
# 
#' Función para simular distribución normal de -5 a 5
#'
#' @param mean media
#' @param sd desviación estándar
#'
#' @return tabla con valores X (abcsisas) y Y (ordenadas, función de densidad de 
#' probabilidad)
#' @export 
#' 
dnorm1 <- function(mean, sd) {
  v = seq(-5, 5, 0.01)
  tibble(x = v, y = dnorm(v, mean, sd))
}

#' Gráfico de residuales
#'
#' @param data tabla con valores para modelo
#' @param param parámetros a graficar
#' @param xlab etiqueta eje X
#'
#' @return 
#' Gráfico de ggplot2 con barras, línea de eje X, y distribución normal simulada
#' 
#' @export
#' @examples
#' resid_dist_plot(resid_ls1, pwRes, 'WRES')
#' 
resid_dist_plot <- function(data, param, xlab) {
  #................................................................................
  # 1 Agrupar data.frame de acuerdo a "Modelo"
  # 2 Resumir por media y desviación estándar del parámetro
  # 3 Crear columna-lista "Densidad" con tabla de valores de acuerdo a dnorm1
  # 4 Desanidar a "Densidad".
  #................................................................................

  Z = group_by(data, Modelo) %>%
    summarise(mn = mean({{ param }}), sd = sd({{ param }})) %>%  
    mutate(Densidad = pmap(list(mean = mn, sd = sd), dnorm1)) %>% 
    unnest(Densidad)
  # Graficar datos  
  g1 <- data %>% 
    ggplot(aes(x = {{ param }}, group = Modelo)) +
    geom_vline(xintercept = 0) +
    geom_histogram(aes(y = ..density.., fill = Modelo),
                   bins = 20, col = 'black', alpha = 0.9) +
    geom_line(Z, mapping = aes(x = x, y = y)) +
    facet_wrap(. ~ Modelo, ncol = 2) +
    coord_cartesian(xlim = c(-4, 4)) +
    ylab('Densidad') + xlab(xlab) + 
    scale_fill_viridis_d() +
    theme(legend.position = 'none')
  
  return(g1)
}

#-------------------------------------------------------------------------------#
# Almacenamiento de los gráficos en formato pdf
resid_dist_plot(resid_ls1, pwRes, 'WRES')
ggsave('figures/dist_residuales/Res_PWRES_Dist.pdf', device = 'pdf', 
       width = 5, height = 4)

resid_dist_plot(resid_ls1, iwRes_mean, 'IWRES')
ggsave('figures/dist_residuales/Res_IWRES_Dist.pdf', device = 'pdf', 
       width = 5, height = 4)

resid_dist_plot(resid_ls1, npde, 'NPDE')
ggsave('figures/dist_residuales/Res_NPDE_Dist.pdf', device = 'pdf', 
       width = 5, height = 4)

#-------------------------------------------------------------------------------#
# Gráfico de residuales -----------------------------------------------------
#-------------------------------------------------------------------------------#
#' Obtener residuales
#'
#' @param data datos de residuales 
#' @param var tipo de residual del cual se quiere obtener autocorrelación
#'
#' @return tibble con lags, acf, y pacf para residuales
#' @export
#'
#' @examples
#' extra_acf(res_3, 'pwRes')
#' 
extra_acf <- function(data, var) {
  
  data1 <- arrange(data, TAD) %>% 
    pull({{ var }})
  
  acf1 <- acf(data1, plot = FALSE)

  acf2 <- pacf(data1, plot = FALSE)

  tibble(lag = acf1$lag[1:20,,], 
         acf = acf1$acf[1:20,,],
         pacf = c(NA, pacf(res_2$pwRes, plot = FALSE)$acf[,,1]))
}

#-------------------------------------------------------------------------------#
#' Graficar la autocorrelación de residuales 
#'
#' @param data modelo con datos de autocorrelación ACF y PACF
#' @param var variable deseada puede ser ACF o PACF
#' @param tipo anotación de tipo de residual del cual parte gráfico
#'
#' @return
#' objeto ggplot2 con gráfico de correlación
#' @export
#'
#' @examples
#' acf.plot(extra_acf(res_3, 'pwRes') pacf)
#' 
acf.plot <- function(data, var, tipo) {
  data %>% 
    ggplot(aes(x=lag, y={{var}}, fill=lag)) +
    geom_hline(yintercept = 0) +
    # Acá suponemos que los datos iniciales son 90, cambiar a necesidad
    geom_hline(yintercept = c(+1.96/sqrt(90),
                              -1.96/sqrt(90)), lty='dashed',
               col='blue4') +
    geom_bar(stat = 'identity') +
    scale_fill_viridis_c() +
    annotate(geom = 'label', x = Inf, y = Inf, hjust = 1, vjust = 1,
             label = paste0('Tipo Residual: ',tipo)) +
    theme(legend.position = 'none')
}

#-------------------------------------------------------------------------------#
# Creación de una lista con gráficos de autocorrelación para cada modelo, tipo 
# de residual y tipo de autocorrelación (acf)

acfls <- list()

acfls[[1]] <- acf.plot(extra_acf(res_1, 'pwRes'), acf, 'pwRes')
acfls[[2]] <- acf.plot(extra_acf(res_1, 'pwRes'), pacf, 'pwRes')
acfls[[3]] <- acf.plot(extra_acf(res_1, 'iwRes_mean'), acf, 'iwRes_mean')
acfls[[4]] <- acf.plot(extra_acf(res_1, 'iwRes_mean'), pacf, 'iwRes_mean')
acfls[[5]] <- acf.plot(extra_acf(res_1, 'npde'), acf, 'npde')
acfls[[6]] <- acf.plot(extra_acf(res_1, 'npde'), pacf, 'npde')

acfls[[7]] <- acf.plot(extra_acf(res_2, 'pwRes'), acf, 'pwRes')
acfls[[8]] <- acf.plot(extra_acf(res_2, 'pwRes'), pacf, 'pwRes')
acfls[[9]] <- acf.plot(extra_acf(res_2, 'iwRes_mean'), acf, 'iwRes_mean')
acfls[[10]] <- acf.plot(extra_acf(res_2, 'iwRes_mean'), pacf, 'iwRes_mean')
acfls[[11]] <- acf.plot(extra_acf(res_2, 'npde'), acf, 'npde')
acfls[[12]] <- acf.plot(extra_acf(res_2, 'npde'), pacf, 'npde')


acfls[[13]] <- acf.plot(extra_acf(res_3, 'pwRes'), acf, 'pwRes')
acfls[[14]] <- acf.plot(extra_acf(res_3, 'pwRes'), pacf, 'pwRes')
acfls[[15]] <- acf.plot(extra_acf(res_3, 'iwRes_mean'), acf, 'iwRes_mean')
acfls[[16]] <- acf.plot(extra_acf(res_3, 'iwRes_mean'), pacf, 'iwRes_mean')
acfls[[17]] <- acf.plot(extra_acf(res_3, 'npde'), acf, 'npde')
acfls[[18]] <- acf.plot(extra_acf(res_3, 'npde'), pacf, 'npde')

acfls[[19]] <- acf.plot(extra_acf(res_4, 'pwRes'), acf, 'pwRes')
acfls[[20]] <- acf.plot(extra_acf(res_4, 'pwRes'), pacf, 'pwRes')
acfls[[21]] <- acf.plot(extra_acf(res_4, 'iwRes_mean'), acf, 'iwRes_mean')
acfls[[22]] <- acf.plot(extra_acf(res_4, 'iwRes_mean'), pacf, 'iwRes_mean')
acfls[[23]] <- acf.plot(extra_acf(res_4, 'npde'), acf, 'npde')
acfls[[24]] <- acf.plot(extra_acf(res_4, 'npde'), pacf, 'npde')

# Creación de gráficos compuestos
plot_resid_1 <- 
(acfls[[1]] + acfls[[2]]) /
(acfls[[3]] + acfls[[4]]) /
(acfls[[5]] + acfls[[6]])

plot_resid_2 <- 
  (acfls[[07]] + acfls[[08]]) /
  (acfls[[09]] + acfls[[10]]) /
  (acfls[[11]] + acfls[[12]])

plot_resid_3 <- 
  (acfls[[13]] + acfls[[14]]) /
  (acfls[[15]] + acfls[[16]]) /
  (acfls[[17]] + acfls[[18]])

plot_resid_4 <- 
  (acfls[[19]] + acfls[[20]]) /
  (acfls[[21]] + acfls[[22]]) /
  (acfls[[23]] + acfls[[24]])

# Almacenamiento
ggsave('figures/correlacion/modelo_1.pdf', plot_resid_1, 'pdf', 
       width = 7, height = 7)
ggsave('figures/correlacion/modelo_2.pdf', plot_resid_2, 'pdf', 
       width = 7, height = 7)
ggsave('figures/correlacion/modelo_3.pdf', plot_resid_3, 'pdf', 
       width = 7, height = 7)
ggsave('figures/correlacion/modelo_4.pdf', plot_resid_4, 'pdf', 
       width = 7, height = 7)


#-------------------------------------------------------------------------------#
# Boxplot de residuales ---------------------------------------------------------
#-------------------------------------------------------------------------------#
res_t <- list(res_1, res_2, res_3, res_4) %>% 
  map_dfr(~.x, .id = 'Modelo') %>% 
  mutate(Modelo = fct_relabel(Modelo, ~glue::glue('Modelo {.x}'))) 

  
#-------------------------------------------------------------------------------#
#' Función de obtención de boxplot de residuales
#'
#' @param data 
#' @param xvar variable eje X
#' @param yvar variable eje Y 
#' @param ylab etiqueta de eje Y
#'
#' @return
#' @export
#'
#' @examples
box.plot.1 <- function(data, xvar, yvar, ylab) {
  data %>% 
    ggplot(aes(x = {{ xvar }}, y = {{ yvar }})) +
    geom_violin() +
    stat_boxplot(col='black', geom = "errorbar", width = 0.2) +
    geom_boxplot(col='black', fill='white',
                 width=0.3, outlier.shape = NULL) +
    geom_point(aes(col = factor(ID)), 
               position = position_jitter(width = 0.1, height = 0.1)) +
    geom_hline(yintercept = 0, lty='solid') +
    scale_y_continuous(breaks = -3:3) +
    scale_colour_viridis_d() +
    ylab(ylab) + 
    coord_cartesian(ylim = c(-3, 3)) +
    theme(legend.position = 'none', 
          panel.grid.major.y = element_line(colour = 'gray70'),
          panel.grid.minor.y = element_line(colour = 'gray90'),
          axis.title.x = element_blank())
}

# Almacenamiento
box.plot.1(res_t, Modelo, pwRes, 'PWRES')
ggsave('figures/boxplot/box_PWRES.pdf', device = 'pdf', 
       width = 5, height = 4)

box.plot.1(res_t, Modelo, iwRes_mean, 'IWRES')
ggsave('figures/boxplot/box_IWRES.pdf', device = 'pdf', 
       width = 5, height = 4)

box.plot.1(res_t, Modelo, npde, 'NPDE')
ggsave('figures/boxplot/box_NPDE.pdf', device = 'pdf', 
       width = 5, height = 4)

3.3 * log2(90)


require(gganimate)

for (i in 1:20) {
  A <- res_t %>% 
    ggplot(aes(x = pwRes)) +
    geom_histogram(bins = i) + 
    facet_wrap(~ Modelo)
  ggsave(paste0('A',i,'.png'), A, device = 'png')
}


