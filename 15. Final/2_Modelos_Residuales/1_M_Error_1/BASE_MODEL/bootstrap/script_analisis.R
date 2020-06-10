##------------------------------------------------------------------------------#
## Nombre del Script: Análisis de resultados de Bootstrap para modelo -----------
## final TBS con error aditivo. 
##  
## Propósito del Script:  Realizar un resumen de los resultados de bootstrap 
## para el modelo base 2 - refinado. 
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  17-03-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------#
# Introducción -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Carga de paquetes
require(magrittr)
require(glue)
require(rlang)
require(tidyverse)
require(Rcpp)
require(RcppArmadillo)

# Cambio en directorio de trabajo
setwd(file.path(getwd(), 'tbs_folder_2'))

#-------------------------------------------------------------------------------#
# Abrir archivo de parámetros originales (t0) -----------------------------------
#-------------------------------------------------------------------------------#
# Ajustar t0
# Se abre el archivo con los parámetros calculados con los datos originales, 
# esto corresponde a t0 en cada parámetro.

popParameters <- read_csv('lambda_-0.20/model_TBS/populationParameters.txt')

# Se calculan los parámetros secundarios a partir de los datos de t0 con el 
# fin de incluir estos datos en el cálculo de intervalos de confianza. 
popParameters1 <-  popParameters %>% 
  select(parameter, value) %>% 
  pivot_wider(names_from = parameter, values_from = value) %>% 
  mutate(
    List = pmap(list(2000, Cl_pop, Q_pop, V1_pop, V2_pop), constants_fun)
  ) %>%
  unnest(cols = c(List)) %>%
  pivot_longer(-contains('ID'), names_to = 'parameter', values_to = 'value') %>% 
  mutate(parameter = 
           factor(parameter, 
                  level = c('V1_pop', 'V2_pop', 'Cl_pop', 'Q_pop','beta_Cl_tSCRMGDL',
                            'omega_V1', 'omega_V2', 'omega_Cl', 'omega_Q',
                            'a', 'b', 'k10', 'k12', 'k21', 'alpha', 'beta',
                            't_alpha', 't_beta', 'A', 'B')))


#-------------------------------------------------------------------------------#
# Apertura de archivos de parámetros poblacionales ------------------------
#-------------------------------------------------------------------------------#
# Esta sección tiene como objetivo leer todos los archivos de parámetros 
# poblacionales en todos los subdirectorios. 
#................................................................................
#  1 Crear un objeto de tipo lista 
#  2 Almacenar en cada uno de sus elementos con una tabla leída de una 
#  dirección específica en el árbol de directorios.
#  3 Convertir el archivo lista en un data.frame con una columna para 
#  identificar de que carpeta proviene cada set de par?metros.
#................................................................................

param_list = list()

for (i in 1:1000) {
  param_list[[i]] = read_csv(glue('boot/B{i}/final_model/populationParameters.txt'),
                             col_types = cols())
}

param_df <- param_list %>%
  map_dfr( ~ as_tibble(.x), .id = 'ID')

#-------------------------------------------------------------------------------#
# Manipulación de archivo de datos leído ----------------------------------
#-------------------------------------------------------------------------------#
# Creación de funciones complejas para cálculo de parámetros secundarios de 
# modelo de dos compartimento.
#................................................................................
#  1 Crear función para el cálculo de microconstantes: se crea una función 
#  que toma el valor de las constantes de velocidad de transferencia, y 
#  realiza una resolución de una ecuación de segundo orden para calcular 
#  macroconstantes de decaimiento exponencial (micro.fun).
#  2. Crear función para el cálculo de macroconstantes: se crea una función 
#  que toma el valor de las microconstantes y realiza cálculos para 
#  encontrar las macroconstantes propias de una función de decaimiento 
#  biexponencial (macro.fun). 
#................................................................................

Rcpp::sourceCpp('boot/funcion_definitiva.cpp')

#-------------------------------------------------------------------------------#
# Manipulación de tabla con resultados de bootstrap no paramétrico
#................................................................................
#  1 Tomar el archivo de datos
#  2 Transformar el layout de los datos a una forma expandida con columnas 
#  para cada parámetro.
#  3 Calcular k10, k12, k21, alpha, beta, t_alpha, t_beta, A, y B. Con 
#  funciones definidas en situ, y funciones complejas micro.fun, y macro.fun
#  4 Transformar el layout de los datos a una forma compacta similar a la 
#  original.
#  5 Ordenar la tabla por ID.
#  6 Transformar la variable parameter en forma de un factor con ordenación 
#  específica.
#................................................................................

param_df1 <- param_df %>%
  pivot_wider(names_from = parameter, values_from = value) %>%
  mutate(
    List = pmap(list(2000, Cl_pop, Q_pop, V1_pop, V2_pop), constants_fun)
  ) %>%
  unnest(cols = c(List)) %>%
  pivot_longer(-contains('ID'), names_to = 'parameter', values_to = 'value') %>%
  arrange(ID) %>% 
  mutate(parameter = 
     factor(parameter, 
          level = c('V1_pop', 'V2_pop', 'Cl_pop', 'Q_pop', 'beta_Cl_tSCRMGDL',
                'omega_V1', 'omega_V2', 'omega_Cl', 'omega_Q',
                'a', 'b', 'k10', 'k12', 'k21', 'alpha', 'beta',
                't_alpha', 't_beta', 'A', 'B')))

#-------------------------------------------------------------------------------#
# Creación de gráficos-----------------------------------------------------
#-------------------------------------------------------------------------------#
# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

#-------------------------------------------------------------------------------#
# Crear una tabla con valores de media, mediana e intervalo de confianza del
# 95% para cada parámetro en param_df1.
lines <- param_df1 %>%
  group_by(parameter) %>%
  summarise(mn = mean(value),
            me = quantile(value, probs = 0.5),
            li = quantile(value, probs = (0.05 / 2)),
            ls = quantile(value, probs = 1 - (0.05 / 2)) )

#-------------------------------------------------------------------------------#
# Crear un gráfico con las distribuciones de los parámetros en cada muestreo 
# de bootstrap, se colocan los valores de densidad, y estadísticos de 
# resumen (incluyendo los que se encuentran en la tabla *lines*)

G1 <- param_df1 %>%
  ggplot(mapping = aes(value), colour = 'black') +
  geom_histogram(aes(y = ..density..,  fill = parameter), bins = 12, col='black') +
  geom_density() + 
  geom_vline(data = popParameters1, aes(xintercept = value), col = 'green2') + 
  geom_vline(lines, mapping = aes(xintercept = me), col = 'red3', lty = 'dashed') +
  geom_vline(lines, mapping = aes(xintercept = li), col = 'red3', lty = 'dashed') +
  geom_vline(lines, mapping = aes(xintercept = ls), col = 'red3', lty = 'dashed') +
  facet_wrap(. ~ parameter, ncol = 4, scales = 'free') +
  theme(legend.position = 'none',
        axis.title.x = element_blank()) + ylab('Densidad')

# Almacenamiento del gráfico
ggsave('boot/figures/1_resultados_bootstrap.pdf', G1, 'pdf',
       width = 9, height = 8)

#-------------------------------------------------------------------------------#
# Crear un gráfico con las distribuciones para el set de parámetros estimados 
# por modelamiento.
#................................................................................
#  1 Utilizar los mismos niveles de param_df1 en param_df
#  2 Eliminar niveles no usados 
#................................................................................

param_df %<>%
  mutate(parameter = factor(parameter, levels = levels(param_df1$parameter)),
         parameter = fct_drop(parameter)) 

G2 <- param_df %>%
  ggplot(mapping = aes(value), colour = 'black') +
  geom_histogram(aes(y = ..density..,  fill = parameter), bins = 34, col='black') +
  geom_vline(data = popParameters %>% 
               mutate(parameter = factor(parameter, levels(param_df$parameter))), 
             aes(xintercept = value), col = 'green2') + 
  geom_density() + 
  facet_wrap(. ~ parameter, ncol = 4, scales = 'free') +
  theme(legend.position = 'none',
        axis.title.x = element_blank()) +
  ylab('Densidad')

# Almacenamiento del gráfico
ggsave('boot/figures/2_originales_bootstrap.pdf', G2, 'pdf',
       width = 8, height = 6)

#-------------------------------------------------------------------------------#
# Crear un gráfico de cajas y bigotes con los valores estimados para cada 
# modelo. 

G3 <- param_df %>%
  ggplot(mapping = aes(x = parameter, y = value), colour = 'black') +
  geom_jitter(aes(colour = parameter), shape = '.') +
  geom_violin(aes(fill = parameter), alpha = 0.5) +
  ylab('Valor') +
  facet_wrap(. ~ parameter, ncol = 4, scales = 'free') +
  theme(legend.position = 'none', 
        axis.title.x = element_blank())

# Almacenamiento del gráfico
ggsave('boot/figures/3_bootstrap_dots.pdf', G3, 'pdf',
       width = 6, height = 4)

#-------------------------------------------------------------------------------#
# Almacenar la tabla como un archivo de CSV de tipo Excel -----------------
#-------------------------------------------------------------------------------#
# Crear una tabla de parámetros con valores de resumen del análisis de bootstrap
#................................................................................
#  1 Tomar la tabla param_df1
#  2 Agrupar de acuerdo al tipo de par?metro.
#  3 Resumir de acuerdo a:
#    a Media del valor de cada par?metro
#    b Desviación estándar del valor de cada parámetro
#    c Coeficiente de variación del valor de cada parámetro
#    d Límite inferior IC95% del valor de cada parámetro
#    e Límite superior IC95% del valor de cada parámetro
#................................................................................

param_table <- param_df1 %>%
  group_by(parameter) %>%
  summarise(
    mn = mean(value),
    sd = sd(value),
    rsd = sd(value) / mean(value),
    li = quantile(value, probs = (0.05 / 2)),
    ls = quantile(value, probs = 1 - (0.05 / 2))
  )  

param_table_gt <- param_table %>%
  mutate(
    parameter = fct_recode(parameter,
     'V1' = 'V1_pop', 'V2' = 'V2_pop', 'Cl' = 'Cl_pop', 'Q' = 'Q_pop',
     'beta_Cl' = 'beta_Cl_tSCRMGDL')) %>%  
  gt::gt() %>%
  gt::fmt_number(columns = vars('mn', 'sd', 'rsd', 'li', 'ls'),
                 decimals = 2) %>%
  gt::cols_merge_range(col_begin = vars(li), col_end = vars(ls)) %>%
  gt::cols_label(
    parameter = 'Parámetro',
    mn = 'Media',
    sd = 'SE',
    rsd = 'RSD (%)',
    li = 'IC 95%'
  ) %>%
  gt::tab_header(gt::md('**Resultados bootstrap no parámétrico**'),
                 gt::md('**Modelo Final TBS con error aditivo**')) %>% 
  gt::tab_options(table.font.size = "smaller",
                  data_row.padding = gt::px(3)) %>% 
  gt::tab_footnote('Obtenido mediante bootstrap no paramétrico', 
                   gt::cells_column_labels(columns = vars(li))) %>% 
  gt::tab_footnote('Obtenido mediante matriz de info. de Fisher (FIM)', 
                   gt::cells_column_labels(columns = vars(sd, rsd)))


# Almacenamiento de tabla de datos de HTML
gt::gtsave(param_table_gt, '4_resultados_bootstrap.html', 
           glue('{getwd()}/boot/figures/'))

  #-------------------------------------------------------------------------------#
# Apertura de registro de convergencia de cada modelo ---------------------
#-------------------------------------------------------------------------------#
# Esta sección tiene como objetivo leer todos los archivos de par?metros 
# poblacionales en todos los subdirectorios. 
#................................................................................
#  1 Crear un objeto de tipo lista 
#  2 Almacenar en cada uno de sus elementos con una tabla le?da de una 
#  dirección específica en el árbol de directorios.
#  3 Convertir el archivo lista en un data.frame con una columna para 
#  identificar de que carpeta proviene cada set de par?metros.
#................................................................................
convergence_list = list()

for (i in 1:1000) {
  convergence_list[[i]] =
    read_csv(glue('boot/B{i}/final_model/ChartsData/Saem/CvParam.txt'),
             col_types = cols())
}

#-------------------------------------------------------------------------------#
# Cálculo de convergencia de modelos
#-------------------------------------------------------------------------------#
# Apertura de script de C++ con función para cálculo de modelo de regresión 
# con entrada de un vector de datos, posición, y lag.
#................................................................................
#  1 Si pos. de cálculo (j) es mayor que el lag (l) para la evaluación de 
#  convergencia, se realiza el cálculo de valor T de la pendiente de 
#  regresión. De otra manera, se regresa un NA numérico.
#................................................................................

slop_signC1 = function(A, j, l) {
  if (j > l) {
    return(slop_signC(A, j, l)[2])
  } else {
    return(NA_real_)
  }
}

#-------------------------------------------------------------------------------#
# Aplicación de la función en los datos
#................................................................................
# 1 Creación de una lista
# 2 Con la lista de los data.frames de convergencia, se toma el j-ésimo,
# se selecciona la columna del indicador de convergencia, y se extrae como 
# un vector numérico.
# 3 Crear una matriz con dos columnas y un tamaño igual.
# 4 LLenar la primera columna con el valor de iteración, y llenar la segunda 
# con el estadístico T de la pendiente calculado para cada iteración.
# 5 Almacenar en el subcomponente de la lista
#................................................................................

Z = list()

for (j in 1:1000) {
  A = convergence_list %>%
    magrittr::extract2(j) %>%
    select(convergenceIndicator) %>% 
    pull(.)
  
  Y = matrix(nrow = length(A), ncol = 2)
  
  for (i in 1:length(A)) {
    Y[i, 1] = i
    Y[i, 2] = slop_signC1(A, i, 100)
  }
  
  Z[[j]] = Y
}

#-------------------------------------------------------------------------------#
# Convertir lista en data.frame para manipulación
#................................................................................
#  1 Convertir lista en un tibble con identificación de la replicación Bootstrap 
#  3 Renombrar las columnas del data.frame
#................................................................................

Z %<>% 
  map_dfr( ~ as_tibble(.x), .id = 'B') %>%
  rename(iteration = V1, estadistico_T = V2)

#-------------------------------------------------------------------------------#
# Crear listas con convergencia de todos los parámetros
#................................................................................
#  1 Seleccionar lista de datos de convergencia
#  2 Adicionar la columna con los valores de indicador de convergencia
#  3 Convertir en tibble
#  4 Calcular valor p
#  5 Agrupar por ID
#  6 Seleccionar puntos en los que se ha alcanzado estacionalidad
#  7 Seleccionar el primero punto 
#................................................................................

convergence_df <- convergence_list %>%
  map_dfr(~ as.data.frame(.x), .id = 'B') %>%
  left_join(Z, by = c('B', 'iteration')) %>% 
  as_tibble(.) %>%
  mutate(pval = 2 * pt(estadistico_T, 100 - 2, 0, 0), 
         pval_log = if_else(pval <= 0.01, T, F))

convergence_df1 <- group_by(convergence_df, B) %>% 
  filter(pval_log == F) %>%  slice(1)

#-------------------------------------------------------------------------------#
# Crear gráfico y almacenar con convergencia de indicador
GITER <- convergence_df %>%
  ggplot(aes(x = iteration, y = convergenceIndicator, group = B, col = phase)) +
  geom_line(alpha = 0.1) +
  geom_point(data = convergence_df1, shape = 4, colour = 'red') +
  labs(x = 'Iteraciones', y = 'Indicador Convergencia') +
  theme(panel.background = element_rect(colour = 'black'),
        legend.position = 'none')

# Almacenamiento del gráfico
ggsave('boot/figures/5_trayectorias_Iter.pdf', GITER, 'pdf',
       width = 11, height = 10)

#-------------------------------------------------------------------------------#
# Eliminación de corridas sin convergencia---------------------------------------
#-------------------------------------------------------------------------------#
# El criterio de convergencia escogido fue que la pendiente de una regresión 
# del criterio de convergencia vs iteraciones. Para cada iteración se 
# calcula un ajuste de reg.lineal con las 100 anteriores iteraciones y se 
# observa la significancia de la pendiente. 
# 
# No todas las corridas de Bootstrap cumplieron este criterio durante la fase
# de exploración. Este criterio debería darse durante la fase de exploración del 
# algoritmo SAEM. Se eliminan aquellas que no convergen de manera adecuada.
# 
# De manera ideal, si existe convergencia el criterio de convergencia debería
# presentarse con las iteraciones siguientes. Esto indica estabilidad a largo 
# plazo del algoritmo, y permite asegurar convergencia. 
# 
# Se crea una función para evaluar cuantas iteraciones deben evaluarse para 
# encontrar convergencia sin deriva a largo plazo *followed_function*. 

followed_function <- function(n) {
  #................................................................................
  #  1 Agrupar por variable ID que representa la corrida
  #  2 Eliminar las iteraciones que son mayores que 1000 ya que están en fase
  #  de alisamiento.
  #  3 Eliminar las iteraciones que tengan NA como valor p en un t-test con 
  #  el valor de la pendiente de la eliminación, estas corresponden a las 100 
  #  iteraciones iniciales.
  #  4 Seleccionar sólo 5 columnas que se consideran como las más importantes
  #................................................................................

  A = convergence_df %>%
    group_by(B) %>%
    filter(iteration < 1000) %>% 
    filter(!is.na(pval)) %>% 
    select(B, iteration, convergenceIndicator, estadistico_T, pval, pval_log)
  
  #  6 Filtrar aquellas iteraciones con convergencia en las n-siguientes medidas 
  #  por el criterio.
  for (i in 0:n) {
    A %<>% filter(lead(pval_log, i) == FALSE)
  }
  #  7 Con la matriz filtrada, seleccionar la primera iteración que cumple con el
  #  requisito de convergencia en las n-siguiente iteraciones.
  #  8 Desagrupar la tabla
  #  9 Resumir el número de corridas que se aceptan bajo el supuesto de n-siguientes
  #  10  Extraer los resultados como un valor numérico simple
  Z <- slice(A, 1) %>%
    ungroup(.) %>%
    summarise(count = n()) %>%
    magrittr::extract2(1)
  
  return(Z)
}

#-------------------------------------------------------------------------------#
# Realizar un estudio del número adecuado de n valores para evaluación de 
# convergencia a larga distancia.
# Matriz de convergencia
df_conv = data.frame(x = 1:60)

# Aplicar la función teniendo en cuenta de 1 hasta 60 valores siguientes a
# evaluar, para la determinación del valor óptimo.
for (i in 1:60) {
  df_conv[i, 'y'] = followed_function(df_conv[i, 'x'])
}

# Graficar el número de lags de convergencia vs la fracción de iteraciones
# escogidas.

G_conv_filt <- as_tibble(df_conv) %>% 
  ggplot(aes(x, y)) +  geom_point(col='green4') + 
  theme_bw() + ylab('Fracción de Iteraciones') + 
  xlab('Lags de Convergencia') + coord_cartesian(ylim = c(0, 1E3)) + 
  labs(title = 'Resultados de Bootstrap',
subtitle = 'Fracción de iteraciones remanentes por filtro de convergencia') 

ggsave('boot/figures/6_filtro_convergencia.pdf', G_conv_filt, 'pdf', 
         width = 7, height = 5)

#-------------------------------------------------------------------------------#
# Eliminación de corridas con problemas de convergencia
# 
# Evaluando el diagrama de iteraciones vs indicador de convergencia con marcación 
# de valores que cumplen el criterio de convergencia, se muestra que la mayor?a 
# de corridas alcanzan convergencia entre 200 y 300 iteraciones, se deber?an 
# dejar estas. Las corridas con problemas de convergencia tienen el indicador 
# positivo después de 500 iteraciones como opinión personal. El objetivo es 
# descartar estas corridas, que son aproximadamente 300 iteraciones. 

# Examinando df_conv se observa que tomando n=5 lags se elimina aprox. 250 iteraciones, 
# estas son las que más cuentan con problemas de convergencia a largo plazo. 
# Se escoge este número de lags para contar con un número adecuado de corridas.   

clean_df = convergence_df %>%
  group_by(B) %>%
  filter(iteration < 1000) %>%
  filter(!is.na(pval))

for (i in 0:5) {
  clean_df = clean_df %>%
    filter(lead(pval_log, i) == FALSE)
}

clean_df <- slice(clean_df, 1) %>%
  ungroup(.) %>%
  magrittr::extract2('B') %>%
  unique() %>%
  as.double() %>%
  sort()

#-------------------------------------------------------------------------------#
# Gráficos con corrección de convergencia ---------------------------------
#-------------------------------------------------------------------------------#
# Manipulación de tabla con resultados de bootstrap no paramétrico
#................................................................................
#  1 Tomar el archivo de datos
#  2 Seleccionar sólo las corridas que hallan pasado el criterio de convergencia
#  3 Transformar el layout de los datos a una forma expandida con columnas 
#  para cada parámetro.
#  4 Calcular k10, k12, k21, alpha, beta, t_alpha, t_beta, A, y B. Con 
#  funciones definidas en situ, y funciones complejas micro.fun, y macro.fun
#  5 Transformar el layout de los datos a una forma compacta similar a la 
#  original.
#  6 Ordenar la tabla por ID.
#  7 Transformar la variable parameter en forma de un factor con ordenación 
#  específica.
#................................................................................

param_df2 <- param_df %>%
  filter(ID %in% as.character(clean_df)) %>% 
  pivot_wider(names_from = parameter, values_from=value) %>% 
  mutate(
    List = pmap(list(2000, Cl_pop, Q_pop, V1_pop, V2_pop), constants_fun)
  ) %>%
  unnest(cols = c(List)) %>%
  pivot_longer(-contains('ID'), names_to = 'parameter', values_to = 'value') %>% 
  arrange(ID) %>% 
  mutate(parameter = 
         factor(parameter, 
          level = c('V1_pop', 'V2_pop', 'Cl_pop', 'Q_pop','beta_Cl_tSCRMGDL',
                    'omega_V1', 'omega_V2', 'omega_Cl', 'omega_Q',
                    'a', 'b', 'k10', 'k12', 'k21', 'alpha', 'beta',
                    't_alpha', 't_beta', 'A', 'B')))

#-------------------------------------------------------------------------------#
# Creación de gráficos-----------------------------------------------------
#-------------------------------------------------------------------------------#
# Crear una tabla con valores de media, mediana e intervalo de confianza del
# 95% para cada parámetro en param_df2.

lines_2 <- param_df2 %>%
  group_by(parameter) %>%
  summarise(mn = mean(value),
            me = quantile(value, probs = 0.5),
            li = quantile(value, probs = (0.05 / 2)),
            ls = quantile(value, probs = 1 - (0.05 / 2)) )

#-------------------------------------------------------------------------------#
# Crear un gráfico con las distribuciones de los par?metros en cada muestreo 
# de bootstrap, se colocan los valores de densidad, y estadísticos de 
# resumen (incluyendo los que se encuentran en la tabla *lines*)

G1b <- param_df2 %>%
  ggplot(mapping = aes(value), colour = 'black') +
  geom_histogram(aes(y = ..density..,  fill = parameter), bins = 34) +
  geom_density() + 
  geom_vline(lines_2, mapping = aes(xintercept = mn), colour = 'green2') +
  geom_vline(lines_2, mapping = aes(xintercept = me), colour = 'red3', lty = 'dashed') +
  geom_vline(lines_2, mapping = aes(xintercept = li), colour = 'red3', lty = 'dashed') +
  geom_vline(lines_2, mapping = aes(xintercept = ls), colour = 'red3', lty = 'dashed') +
  facet_wrap(. ~ parameter, ncol = 4,
             scales = 'free') +
  theme(legend.position = 'none', 
        axis.title.x = element_blank()) +
  ylab('Densidad')

# Almacenamiento del gráfico
ggsave('boot/figures/8_resultados_bootstrap_limpio.pdf', G1b, 'pdf',
       width = 8, height = 6)

#-------------------------------------------------------------------------------#
# Crear un gráfico de cajas y bigotes con los valores estimados para cada 
# modelo. 

G3b <- param_df2 %>%
  filter(!(parameter %in% c('Q_pop', 'A', 'B'))) %>% 
  ggplot(mapping = aes(x = parameter, y = value), colour = 'black') +
  geom_jitter(aes(colour = parameter), shape = '.') +
  geom_violin(aes(fill = parameter), alpha = 0.5) +
  xlab('') + ylab('Valor') +
  facet_wrap(. ~ parameter, ncol = 4, scales = 'free') +
  theme(legend.position = 'none', 
        axis.text.x = element_blank())

# Almacenamiento del gráfico
ggsave('boot/figures/9_bootstrap_dots_limpia.pdf', G3b, 'pdf',
       width = 6, height = 4)

#-------------------------------------------------------------------------------#
# Crear una tabla de parámetros con valores de resumen del análisis de bootstrap
#................................................................................
# 1 Seleccionar las columnas "parameter" y "value" de *param_df2*, este contiene 
# todos los replicados del Bootstrap.
# 2 Agrupar y anidar por "parameter", esto deja a los replicados como un vector 
# dentro de cada parámetro.
# 3 Calcular media, y mediana de los replicados mediante un mapeo.
# 4 Aplicar la función "confints" en C++ que calcula varios tipos de intervalos 
# de confianza, se aplica nivel de significancia de 0.05.
# 5 Calcular mediana de replicados.
# 6 Desanidar a B (que es el data.frame con los IC95%), al hacer esto se vuelven 
# columnas.
# 7 Unir a la tabla con parámetros estimados originales, se y rse (estimados con 
# Matriz de Información de Fisher).
# 8 Mover las columnas de valor estimado, se, y rse al frente.
#................................................................................

param_df3 <- param_df2 %>% 
  select(parameter, value) %>% 
  nest_by(parameter) %>%
  left_join(popParameters1, by = 'parameter') %>% 
  mutate(B = pmap(list(param=data, t0=value, alpha=0.05), confints)) %>% 
  mutate(median = map_dbl(data, ~ median(.x)),
         n      = map_dbl(data, ~length(.x))) %>% 
  unnest(B) %>% 
  left_join(popParameters, by = 'parameter') %>%
  relocate(c(value.x, se_sa, rse_sa), .before=data)

#-------------------------------------------------------------------------------#
# Reportar Resultados de Bootstrap en formato HTML ------------------------------
#------------------------------------------------------------------------------#

param_df3_gt <- param_df3 %>%
  select(-data,-value.y, -n) %>% 
  ungroup() %>% 
    mutate(parameter = fct_recode(parameter,'V1'='V1_pop','V2'='V2_pop',
                                  'Cl'='Cl_pop', 'Q'='Q_pop',
                                  'beta_Cl'='beta_Cl_tSCRMGDL')) %>%
  rename(value = value.x) %>% 
  gt::gt() %>%
  gt::fmt_number(columns = vars(value,se_sa,rse_sa,
                  median,classic.LI,classic.LS,percent.LI,percent.LS,
                  normal.LI,normal.LS,pivote.LI,pivote.LS,BCa.LI,
                  BCa.LS),
                 decimals = 2) %>% 
  gt::cols_merge_range(vars(classic.LI), vars(classic.LS), sep = I(", ")) %>% 
  gt::cols_merge_range(vars(percent.LI), vars(percent.LS), sep = I(", ")) %>% 
  gt::cols_merge_range(vars(normal.LI), vars(normal.LS), sep = I(", ")) %>% 
  gt::cols_merge_range(vars(pivote.LI), vars(pivote.LS), sep = I(", ")) %>% 
  gt::cols_merge_range(vars(BCa.LI), vars(BCa.LS), sep = I(", ")) %>% 
  gt::tab_spanner('Resultados Bootstrap (n = 800)', 
                  vars(median, classic.LI,classic.LS,percent.LI,percent.LS,
                       normal.LI,normal.LS,pivote.LI,pivote.LS,BCa.LI,
                       BCa.LS)) %>% 
  gt::cols_label(
    parameter = 'Parámetro',
    value = 'Estimado',
    se_sa = 'SE',
    rse_sa = 'RSD (%)',
    median = 'Mediana',
    classic.LI = 'Clásico',
    percent.LI = 'Percentil',
    normal.LI = 'Normal',
    pivote.LI = 'Pivote',
    BCa.LI = 'BCa'
  ) %>%
  gt::fmt_missing(1:15, 1:19, missing_text = '-') %>% 
  gt::tab_header(gt::md('**Resultados bootstrap no parámétrico**'),
                 gt::md('**Modelo Final TBS con error aditivo**')) %>% 
  gt::tab_options(table.font.size = "smaller",
                  data_row.padding = gt::px(3)) %>% 
  gt::tab_footnote('Obtenido mediante bootstrap no paramétrico', 
                   gt::cells_column_spanners(spanners = "Resultados Bootstrap (n = 800)")) %>% 
  gt::tab_footnote('Obtenido mediante matriz de info. de Fisher (FIM)', 
                   gt::cells_column_labels(columns = vars(se_sa, rse_sa)))


# Almacenamiento de tabla de datos de HTML
gt::gtsave(param_df3_gt, '10_IC_bootstrap.html', 
           glue('{getwd()}/boot/figures/'))




