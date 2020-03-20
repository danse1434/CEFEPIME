##------------------------------------------------------------------------#
## Nombre del Script: Análisis de resultados de Bootstrap para modelo 
##  Base Refinado.
##  
## Proposito del Script:  Realizar un resumen de los resultados de bootstrap 
## para el modelo base 2 - refinado. 
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion:  17-03-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
##########################################################################-
# Introducción -----------------------------------------------------
##########################################################################-
# Carga de paquetes
require(rlang)
require(tidyverse)
require(Rcpp)
require(RcppArmadillo)

# Selección de directorio de trabajo
Equipo = 2

if (Equipo == 1) {
  homd <- file.path('C:', 'Users', 'Daniel', 'OneDrive', 'Documents', 
                    '(Proyecto)_Estudio_PKPD')
} else {
  homd <- file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD')
  
}

setwd(file.path(homd, 'CEFEPIME', '10. Base_Refin', 
                'ANALISIS'))

##########################################################################-
# Apertura de archivos de parámetros poblacionales ------------------------
##########################################################################-
# Esta sección tiene como objetivo leer todos los archivos de parámetros 
# poblacionales en todos los subdirectorios. 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Crear un objeto de tipo lista 
##  2 Almacenar en cada uno de sus elementos con una tabla leída de una 
##  dirección específica en el arbol de directorios.
##  3 Convertir el archivo lista en un data.frame con una columna para 
##  identificar de que carpeta proviene cada set de parámetros.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
param_list = list()

for (i in 1:1000) {
  param_list[[i]] = read_csv(
    file.path(
      '..',
      'BOOT',
      paste0('B', i),
      'BASE_MODEL',
      'populationParameters.txt'
    ),
    col_types = cols()
  )
}

param_df <- param_list %>%
  map_dfr( ~ as.data.frame(.x), .id = 'ID') %>% 
  as_tibble(.)

##########################################################################-
# Manipulación de archivo de datos leído ----------------------------------
##########################################################################-
# Creación de funciones complejas para cálculo de parámetros secundarios de 
# modelo de dos compartimento.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Crear función para el cálculo de microconstantes: se crea una función 
##  que toma el valor de las constantes de velocidad de transferencia, y 
##  realiza una resolución de una ecuación de segundo órden para calcular 
##  macroconstantes de decaímiento exponencial (micro.fun).
##  2. Crear función para el cálculo de macronostantes: se crea una función 
##  que toma el valor de las microconstantes y realiza cálculos para 
##  encontrar las macroconstantes propias de una función de decaimiento 
##  biexponencial (macro.fun). 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
micro.fun <- function(x, y, z) {
  alpha = (x + y + z + sqrt(((x + y + z) ^ 2) - (4 * x * z))) / 2
  beta =  (x + y + z - sqrt(((x + y + z) ^ 2) - (4 * x * z))) / 2
  return(list(alpha, beta))
}

macro.fun <- function(D, alpha, k21, V1, beta) {
  A = D * (alpha - k21) / (V1 * (alpha - beta))
  B = D * (k21 - beta) / (V1 * (alpha - beta))
  return(list(A, B))
}
##########################################################################-
# Manipulación de tabla con resultados de bootstrap no paramétrico
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Tomar el archivo de datos
##  2 Transformar el layout de los datos a una forma expandida con columnas 
##  para cada parámetro.
##  3 Calcular k10, k12, k21, alpha, beta, t_alpha, t_beta, A, y B. Con 
##  funciones definidas en situ, y funciones complejas micro.fun, y macro.fun
##  4 Transformar el layout de los datos a una forma compacta similar a la 
##  original.
##  5 Ordenar la tabla por ID.
##  6 Transformar la variable parameter en forma de un factor con ordenación 
##  específica.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
param_df1 <- param_df %>% 
  spread(., parameter, value) %>% 
  mutate(
    k10 = Cl_pop / V1_pop,
    k12 = Q_pop / V1_pop,
    k21 = Cl_pop / V2_pop,
    alpha = micro.fun(k10, k12, k21)[[1]],
    beta =  micro.fun(k10, k12, k21)[[2]],
    t_alpha = log(2) / alpha,
    t_beta = log(2) / beta,
    A = macro.fun(2000, alpha, k21, V1_pop, beta)[[1]],
    B = macro.fun(2000, alpha, k21, V1_pop, beta)[[2]]
  ) %>%
  gather(., -contains('ID'), key = 'parameter', value = 'value') %>%
  arrange(ID) %>% 
  mutate(parameter = 
           factor(parameter, 
                  level = c('V1_pop', 'V2_pop', 'Cl_pop', 'Q_pop',
                            'omega_V1', 'omega_V2', 'omega_Cl', 'omega_Q',
                            'a', 'b',
                            'k10', 'k12', 'k21', 'alpha', 'beta',
                            't_alpha', 't_beta', 'A', 'B')))

##########################################################################-
# Creación de gráficos-----------------------------------------------------
##########################################################################-
# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

##########################################################################-
# Crear una tabla con valores de media, mediana e intervalo de confianza del
# 95% para cada parámetro en param_df1.
lines <- param_df1 %>%
  group_by(parameter) %>%
  summarise(mn = mean(value),
            me = quantile(value, probs = 0.5),
            li = quantile(value, probs = (0.05 / 2)),
            ls = quantile(value, probs = 1 - (0.05 / 2)) )

##########################################################################-
# Crear un gráfico con las distribuciones de los parámetros en cada muestreo 
# de bootstrap, se colocan los valores de densidad, y estadísticos de 
# resumen (incluyendo los que se encuentran en la tabla *lines*)

G1 <- param_df1 %>%
  ggplot(mapping = aes(value), colour = 'black') +
  geom_histogram(aes(y = ..density..,  fill = parameter), bins = 34) +
  geom_density() + 
  geom_vline(lines, mapping = aes(xintercept = mn), colour = 'green2') +
  geom_vline(lines, mapping = aes(xintercept = me), colour = 'red3', lty = 'dashed') +
  geom_vline(lines, mapping = aes(xintercept = li), colour = 'red3', lty = 'dashed') +
  geom_vline(lines, mapping = aes(xintercept = ls), colour = 'red3', lty = 'dashed') +
  facet_wrap(. ~ parameter, ncol = 4,
             scales = 'free') +
  theme(legend.position = 'none') +
  xlab('') + ylab('Densidad')

  # Almacenamiento del gráfico
ggsave(filename = './FIGURAS/resultados_bootstrap.pdf',
       plot = G1,
       device = 'pdf',
       width = 11, height = 11)

##########################################################################-
# Crear un gráfico con las distribuciones para el set de parámetros estimados 
# por modelamiento. 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Utilizar los mismos niveles de param_df1 en param_df
##  2 Eliminar niveles no usados 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

param_df <- param_df %>%
  mutate(parameter = factor(parameter, levels = levels(param_df1$parameter))) %>%
  mutate(parameter = fct_drop(parameter))

G2 <- param_df %>%
  ggplot(mapping = aes(value), colour = 'black') +
  geom_histogram(aes(y = ..density..,  fill = parameter), bins = 34) +
  geom_density() + 
  facet_wrap(. ~ parameter, ncol = 4, scales = 'free') +
  theme(legend.position = 'none') +
  xlab('') + ylab('Densidad')

# Almacenamiento del gráfico
ggsave(filename = './FIGURAS/originales_bootstrap.pdf',
       plot = G2,
       device = 'pdf',
       width = 11, height = 11)

##########################################################################-
# Crear un grafico de cajas y bigotes con los valores estimados para cada 
# modelo. 

G3 <- param_df %>%
  ggplot(mapping = aes(x = parameter, y = value), colour = 'black') +
  geom_jitter(aes(colour = parameter), shape = '.') +
  geom_violin(aes(fill = parameter), alpha = 0.5) +
  xlab('') + ylab('Valor') +
  facet_wrap(. ~ parameter, ncol = 4, scales = 'free') +
  theme(legend.position = 'none', axis.title.x = element_blank())

# Almacenamiento del gráfico
ggsave(filename = './FIGURAS/bootstrap_dots.pdf',
       plot = G3,
       device = 'pdf',
       width = 6, height = 4)

##########################################################################-
# Almacenar la tabla como un archivo de CSV de tipo Excel -----------------
##########################################################################-
# Crear una tabla de parámetros con valores de resumen del análisis de bootstrap
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Tomar la tabla param_df1
##  2 Agrupar de acuerdo al tipo de parámetro.
##  3 Resumir de acuerdo a:
##    a Media del valor de cada parámetro
##    b Desviación estándar del valor de cada parámetro
##    c Coeficiente de variación del valor de cada parámetro
##    d Límite inferior IC95% del valor de cada parámetro
##    e Límite superior IC95% del valor de cada parámetro
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

param_table <- param_df1 %>%
  group_by(parameter) %>%
  summarise(
    mn = mean(value),
    sd = sd(value),
    rsd = sd(value) / mean(value),
    li = quantile(value, probs = (0.05 / 2)),
    ls = quantile(value, probs = 1 - (0.05 / 2))
  )  

# Almacentamiento de tabla de datos de CSV
param_table %>%
  write_csv(., path = './FIGURAS/Resultados_bootstrap.csv')


##########################################################################-
# Apertura de registro de convergencia de cada modelo ---------------------
##########################################################################-
# Esta sección tiene como objetivo leer todos los archivos de parámetros 
# poblacionales en todos los subdirectorios. 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Crear un objeto de tipo lista 
##  2 Almacenar en cada uno de sus elementos con una tabla leída de una 
##  dirección específica en el arbol de directorios.
##  3 Convertir el archivo lista en un data.frame con una columna para 
##  identificar de que carpeta proviene cada set de parámetros.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
convergence_list = list()

for (i in 1:1000) {
  convergence_list[[i]] = read_csv(
    file.path('..', 'BOOT', paste0('B', i), 'BASE_MODEL', 
              'ChartsData', 'Saem', 'CvParam.txt'),
    col_types = cols())
}

##########################################################################-
# Cálculo de convergencia de modelos
##########################################################################-
# Apertura de script de C++ con función para cálculo de modelo de regresión 
# con entrada de un vector de datos, posición, y lag.
Rcpp::sourceCpp('funcion_definitiva.cpp')
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Si la posición de cálculo es menor que el lag para la evaluación de 
##  convergencia, se realiza el cálculo de valor T de la pendiente de 
##  regresión. De otra manera, se regresa un NA numérico.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
slop_signC1 = function(A, j, l) {
  if (i > l) {
    return(slop_signC(A, j, l)[2])
  } else {
    return(NA_real_)
  }
}
##########################################################################-
# Aplicación de la función en los datos
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Creación de una lista
##  2 Con la lista de los data.frames de convergencia, se toma el j-ésimo,
##  se selecciona la columna del indicador de convergencia, y se extrae como 
##  un vector numérico.
##  3 Crear una matriz con dos columnas y un tamaño igual.
##  4 LLenar la primera columna con el valor de iteración, y llenar la segunda 
##  con el estadístico T de la pendiente calculado para cada iteración.
##  5 Almacenar en el subcomponente de la lista
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Z = list()

for (j in 1:1000) {
  A = convergence_list %>%
    magrittr::extract2(j) %>%
    select(convergenceIndicator) %>% 
    pull(.)
  
  Y = matrix(nrow = length(A), ncol = 2)
  
  for (i in 1:length(A)) {
    Y[i,1] = i
    Y[i,2] = slop_signC1(A, i, 100)
  }
  
  Z[[j]] = Y
}

##########################################################################-
# Convertir lista en data frame para manipulación
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar lista Z
##  2 Convertir lista en un data.frame con identificación de los datos 
##  remuestreados correspondientes.
##  3 Convertir el data frame en un tibble
##  4 Renombrar algunas columnas del df
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
Z <- Z %>% 
  map_dfr( ~ as.data.frame(.x), .id = 'ID') %>%
  as_tibble(.) %>% 
  rename(iteration = V1,
         estadistico_T = V2)

##########################################################################-
# Crear listas con convergencia de todos los parámetros
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar lista de datos de convergencia
##  2 Adicionar la columna con los valores de indicador de convergencia
##  3 Convertir en tibble
##  4 Calcular valor p
##  5 Agrupar por ID
##  6 Seleccionar puntos en los que se ha alcanzado estacionalidad
##  7 Seleccionar el primero punto 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
convergence_df <- 
  convergence_list %>%
  map_dfr(~ as.data.frame(.x), .id = 'ID') %>%
  left_join(., Z, by = c('ID', 'iteration')) %>% 
  as_tibble(.) %>% 
  mutate(pval = 2*pt(q = estadistico_T, df = 100-2, lower.tail = F),
         pval_log = if_else(pval <= 0.01, T, F))

convergence_df1 <- convergence_df %>% 
  group_by(ID) %>% 
  filter(pval_log == F) %>% 
  slice(1)
  
##########################################################################-
# Crear gráfico y almacenar con convergencia de indicador
GITER <-
convergence_df %>%
  ggplot(aes(x = iteration, y = convergenceIndicator, group = ID, col = phase)) +
  geom_line(alpha = 0.1) +
  geom_point(data = convergence_df1, shape = 4, colour = 'red') +
  labs(x = 'Iteraciones', y = 'Indicador Convergencia') +
  theme(panel.background = element_rect(colour = 'black'),
        legend.position = 'none')

# Almacenamiento del gráfico
ggsave(filename = './FIGURAS/trayectorias_Iter.pdf',
       plot = GITER,
       device = 'pdf',
       width = 11, height = 10)


  




















##########################################################################-
# Limpieza de corridas sin convergencia -----------------------------------
##########################################################################-
# El criterio de convergencia escogido fue que la pendiente de una regresión 
# del criterio de convergencia vs iteraciones. Para cada iteración se 
# calcula un ajuste de reg.lineal con las 100 anteriores iteraciones y se 
# observa la significancia de la pendiente. 
# 
# Todas las corridas de Bootstrap cumplieron este criterio durante la fase
# de exploración (N = 1000). Este criterio debería darse durante la fase de 
# exploración del algoritmo SAEM. Se eliminan aquellas que no convergen de 
# manera adecuada.
# 
# De manera ideal, si existe convergencia el criterio de convergencia debería
# presentarse con las iteraciones siguientes. Esto indica estabilidad a largo 
# plazo del algoritmo, y permite asegurar convergencia. 
# 
# Se crea una función para evaluar cuantas iteraciones deben evaluarse para 
# encontrar convergencia sin deriva a largo plazo *followed_function*. 

followed_function <- function(n) {
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ##  1 Seleccionar tabla
  ##  2 Agrupar por variable ID que representa la corrida
  ##  3 Eliminar las iteraciones que son mayores que 1000 ya que están en fase
  ##  de alisamiento.
  ##  4 Eliminar las iteraciones que tengan NA como valor p en un t-test con 
  ##  el valor de la pendiente de la eliminación, estas corresponden a las 100 
  ##  iteraciones iniciales.
  ##  5 Seleccionar sólo 5 columnas que se consideran como las más importantes
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  
  A = convergence_df %>%
    group_by(ID) %>%
    filter(iteration < 1000) %>% 
    filter(!is.na(pval)) %>% 
    select(ID, iteration, convergenceIndicator, estadistico_T, pval, pval_log)
  
  ##  6 Fitrar aquellas iteraciones con convergencia en las n-siguientes medidas 
  ##  por el criterio.
  for (i in 0:n) {
    A = A %>%
      filter(lead(pval_log, i) == FALSE)
  }
  ##  7 Con la matriz filtrada, seleccionar la primera iteración que cumple con el
  ##  requisito de convergencia en las n-siguiente iteraciones.
  ##  8 Desagrupar la tabla
  ##  9 Resumir el numero de corridas que se aceptan bajo el supuesto de n-siguientes
  ##  10  Extraer los resultados como un valor numérico simple
  Z <- A %>% 
    slice(1) %>% 
    ungroup(.) %>%
    summarise(count = n()) %>% 
    magrittr::extract2(1)
  
  return(Z)
}

##########################################################################-
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

plot(x = df_conv[,1], y = df_conv[,2], type = 'p', 
     xlab = 'Lag_Convergencia', ylab = 'Fracción de Iteraciones')

##########################################################################-
# Eliminación de corridas con problemas de convergencia
# 
# Evaluando el diagrama de iteraciones vs indicador de convergencia con marcación 
# de valores que cumplen el criterio de convergencia, se muestra que la mayoría 
# de corridas alcanzan convergencia entre 200 y 300 iteraciones, se deberían 
# dejar estas. Las corridas con problemas de convergencia tienen el indicador 
# positivo después de 500 iteraciones como opinión personal. El objetivo es 
# descartar estas corridas, que son aproximadamente 300 iteraciones. 

# Examinando df_conv se observa que tomando n = 20 lags se eliminan aproximadamente 
# 97 iteraciones, estas podrían corresponder a aquellas que tienen problemas de 
# convergencia, de manera arbitraria se escoge este número de lags.
# 

clean_df = convergence_df %>%
  group_by(ID) %>%
  filter(iteration < 1000) %>%
  filter(!is.na(pval))

for (i in 0:20) {
  clean_df = clean_df %>%
    filter(lead(pval_log, i) == FALSE)
}

clean_df <- clean_df %>%
  slice(1) %>%
  ungroup(.) %>%
  magrittr::extract2('ID') %>%
  unique() %>%
  as.double() %>%
  sort()


##########################################################################-
# Gráficos con corrección de convergencia ---------------------------------
##########################################################################-
# Manipulación de tabla con resultados de bootstrap no paramétrico
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Tomar el archivo de datos
##  2 Seleccionar sólo las corridas que hallan pasado el criterio de convergencia
##  3 Transformar el layout de los datos a una forma expandida con columnas 
##  para cada parámetro.
##  4 Calcular k10, k12, k21, alpha, beta, t_alpha, t_beta, A, y B. Con 
##  funciones definidas en situ, y funciones complejas micro.fun, y macro.fun
##  5 Transformar el layout de los datos a una forma compacta similar a la 
##  original.
##  6 Ordenar la tabla por ID.
##  7 Transformar la variable parameter en forma de un factor con ordenación 
##  específica.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
param_df2 <- param_df %>% 
  filter(ID %in% as.character(clean_df)) %>% 
  spread(., parameter, value) %>% 
  mutate(
    k10 = Cl_pop / V1_pop,
    k12 = Q_pop / V1_pop,
    k21 = Cl_pop / V2_pop,
    alpha = micro.fun(k10, k12, k21)[[1]],
    beta =  micro.fun(k10, k12, k21)[[2]],
    t_alpha = log(2) / alpha,
    t_beta = log(2) / beta,
    A = macro.fun(2000, alpha, k21, V1_pop, beta)[[1]],
    B = macro.fun(2000, alpha, k21, V1_pop, beta)[[2]]
  ) %>%
  gather(., -contains('ID'), key = 'parameter', value = 'value') %>%
  arrange(ID) %>% 
  mutate(parameter = 
           factor(parameter, 
                  level = c('V1_pop', 'V2_pop', 'Cl_pop', 'Q_pop',
                            'omega_V1', 'omega_V2', 'omega_Cl', 'omega_Q',
                            'a', 'b',
                            'k10', 'k12', 'k21', 'alpha', 'beta',
                            't_alpha', 't_beta', 'A', 'B')))

##########################################################################-
# Creación de gráficos-----------------------------------------------------
##########################################################################-
# Crear una tabla con valores de media, mediana e intervalo de confianza del
# 95% para cada parámetro en param_df2.
lines_2 <- param_df2 %>%
  group_by(parameter) %>%
  summarise(mn = mean(value),
            me = quantile(value, probs = 0.5),
            li = quantile(value, probs = (0.05 / 2)),
            ls = quantile(value, probs = 1 - (0.05 / 2)) )

##########################################################################-
# Crear un gráfico con las distribuciones de los parámetros en cada muestreo 
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
  theme(legend.position = 'none') +
  xlab('') + ylab('Densidad')

# Almacenamiento del gráfico
ggsave(filename = './FIGURAS/resultados_bootstrap_limpio.pdf',
       plot = G1b,
       device = 'pdf',
       width = 11, height = 11)


##########################################################################-
# Crear un grafico de cajas y bigotes con los valores estimados para cada 
# modelo. 

G3b <- 
  param_df2 %>%
  filter(!(parameter %in% c('Q_pop', 'A', 'B'))) %>% 
  ggplot(mapping = aes(x = parameter, y = value), colour = 'black') +
  geom_jitter(aes(colour = parameter), shape = '.') +
  geom_violin(aes(fill = parameter), alpha = 0.5) +
  xlab('') + ylab('Valor') +
  facet_wrap(. ~ parameter, ncol = 4, scales = 'free') +
  theme(legend.position = 'none', 
        axis.text.x = element_blank())

# Almacenamiento del gráfico
ggsave(filename = './FIGURAS/bootstrap_dots_limpia.pdf',
       plot = G3b,
       device = 'pdf',
       width = 6, height = 5)

##########################################################################-
# Almacenar la tabla como un archivo de CSV de tipo Excel -----------------
##########################################################################-
# Ajustar t0
# Se abre el archivo con los parámetros calculados con los datos originales, 
# esto corresponde a t0 en cada parámetro.
popParameters <-
  read_csv(file.path('..', 'BASE', 'BASE_MODEL', 'populationParameters.txt'))

# Se calculan los parámetros secundarios a partir de los datos de t0 con el 
# fin de incluir estos datos en el cálculo de intervalos de confianza. 
popParameters <- popParameters %>%  
  select(parameter, value) %>% 
  spread(., parameter, value) %>% 
  mutate(
    k10 = Cl_pop / V1_pop,
    k12 = Q_pop / V1_pop,
    k21 = Cl_pop / V2_pop,
    alpha = micro.fun(k10, k12, k21)[[1]],
    beta =  micro.fun(k10, k12, k21)[[2]],
    t_alpha = log(2) / alpha,
    t_beta = log(2) / beta,
    A = macro.fun(2000, alpha, k21, V1_pop, beta)[[1]],
    B = macro.fun(2000, alpha, k21, V1_pop, beta)[[2]]
  ) %>%
  gather(., key = 'parameter', value = 'value')

##########################################################################-
# Crear una tabla de parámetros con valores de resumen del análisis de bootstrap
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Tomar la tabla param_df2
##  2 Agrupar de acuerdo al tipo de parámetro.
##  3 Resumir de acuerdo a:
##    a Media del valor de cada parámetro
##    b Desviación estándar del valor de cada parámetro
##    c Coeficiente de variación del valor de cada parámetro
##    d Límite inferior IC95% del valor de cada parámetro
##    e Límite superior IC95% del valor de cada parámetro
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

param_table <- param_df2 %>%
  group_by(parameter) %>%
  summarise(
    mn = mean(value),
    sd = sd(value),
    rsd = sd(value) / mean(value),
    li = quantile(value, probs = (0.05 / 2)),
    ls = quantile(value, probs = 1 - (0.05 / 2))
  )  

##########################################################################-
# Cálculo de Intervalos de Confianza del 95% (Método Pivote)
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Incluir en el cálculo el valor de t0, valor estimado a partir del set 
##  de datos individuales.
##  2 Renombrar el valor pegado de la tabla como t0.
##  3 Calcular LI y LS a partir del método PIVOTE.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

param_table <- param_table %>%
  inner_join(popParameters, by = c("parameter")) %>%
  rename(t0 = value) %>%
  mutate(LI = 2 * t0 - ls,
         LS = 2 * t0 - li)

# Almacentamiento de tabla de datos de CSV
param_table %>%
  write_csv(., path = './FIGURAS/Resultados_bootstrap.csv')



