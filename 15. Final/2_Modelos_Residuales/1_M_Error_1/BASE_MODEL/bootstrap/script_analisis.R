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
setwd(file.path(getwd(), 'BASE_MODEL', 'bootstrap', 'nonParametric'))

#-------------------------------------------------------------------------------#
# Función definitiva
# 
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
Rcpp::sourceCpp('../funcion_definitiva.cpp')

#-------------------------------------------------------------------------------#
# Abrir archivo de parámetros originales (t0) -----------------------------------
#-------------------------------------------------------------------------------#
# Ajustar t0
# Se abre el archivo con los parámetros calculados con los datos originales, 
# esto corresponde a t0 en cada parámetro.

popParameters <- read_csv('../../populationParameters.txt')

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
  param_list[[i]] = read_csv(
    glue(
      'M_Error_bootstrap_bootstrap_{i}/populationParameters.txt'
    ),
    col_types = cols()
  )
}

param_df <- param_list %>%
  map_dfr( ~ as_tibble(.x), .id = 'ID')

#-------------------------------------------------------------------------------#
# Manipulación de archivo de datos leído ----------------------------------
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
ggsave('../figures/1_resultados_bootstrap.pdf', G1, 'pdf',
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
ggsave('../figures/2_originales_bootstrap.pdf', G2, 'pdf', width = 8, height = 6)

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
ggsave('../figures/3_bootstrap_dots.pdf', G3, 'pdf', width = 6, height = 4)

#-------------------------------------------------------------------------------#
# Almacenar la tabla como un archivo de CSV de tipo Excel -----------------
#-------------------------------------------------------------------------------#
# Crear una tabla de parámetros con valores de resumen del análisis de bootstrap
#................................................................................
#  1 Tomar la tabla param_df1
#  2 Agrupar de acuerdo al tipo de parámetro.
#  3 Resumir de acuerdo a:
#    a Media del valor de cada parámetro
#    b Desviación estándar del valor de cada parámetro
#    c Coeficiente de variación del valor de cada parámetro
#    d Límite inferior IC95% del valor de cada parámetro
#    e Límite superior IC95% del valor de cada parámetro
#................................................................................

param_df3 <- param_df1 %>% 
  select(-ID) %>% 
  nest_by(parameter) %>% 
  left_join(popParameters1, by = 'parameter') %>%
  rename(t0 = value) %>% 
  mutate(median = map_dbl(data, ~ median(.x)),
         n      = map_dbl(data, ~length(.x))) %>% 
  mutate(B = pmap(list(param=data, t0=t0, alpha=0.05), confints)) %>% 
  unnest(B) %>% 
  left_join(popParameters, by = 'parameter') %>%
  relocate(c(value), .before=data)

#-------------------------------------------------------------------------------#
# Tabla presentación
param_df3_gt <- param_df3 %>%
  select(-data, -n) %>% 
  ungroup() %>% 
  mutate(parameter = fct_recode(parameter,'V1'='V1_pop','V2'='V2_pop',
                                'Cl'='Cl_pop', 'Q'='Q_pop',
                                'beta_Cl'='beta_Cl_tSCRMGDL')) %>%
  gt::gt() %>%
  gt::fmt_number(columns = vars(value,median,t0,classic.LI,classic.LS,percent.LI,percent.LS,
                                normal.LI,normal.LS,pivote.LI,pivote.LS,BCa.LI,
                                BCa.LS), decimals = 2) %>% 
  gt::cols_merge_range(vars(classic.LI), vars(classic.LS), sep = I(", ")) %>% 
  gt::cols_merge_range(vars(percent.LI), vars(percent.LS), sep = I(", ")) %>% 
  gt::cols_merge_range(vars(normal.LI), vars(normal.LS), sep = I(", ")) %>% 
  gt::cols_merge_range(vars(pivote.LI), vars(pivote.LS), sep = I(", ")) %>% 
  gt::cols_merge_range(vars(BCa.LI), vars(BCa.LS), sep = I(", ")) %>% 
  gt::tab_spanner('Resultados Bootstrap (n = 1000)', 
                  vars(median, classic.LI,classic.LS,percent.LI,percent.LS,
                       normal.LI,normal.LS,pivote.LI,pivote.LS,BCa.LI,
                       BCa.LS)) %>% 
  gt::cols_label(
    parameter = 'Parámetro',
    value = 'Estimado',
    median = 'Mediana',
    classic.LI = 'Clásico',
    percent.LI = 'Percentil',
    normal.LI = 'Normal',
    pivote.LI = 'Pivote',
    BCa.LI = 'BCa'
  ) %>%
  gt::fmt_missing(1:14, 1:19, missing_text = '-') %>% 
  gt::tab_header(gt::md('**Resultados bootstrap no parámétrico**'),
                 gt::md('**Modelo Final con error aditivo**')) %>% 
  gt::tab_options(table.font.size = "smaller",
                  data_row.padding = gt::px(3)) %>% 
  gt::tab_footnote('Obtenido mediante bootstrap no paramétrico', 
                   gt::cells_column_spanners(spanners = "Resultados Bootstrap (n = 1000)")) 

# Almacenamiento de tabla de datos de HTML
gt::gtsave(param_df3_gt, '4_resultados_bootstrap.html', glue('{getwd()}/../figures/'))
