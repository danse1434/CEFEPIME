##========================================================================#
## Nombre del Script: Análisis de influencia de valores iniciales con un 
##  diseño factorial completo ----------------------------------------------
##  
## Proposito del Script: analizar la influencia de los valores iniciales de 
##  efectos fijos del modelo, en la verosimilitud o los parámetros que se 
##  alcanzan. Los valores iniciales de efectos aleatorios fueron 1.0 para 
##  todos excepto, b que contaba con 0.3 de manera inicial. Las configuraciones 
##  del algoritmo SAEM: fase exploratoria 500/1000, fase alisamiento 
##  500/1000, anillamiento simulado 500 - condiciones similares a bootstrap.  
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion: 05-03-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##========================================================================#
# Carga de paquetes
require(rlang)
require(tidyverse)
# Definición de directorio principal
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '07. Minimizacion', '2_Efecto_Evaluacion'))

##########################################################################-
# Especificación de variación en condiciones iniciales --------------------
##########################################################################-
# Apertura de archivo de datos de parámetros de modelo base
populationParameters <-
  read_csv(
    "../../04. Residual/Modelo_8/RES_M1/populationParameters.txt"
  )

##########################################################################-
# Cálculo de valor mínimo (50%) y valor máximo (150%) respecto al valor 
# nominal estimado en el modelo base.
pop_par <- populationParameters %>% 
  mutate(minvalue = 0.5 * value,
         maxvalue = 1.5 * value)

##########################################################################-
# Creación de matriz con parámetros (efectos fijos) y sus perturbaciones.
A <- pop_par %>%
  filter(parameter %in% c('Cl_pop', 'V1_pop', 'Q_pop', 'V2_pop')) %>%
  select(parameter, value, minvalue, maxvalue) %>%
  t(.)
# Seleccionar encabezados a partir de primera fila, eliminar la primera 
# fila, y convertir la matriz resultante en numérica.
h <- A[1,]
A <- A[-1,]
A <- apply(A, 2, as.numeric)
# Asignación de nombres de columna a matriz
colnames(A) <- h
A # Visualizar la matriz creada

##########################################################################-
# Crear una expansión a modo de diseño factorial completo, se toman como 
# factores a cada una de las columnas.
B <- expand.grid(
  Factor1 = A[, 1],
  Factor2 = A[, 2],
  Factor3 = A[, 3],
  Factor4 = A[, 4]
)
# Asignación de nombres de columna a matriz expandida
colnames(B) <- h
B # Visualizar la matriz creada

# # NO CORRER AQUÍ
# ##########################################################################-
# # Creación de archivos individuales ---------------------------------------
# ##########################################################################-
# # Apertura del archivo de control
# ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# ##  1 Abrir archivo de texto de control de Monolix.
# ##  2 Modificar el archivo para ser leído como una string de R.
# ##  3 Reemplazar el nombre de archivo a leer por el control de Monolix por 
# ##  el correspondiente en cada carpeta.
# ##  Se cambian los valores de parámetro inicial, se conserva la misma semilla 
# ##  de simulación para cada iteración.
# ##  4 Crear carpetas para contener los archivos.
# ##  5 Almacenar el archivo en el directorio correspondiente con el nombre 
# ##  de "BASE_NEW.mlxtran", con formato de texto.
# ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# fileName <- 'MODELO_BASE_efecto.mlxtran'
# Z = readChar(fileName, file.info(fileName)$size)
# 
# for (i in 1:dim(B)[1]) {
#   dir.create(file.path('assessment', paste0('A', i)))
#   
#   Y <- Z %>%
#     str_replace("(?<=Cl\\_pop\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)",
#                 paste(B[i, 1])) %>%
#     str_replace("(?<=V1\\_pop\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)",
#                 paste(B[i, 2])) %>%
#     str_replace("(?<=Q\\_pop\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)",
#                 paste(B[i, 3])) %>%
#     str_replace("(?<=V2\\_pop\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)",
#                 paste(B[i, 4]))
#   
#   write_lines(Y,
#               file.path('assessment',
#                         paste0('A', i),
#                         'MODELO_BASE_efecto.mlxtran'),
#               sep = '\n')
#   
#   file.copy(
#     from = file.path('data.csv'),
#     to = file.path('assessment', paste0('A', i)))
# }

##########################################################################-
# Análisis de iteraciones -------------------------------------------------
##########################################################################-
# Creación de una lista vacía
LL_ls <- vector(mode = 'list', length = dim(B)[1])
# Colocar en cada item un data table abriendo cada carpeta de manera iterativa
for (i in 1:dim(B)[1]) {
  LL_ls[[i]] <-
    read_csv(
      file.path('assessment_2', paste0('A', i), 'MODELO_BASE_efecto', 
                'LogLikelihood', 'logLikelihood.txt'), col_types = cols())
}

# Cambio en la matriz con el diseño factorial para incluir el identificador 
# de iteración.
B <- B %>%
  rownames_to_column(var = 'Iteration')

##########################################################################-
# Cambiar la matriz DF
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Cambiar la lista a un data.frame con identificación de la iteration
##  2 Convertir en tibble
##  3 Unir la matriz de diseño factorial
##  4 Renombrar la tabla creada en la variable IS
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

LL_df <- LL_ls %>%
  map_dfr( ~ as.data.frame(.x), .id = 'Iteration') %>%
  as_tibble(.) %>%
  inner_join(B, by = "Iteration") %>% 
  rename(IS = importanceSampling)
  
##########################################################################-
# Aplicar un ANOVA con iteraciones en el criterio -2LL, las variables 
# dependientes son cuantitativas.
LL_df %>%
  filter(., criteria == "-2LL") %>%
  aov(IS ~ Cl_pop * V1_pop * Q_pop * V2_pop, data = .) %>% 
  summary()

##########################################################################-
# Aplicar un ANOVA sin iteraciones en el criterio -2LL, las variables 
# dependientes son cuantitativas.

LL_df %>%
  filter(., criteria == "-2LL") %>%
  aov(IS ~ Cl_pop + V1_pop + Q_pop + V2_pop, data = .) %>% 
  summary()

##########################################################################-
# Función que crea listas con especificaciones de variable
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Crear función interna:
##    1 Con parámetros P y X, crear un carácter con las partes especificadas 
##    de la matríz A
##  2 Crear una lista con los números generados con los nombres para pasar 
##  como argumento a la función fct_recode mediante evaluación tidy.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

fcrec <- function(x) {
  fcrc <- function(P, X) {
    paste(A[X, P])
  }
  
  list(L = fcrc(x, 2),
       M = fcrc(x, 1),
       S = fcrc(x, 3)) %>% return(.)
}

##########################################################################-
# Arreglar la matríz
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Cambiar las variables de parámetros por factores.
##  2 Recodificar los factores con el tipo de NIVEL para generar una variable 
##  cualtitativa nominal (L, M, H) - evaluación tidy.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

LL_df <- LL_df %>% 
  mutate(Cl_pop = factor(Cl_pop),
         V1_pop = factor(V1_pop),
         Q_pop = factor(Q_pop),
         V2_pop = factor(V2_pop)) %>% 
  mutate(Cl = fct_recode(Cl_pop, !!!fcrec(1)),
         V1 = fct_recode(V1_pop, !!!fcrec(2)),
          Q = fct_recode(Q_pop, !!!fcrec(3)), 
         V2 = fct_recode(V2_pop, !!!fcrec(4)) )

##########################################################################-
# Visualizar la matríz creada
LL_df
##########################################################################-
# Aplicar un ANOVA sin iteraciones en el criterio -2LL, no hay diferencias 
# al utilizar otros criterios. Las variables independientes son cualitativas.
LL_df %>%
  filter(., criteria == "-2LL") %>%
  aov(IS ~ Cl + V1 + Q + V2, data = .) %>% 
  summary(.)

##########################################################################-
# Lista de parámetros -----------------------------------------------------
##########################################################################-
Param_ls <- vector(mode = 'list', length = dim(B)[1])
# Directorio de tablas con parámetros

for (i in 1:dim(B)[1]) {
  Param_ls[[i]] <-
    read_csv(
      file.path('assessment_2', paste0('A', i), 'MODELO_BASE_efecto', 
                'populationParameters.txt'), col_types = cols())
}

##########################################################################-
# Cambiar la matriz DF
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Cambiar la lista a un data.frame con identificación de la iteration
##  2 Convertir en tibble
##  3 Unir la matriz de diseño factorial
##  4 Convertir la variable parameter en un factor ordenado
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Param_df <- Param_ls %>% 
  map_dfr( ~ as.data.frame(.x), .id = 'Iteration') %>%
  as_tibble(.) %>%
  inner_join(B, by = "Iteration") %>% 
  mutate(parameter = factor(parameter, 
                            levels = c('Cl_pop', 'V1_pop', 'Q_pop', 'V2_pop',
                                       'omega_Cl', 'omega_V1', 'omega_Q', 'omega_V2',
                                       'a', 'b')))

##########################################################################-
# Creación de gráfico con iteraciones
# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

##########################################################################-
# Resumen de medianas para cada parámetro en específico.
# Agrupar por parámetro; resumir por mediana de valor en todas las iteraciones
Param_DF <- Param_df %>%
  group_by(parameter) %>%
  summarise(mn = mean(value))

# Crear gráfico de comparación de parámetros en las corridas, los parámetros 
# tienen barras de error correspondientes a SE.
G1 <- Param_df %>%
  mutate(Iteration = as.double(Iteration)) %>%
  ggplot(aes(x = Iteration, col = factor(Iteration))) +
  geom_hline(data = Param_DF, mapping = aes(yintercept = mn), col = 'gray') +
  geom_point(aes(y = value), shape = 20, colour = "#1c86ee") +
  geom_errorbar(aes(ymin = value - se_sa, ymax = value + se_sa), colour = "#1c86ee") +
  facet_wrap( ~ parameter, ncol =  4, scales = "free") +
  theme(legend.position = "none") +
  xlab('Corrida') + ylab('Valor')

G1

# Almacenamiento de gráfico
ggsave(G1, filename = 'G1.pdf', device = 'pdf', width = 8, height = 6)

##########################################################################-
# Lista de datos de convergencia ------------------------------------------
##########################################################################-
# Creación de lista vacía pre-alocada
Conver_ls <- vector(mode = 'list', length = dim(B)[1])

# Apertura de tablas en cada uno de los directorios de manera iterativa
for (i in 1:dim(B)[1]) {
  Conver_ls[[i]] <-
    read_csv(
      file.path('assessment_2', paste0('A', i), 'MODELO_BASE_efecto', 
                'ChartsData', 'Saem', 'CvParam.txt'), col_types = cols())
}

##########################################################################-
# Creación y modificación de tabla con datos de convergencia
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Selección de lista
##  2 Creación de tabla única con identificación de corrida
##  3 Conversión de tibble
##  4 Unión con tabla de identificación de diseño factorial
##  5 Conversión en datos colapsados, para los parámetros del modelo además 
##  de criterio de convergencia.
##  6 Cambio del nombre de parámetro, eliminando el sufijo .x
##  7 Cambiar la variable parameter a un factor ordenado
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Conver_df <- Conver_ls %>% 
  map_dfr( ~ as.data.frame(.x), .id = 'Run') %>%
  as_tibble(.) %>%
  inner_join(B, by = c("Run" = "Iteration")) %>% 
  gather(matches("\\_pop\\.x|omega\\_|a$|b$|conv"), key = "Parameter", value = "Value") %>% 
  mutate(Parameter = str_replace(Parameter, "\\_pop\\.x", "_pop")) %>% 
  mutate(Parameter = factor(Parameter, 
                            levels = c('Cl_pop', 'V1_pop', 'Q_pop', 'V2_pop',
                                       'omega_Cl', 'omega_V1', 'omega_Q', 'omega_V2',
                                       'a', 'b', 'convergenceIndicator')))


##########################################################################-
# Creación de gráfico con trayectorias de convergencia

# Superposición de parámetros estimados de modelo base, se cambian nombres 
# y variables para ajustarse al acople.
pop_par = pop_par %>%
  rename(Parameter = parameter) %>%
  mutate(Parameter = factor(Parameter, levels = levels(Conver_df$Parameter)))

G2 <- Conver_df %>% 
  filter((iteration + 10 + 1) %% 10 == 0) %>%
  mutate(parameter = Parameter) %>% 
  ggplot(aes(x = iteration, y = Value, group = Run, col = Run)) +
  geom_line(alpha = 0.5, col = '#1c86ee') +
  geom_hline(data = pop_par, mapping = aes(yintercept = value), lty = 'dashed') +
  geom_hline(data = pop_par, mapping = aes(yintercept = value - se_sa), 
             lty = 'dashed', size = 0.5) +
  geom_hline(data = pop_par, mapping = aes(yintercept = value + se_sa), 
             lty = 'dashed', size = 0.5) +
  facet_wrap( ~ parameter, ncol = 4, scales = "free") +
  theme(legend.position = "none") +
  xlab('Iteración')  + ylab('Valor')

G2

ggsave(G2, filename = 'G2.pdf', device = 'pdf', width = 8, height = 6)
ggsave(G2, filename = 'G2.png', device = 'png', width = 8, height = 6, dpi = 300)

##########################################################################-
# Identificación de set de valores iniciales atípico ----------------------
##########################################################################-
# Código de diseño factorial en escala cualitativa nominal
Cod_df <- LL_df %>% 
  filter(criteria == "-2LL") %>% 
  select(Iteration, Cl, V1, Q, V2)

# Búsqueda de iteración atípica, formación
Param_df %>% 
  select(Iteration, parameter, value) %>% 
  inner_join(Cod_df, by = "Iteration") %>% 
  mutate(Iteration = as.double(Iteration)) %>% 
  spread(key = parameter, value = value) %>% 
  knitr::kable(.)
  

# El resultado "atípico" proviene del caso N.33 que consiste en valores de 
# Cl - S, V1 - L, Q - M, V2 - L.


Param_df %>%
  filter(parameter == 'omega_V1') %>% 
  select(Iteration, parameter, value, se_sa, rse_sa) %>% 
  knitr::kable(.)

# |Iteration |parameter |     value|     se_sa|    rse_sa|
# |:---------|:---------|---------:|---------:|---------:|
# |2         |omega_V1  | 0.2327993| 0.3646411| 156.63327|
# |42        |omega_V1  | 0.2524595| 0.5331203| 211.17061|
# |52        |omega_V1  | 0.2415468| 1.2453043| 515.55395|
# |67        |omega_V1  | 0.2203358| 0.4544909| 206.27192|

# | Iteration|Cl |V1 |Q  |V2 |   
# |---------:|:--|:--|:--|:--|
# |         2|L  |M  |M  |M  |
# |        42|S  |L  |L  |L  |
# |        52|M  |S  |S  |L  |
# |        67|M  |L  |L  |S  |




