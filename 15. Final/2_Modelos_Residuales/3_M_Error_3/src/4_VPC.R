##------------------------------------------------------------------------------#
## Nombre del Script: Script de generación de VPC -------------------------------
##  
## Propósito del Script: crear gráficos de chequeo predictivo visual mediante  
##  simulación por medio del paquete MLXR de Monolix. 
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 18-02-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------#
# Introducción -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Carga de paquetes
require(tidyverse)
require(rlang)
require(mlxR)

#-------------------------------------------------------------------------------#
auxdir <- file.path('MODELO_FINAL', 'BASE_MODEL')

#-------------------------------------------------------------------------------#
# Apertura y modificación de archivo de datos con observaciones
#-------------------------------------------------------------------------------#
#................................................................................
#  1 Abrir archivo de observaciones TAD
#  2 Filtrar sólo eventos de observación (EVID == 0)
#  3 Renombrar columnas para armonizar con archivo de simulación
#  4 Cambiar el tipo de columna para ID a factor
#  5 Reasignar la variable data_TAD como una nueva variable tipo lista 
#  para el comando SimulX.
#................................................................................

# Lectura de archivo de observaciones
data_TAD <- read_delim(
  file.path('MODELO_FINAL', 'DATA', '2_data_TAD.csv'),
  ",", escape_double = FALSE, locale = locale(), 
  trim_ws = TRUE, na = ".")

data_OBS <- data_TAD %>% 
  filter(EVID == 0) %>% 
  rename(time = TIME, id = ID, y_1 = DV) %>% 
  mutate(id = factor(id))

data_TAD <- 
  file.path('MODELO_FINAL', 'DATA', '2_data_TAD.csv') %>% 
  mlxR::readDatamlx(datafile = ., 
                    header = c("id", "y", "mdv", "evid",	"amt", "rate", 
                               "addl", "ii", "ss", "time", 
                               rep('ignore', 5), 'contcov', rep('ignore', 15)))

#-------------------------------------------------------------------------------#
# Simulación de parámetros teóricos ---------------------------------------
#-------------------------------------------------------------------------------#
# El objetivo es simular 1000 set de datos con el diseño de dosis y 
# covariables original. Con esto se generan varios modelos teóricas que se 
# deben comparar con las observaciones empíricas.
#................................................................................
#  1 Seleccionar el archivo con el Modelo de Monolix obtenido, este fue 
#  modificado con un archivo de datos con TAD en lugar de TSFD. 
#  2 Abrir el archivo de parámetros poblacionales
#  3 Asignar los nombres a un vector
#  4 Ajustar las observaciones para tener 1000 puntos entre el intervalo 
#  de dosificación de 0 a 8 horas, almacenar esta configuración en una lista.
#  5 Configurar un vector con parámetros a simular. 
#  6 Pre-localizar un vector con 1000 posiciones.
#  7 Colocar en cada posición del vector un data.frame con simulaciones 
#  originadas a partir del diseño presente en el set de datos original.
#  8 Realizar la simulación teniendo en cuenta el archivo de proyecto con 
#  *simulx*. Elaborar una lista con 1000 tablas que contienen las simulaciones. 
# Se demora 6 minutos aproxim.
#................................................................................
project.file <- 'MODELO_FINAL/BASE_MODEL.mlxtran'

# Parámetros poblacionales
param <-
  read_csv(file.path(auxdir, "populationParameters.txt"), 
           col_types = cols()) %>%
  select(-se_sa, -rse_sa) %>%
  column_to_rownames(var = "parameter")

param <- setNames(pull(param), as.character(rownames(param)))
out1  <- list(name = 'y_1', time = seq(0, 8, length.out = 1e3))
data_list <- vector(mode = "list", length = 1000)

# ptm <- proc.time()
for (i in 1:1000) {
  data_list[i] <- simulx(
    project = project.file,
    parameter = param,
    treatment = data_TAD$treatment,
    output = list(out1), 
    settings = list(load.design=FALSE)
  )['y_1']
  print(paste('Lista la interacción N.º: ', i))
}
# print(proc.time() - ptm)
# Demora 6 minutos

#-------------------------------------------------------------------------------#
# Modificación de los archivos de datos -----------------------------------
#-------------------------------------------------------------------------------#
# Calcular intervalos de predicción del 80% (P10, P50, P100) para cada una 
# de las simulaciones teóricas (dentro de cada posición de la lista). 
#................................................................................
#  1 Seleccionar lista con data frames *data_list*
#  2 En cada posición, agregar una columna que indique el grupo al que 
#  pertenece, con esto se generan 100 grupos por ordenación de la variable 
#  tiempo. 
#  3 En cada posición, agrupar el data frame por la variable grupo. 
#  4 En cada posición, resumir la variable tiempo por su valor promedio, en 
#  cada grupo calcular la media e IP80% para la variable concentración (IP). 
#  5 Convertir la lista en un data.frame, este contiene una columna que 
#  referencia a la posición.
#................................................................................

options(dplyr.summarise.inform = FALSE) # Elimina unos mensajes molestos

dfr_percs <- data_list %>%
  map(~ mutate(.x, gr = ntile(time, 100))) %>% 
  map(~ group_by(.x, gr)) %>% 
  map( ~ summarise(.x, TIME = mean(time), 
                   ME = quantile(x = y_1, probs = 0.50), 
                   LI = quantile(x = y_1, probs = 0.10),
                   LS = quantile(x = y_1, probs = 0.90))) %>% 
  map_dfr(~ as.data.frame(.x), .id = 'sim_id') 

#-------------------------------------------------------------------------------#
# Calcular intervalos de confianza del 95% para cada intervalo de predicción
# Se debe realizar un resumen del valor de IP dentro de cada bin, entre todas 
# las simulaciones, se obtienen medianas y percentiles P5% y P95% 
#................................................................................
#  1 Tomar el data frame con los intervalos de predicción calculados
#  2 Agrupar por bins (hay 100 x 1000)
#  3 Resumir en cada bin el valor promedio del bin e IC95% para el límite 
#  inferior y el límite superior del IP
#................................................................................

dfr_percs1 <- dfr_percs %>% 
  group_by(gr) %>% 
  summarise(TIME = mean(TIME),
            ME_me = quantile(ME, probs = 0.50),
            ME_li = quantile(ME, probs = 0.05),
            ME_ls = quantile(ME, probs = 0.95),
            
            LI_me = quantile(LI, probs = 0.50),
            LI_li = quantile(LI, probs = 0.05),
            LI_ls = quantile(LI, probs = 0.95),
            
            LS_me = quantile(LS, probs = 0.50),
            LS_li = quantile(LS, probs = 0.05),
            LS_ls = quantile(LS, probs = 0.95) )

#-------------------------------------------------------------------------------#
# Calcular percentiles empíricos de las observaciones
#................................................................................
#  1 Seleccionar el archivo de observaciones
#  2 Agregar una columna con 6 bins, de acuerdo a la ordenación de la variable 
#  time, estos grupos contienen 15 puntos (pues ahí 15 individuos).
#  3 En cada bin, obtener el promedio del tiempo, la mediana de DV, y 
#  percentiles P10-P90.
#................................................................................

data_OBS1 <- data_OBS %>% 
  mutate(gr = ntile(time, 6)) %>% 
  group_by(gr) %>% 
  summarise(TIME = mean(time), 
            ME = quantile(x = y_1, probs = 0.50), 
            LI = quantile(x = y_1, probs = 0.10),
            LS = quantile(x = y_1, probs = 0.90),
            n = n())

#-------------------------------------------------------------------------------#
# Creación del gráfico VPC con percentiles --------------------------------
#-------------------------------------------------------------------------------#
# Función que adiciona lineas y puntos, diseñada para observaciones
linedots <- function(data, x, y) {
  x = rlang::ensym(x)
  y = rlang::ensym(y)
  return(
    list(geom_line(data = data, aes(x = !!x, y = !!y)),
         geom_point(data = data, aes(x = !!x, y = !!y)))
  )
}

# VPC con percentiles predichos
theme_set(theme_bw() +
            theme(panel.border = element_rect(fill = NULL, colour = 'black')))

g_percs <- 
  dfr_percs1 %>%
  ggplot(aes(x = TIME)) +
  geom_ribbon(aes(ymin = ME_li, ymax = ME_ls), alpha = 0.5, fill = 'pink') +
  geom_ribbon(aes(ymin = LI_li, ymax = LI_ls), alpha = 0.5, fill = 'blue') +
  geom_ribbon(aes(ymin = LS_li, ymax = LS_ls), alpha = 0.5, fill = 'blue') +
  linedots(data_OBS1, TIME, ME) +
  linedots(data_OBS1, TIME, LI) +
  linedots(data_OBS1, TIME, LS) +
  coord_cartesian(ylim = c(0, 100)) +
  xlab('TAD') + ylab(expression(C[P]~(mg/L)))

# Almacenamiento en formato PDF
ggsave('figures/VPC_percentil.pdf', g_percs, 'pdf', 
       width = 6, height = 5)

#-------------------------------------------------------------------------------#
# Simulación con pcVPC ----------------------------------------------------
#-------------------------------------------------------------------------------#
# El objetivo es simular un set de datos con el diseño de dosis y 
# covariables original. Este set de datos se comparte en todas las simulaciones 
# ya realizadas por que son predicciones poblacionales. Los parámetros de 
# variabilidad se ajustan a cero para reflejar el valor esperado de acuerdo 
# al diseño y el tiempo seleccionado.
#................................................................................
#  1 Ajustar los parámetros de variabilidad a cero.
#  2 Realizar la simulación teniendo en cuenta el archivo de proyecto con 
#  *simulx*. Este consiste de una sola tabla que se compartirá en todas 
#  las tablas de data_list, y se trata de simulaciones poblacionales (PRED).
#................................................................................

param["omega_Cl"] <- 0
param["omega_V1"] <- 0
param["omega_Q"] <- 0
param["omega_V2"] <- 0
param["a"] <- 0
param["b"] <- 0

data_pred <- simulx(project = project.file,
                    parameter = param, 
                    treatment = data_TAD$treatment,
                    output = out1)[['y_1']]

#-------------------------------------------------------------------------------#
# Cálculo del valor medio de PRED -----------------------------------------
#-------------------------------------------------------------------------------#
# El cálculo del valor medio de PRED es importante para hacer la correción 
# en los valores de simulación.
#................................................................................
#  1 Tomar *data_pred*
#  2 Con el df data_pred se deben generar 100 grupos (bins) de acuerdo al 
#  ordenamiento de tiempo.
#  3 Agrupar el df de acuerdo al grupo.
#  4 Resumir cada bin por el valor medio de tiempo y la mediana de las 
#  concentraciones (en este caso predicciones poblacionales).
#................................................................................

data_pred_1 <- data_pred %>%
  mutate(gr = ntile(time, 100)) %>%
  group_by(gr) %>% 
  summarise(time_bin = mean(time),
            ME_bin = quantile(x = y_1, probs = 0.50))

#-------------------------------------------------------------------------------#
# Adicionar los valores medios de PRED en el bin al data.frame con las 
# predicciones poblacionales de acuerdo al diseño.
#................................................................................
#  1 Seleccionar data_pred
#  2 Calcular los bins, teniendo en cuenta que deben haber 100 de ellos.
#  3 Adicionar los valores medios de PRED a los data_pred en cada grupo 
#  correspondiente.
#  4 Redondear los tiempos en el df a 4 cifras decimales
#  5 Convertir la tabla en tibble.   
#................................................................................

data_pred = data_pred %>%
  mutate(gr = ntile(time, 100)) %>%
  left_join(., data_pred_1, by = 'gr') %>% 
  mutate(time = round(time, 4)) %>% 
  as_tibble(.) 

#-------------------------------------------------------------------------------#
# Adicionar los valores de PRED y PRED_bin_med
#................................................................................
#  1 Seleccionar el archivo de lista con data.frames de cada simulación teórica. 
#  2 En cada df, se debe ajustar el tiempo a 4 cifras significativas.
#  3 En cada df, se debe crear una variable con 100 bins.
#  4 En cada df, se deben adicionar las columnas relacionadas con PRED, 
#  PRED_bin_med
#................................................................................

data_list <-  data_list %>%
  map( ~ mutate(.x, time = round(time, 4))) %>% 
  map(~ mutate(.x, gr = ntile(time, 100))) %>%
  map(~ left_join(.x, data_pred, 
                  by = c('id', 'gr', 'time')))

#-------------------------------------------------------------------------------#
# Modificación de los archivos de datos -----------------------------------
#-------------------------------------------------------------------------------#
# Calcular intervalos de predicción del 80% (P10, P50, P100) para cada una 
# de las simulaciones teóricas (dentro de cada posición de la lista). 
#................................................................................
#  1 Seleccionar la lista con los data frames *data_list*
#  2 Convertir la lista en un data.frame, este contiene una columna que 
#  referencia a la posición *sim_id*.
#  3 Calcular una columna con simulaciones corregidas por predicción PPRED
#  4 Agrupar el data.frame por la variable sim_id y grupo. 
#  5 Resumir la variable tiempo por su valor promedio, en cada grupo calcular 
#  la media e IP80% para la variable concentración (IP). 
#................................................................................

dfr_percs_pcVPC <- data_list %>% 
  map_dfr( ~ as.data.frame(.x), .id = 'sim_id') %>% 
  mutate(pcy_1 = y_1.x * ME_bin / y_1.y) %>% 
  group_by(sim_id, gr) %>% 
  summarise(TIME = mean(time),
            ME = quantile(x = pcy_1, probs = 0.50),
            LI = quantile(x = pcy_1, probs = 0.10),
            LS = quantile(x = pcy_1, probs = 0.90))

#-------------------------------------------------------------------------------#
# Calcular intervalos de confianza del 95% para cada intervalo de predicción
# Se debe realizar un resumen del valor de IP dentro de cada bin, entre todas 
# las simulaciones, se obtienen medianas y percentiles P5% y P95% 
#................................................................................
#  1 Tomar el data frame con los intervalos de predicción calculados
#  2 Agrupar por bins (hay 100 x 1000)
#  3 Resumir en cada bin el valor promedio del bin e IC95% para el límite 
#  inferior y el límite superior del IP
#................................................................................

dfr_percs1_pcVPC <- dfr_percs_pcVPC %>% 
  ungroup(.) %>% 
  group_by(gr) %>% 
  summarise(TIME = mean(TIME),
            ME_me = quantile(ME, probs = 0.50),
            ME_li = quantile(ME, probs = 0.05),
            ME_ls = quantile(ME, probs = 0.95),
            
            LI_me = quantile(LI, probs = 0.50),
            LI_li = quantile(LI, probs = 0.05),
            LI_ls = quantile(LI, probs = 0.95),
            
            LS_me = quantile(LS, probs = 0.50),
            LS_li = quantile(LS, probs = 0.05),
            LS_ls = quantile(LS, probs = 0.95) )

#-------------------------------------------------------------------------------#
# Cálculo de valor de PRED para las observaciones
#  Se simula el valor de PRED teniendo en cuenta el modelo sin variabilidad, +
#   y con el mismo archivo de proyecto. En este caso, se tienen en cuenta 
#   los muestreos de observaciones iguales al set de datos originales. 
#   Se renombra la columna y_1 como PRED.

out2 <- list(name = 'y_1', time = data_TAD$y$time)

data_OBS_PRED <- simulx(project = project.file,
                        parameter = param, 
                        treatment = data_TAD$treatment,
                        output = out2)[['y_1']] %>% 
  rename(PRED = y_1)

#-------------------------------------------------------------------------------#
# Con el archivo de datos PRED obtenidos para el modelo original se calcula 
# PRED_BIN_MED.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar el archivo de datos originales
##  2 Adicionar una columna con el bin correspondiente de 6 posibles valores 
##  (para las observaciones).
##  3 Agrupar el df por la variable gr que corresponde al bin.
##  4 Resumir en cada bin: time_bin que es el promedio del tiempo, y ME_PRED 
##  que el valor medio dela predicción poblacional de y_1.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data_OBS_PRED_tile <- data_OBS_PRED %>%
  mutate(gr = ntile(time, 6)) %>%
  group_by(gr) %>%
  summarise(time_bin = mean(time),
            ME_PRED = quantile(x = PRED, probs = 0.50))

#-------------------------------------------------------------------------------#
# Adicionar los datos de ME_PRED al archivo que contiene los PRED
#................................................................................
#  1 Seleccionar el archivo con ME_PRED
#  2 Agregar una columna con 6 bins, de acuerdo a la ordenación de la variable 
#  time, estos grupos contienen 15 puntos (pues ahí 15 individuos).
#  3 Adicionar el archivo con ME_PRED al archivo con PRED.
#................................................................................

data_OBS_PRED <- data_OBS_PRED %>%
  mutate(gr = ntile(time, 6)) %>%
  left_join(., data_OBS_PRED_tile, by = 'gr') %>% 
  mutate(time = round(time, 4)) 

#-------------------------------------------------------------------------------#
# Adicionar los datos de PRED y ME_PRED al archivo de observaciones originales
#................................................................................
#  1 Seleccionar data_OBS
#  2 Unir el archivo de datos con PRED y ME_PRED al archivo con las 
#  observaciones originales.
#  3 Calcular las predicciones corregidas por PRED para las observaciones
#  4 Agrupar por gr (bins)
#  5 Resumir por tiempo, mediana y percentiles empíricos del 10% y 90% para 
#  las observaciones.
#................................................................................

data_OBS_PRED_sum <- data_OBS %>%
  mutate(time = round(time, 4)) %>% 
  left_join(., data_OBS_PRED, by = c('id', 'time')) %>%
  mutate(pcVPC = y_1 * ME_PRED / PRED) %>%
  group_by(gr) %>%
  summarise(
    TIME = mean(time),
    ME = quantile(x = y_1, probs = 0.50),
    LI = quantile(x = y_1, probs = 0.10),
    LS = quantile(x = y_1, probs = 0.90),
    n = n()
  )


#-------------------------------------------------------------------------------#
# VPC con percentiles predichos corregidos por PRED
g_percs1 <- 
  dfr_percs1_pcVPC %>%
  ggplot(aes(x = TIME)) +
  geom_ribbon(aes(ymin = ME_li, ymax = ME_ls), alpha = 1.0, fill = '#FFF569') +
  geom_ribbon(aes(ymin = LI_li, ymax = LI_ls), alpha = 1.0, fill = '#7595FF') +
  geom_ribbon(aes(ymin = LS_li, ymax = LS_ls), alpha = 1.0, fill = '#7595FF') +
  linedots(data_OBS_PRED_sum, TIME, ME) +
  linedots(data_OBS_PRED_sum, TIME, LI) +
  linedots(data_OBS_PRED_sum, TIME, LS) +
  coord_cartesian(ylim = c(0, 100)) +
  xlab('TAD') + ylab('Concentraciones corregidas por PRED')

g_percs1
# Almacenamiento del archivo PDF
ggsave(filename = 'figures/pcVPC_percentil.pdf', g_percs1, 
       device = 'pdf', width = 6, height = 5)

# Eliminación del objeto data_list que es muy grande como para ser almacenado 
# y transferido
rm(data_list)

#-------------------------------------------------------------------------------#
# Almacenar RDS
saveRDS(g_percs, 'figures/RDS/VPC_percentil.rds')
saveRDS(g_percs1, 'figures/RDS/pcVPC_percentil.rds')
