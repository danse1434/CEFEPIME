##------------------------------------------------------------------------#
## Nombre del Script: Script de generación de VPC  ------------------------
##  
## Proposito del Script: crear gráficos de chequeo predictivo visual mediante  
##  simulación por medio del paquete MLXR de Monolix. 
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion: 18-02-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
##########################################################################-
# Introducción -----------------------------------------------------
##########################################################################-
# Carga de paquetes
require(tidyverse)
require(rlang)
require(mlxR)
##########################################################################-
# Selección de directorio principal
setwd(file.path('C:', 'Users', 'Daniel', 'OneDrive', 'Documents', 
                '(Proyecto)_Estudio_PKPD', 'CEFEPIME', '10. Base_Refin', 
                'BASE', 'BASE_MODEL'))

# Apertura y modificación de archivo de datos con observaciones
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Abrir archivo de observaciones TAD
##  2 Filtrar sólo eventos de observación (EVID == 0)
##  3 Renombrar columnas para armonizar con archivo de simulación
##  4 Cambiar el tipo de columna para ID a factor
##  5 Reasignar la variable data_TAD como una nueva variable tipo lista 
##  para el comando SimulX.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Lectura de archivo de observaciones
data_TAD <- read_delim("../Monolix_data_TAD.csv", 
                       ",", escape_double = FALSE, locale = locale(), 
                       trim_ws = TRUE, na = ".")

data_OBS <- data_TAD %>% 
  filter(EVID == 0) %>% 
  rename(time = TIME, id = ID, y_1 = DV) %>% 
  mutate(id = factor(id))


data_TAD <- 
  mlxR::readDatamlx(datafile = '../Monolix_data_TAD.csv', 
                    header = c("id", "y", "mdv", "evid",	"amt", "rate", 
                               "addl", "ii", "ss", "time", rep('ignore', 22)))

##########################################################################-
# Simulación de parámetros teóricos ---------------------------------------
##########################################################################-
# El objetivo es simular 1000 set de datos con el diseño de dosis y 
# covariables original. Con esto se generan varios modelos teóricas que se 
# deben comparar con las observaciones empíricas.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar el archivo con el Modelo de Monolix obtenido, este fue 
##  modificado con un archivo de datos con TAD en lugar de TSFD. 
##  2 Abrir el archivo de parámetros poblacionales
##  3 Asignar los nombres a un vector
##  4 Ajustar las observaciones para tener 1000 puntos entre el intervalo 
##  de dosificación de 0 a 8 horas, almacenar esta configuración en una lista.
##  5 Configurar un vector con parámetros a simular. 
##  6 Pre-localizar un vector con 1000 posiciones.
##  7 Colocar en cada posición del vector un data.frame con simulaciones 
##  originadas a partir del diseño presente en el set de datos original.
##  8 Realizar la simulación teniendo en cuenta el archivo de proyecto con 
##  *simulx*. Elaborar una lista con 1000 tablas que contienen las simulaciones. 
# Se demora 6 minutos aproxim
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
project.file <- '../BASE_MODEL.mlxtran'

# Parámetros poblacionales
param <-
  read_csv("populationParameters.txt", col_types = cols()) %>%
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

##########################################################################-
# Modificación de los archivos de datos -----------------------------------
##########################################################################-
# Calcular intervalos de predicción del 80% (P10, P50, P100) para cada una 
# de las simulaciones teóricas (dentro de cada posición de la lista). 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar lista con data frames *data_list*
##  2 En cada posición, agregar una columna que indique el grupo al que 
##  pertenece, con esto se generan 100 grupos por ordenación de la variable 
##  tiempo. 
##  3 En cada posición, agrupar el data frame por la variable grupo. 
##  4 En cada posición, resumir la variable tiempo por su valor promedio, en 
##  cada grupo calcular la media e IP80% para la variable concentración (IP). 
##  5 Convertir la lista en un data.frame, este contiene una columna que 
##  referencia a la posición.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

dfr_percs <- data_list %>%
  map(~ mutate(.x, gr = ntile(time, 100))) %>% 
  map(~ group_by(.x, gr)) %>% 
  map( ~ summarise(.x, TIME = mean(time), 
                   ME = quantile(x = y_1, probs = 0.50), 
                   LI = quantile(x = y_1, probs = 0.10),
                   LS = quantile(x = y_1, probs = 0.90))) %>% 
  map_dfr(~ as.data.frame(.x), .id = 'sim_id') 
  
##########################################################################-
# Calcular intervalos de confianza del 95% para cada intervalo de predicción
# Se debe realizar un resumen del valor de IP dentro de cada bin, entre todas 
# las simulaciones, se obtienen medianas y percentiles P5% y P95% 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Tomar el data frame con los intervalos de predicción calculados
##  2 Agrupar por bins (hay 100 x 1000)
##  3 Resumir en cada bin el valor promedio del bin e IC95% para el límite 
##  inferior y el límite superior del IP
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

##########################################################################-
# Calcular percentiles empíricos de las observaciones
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar el archivo de observaciones
##  2 Agregar una columna con 6 bins, de acuerdo a la ordenación de la variable 
##  time, estos grupos contienen 15 puntos (pues ahí 15 individuos).
##  3 En cada bin, obtener el promedio del tiempo, la mediana de DV, y 
##  percentiles P10-P90.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data_OBS1 <- data_OBS %>% 
  mutate(gr = ntile(time, 6)) %>% 
  group_by(gr) %>% 
  summarise(TIME = mean(time), 
            ME = quantile(x = y_1, probs = 0.50), 
            LI = quantile(x = y_1, probs = 0.10),
            LS = quantile(x = y_1, probs = 0.90),
            n = n())

##########################################################################-
# Creación del gráfico VPC con percentiles --------------------------------
##########################################################################-
# Función que adiciona lineas y puntos, diseñada para observaciones
linedots <- function(data, x, y) {
  x = rlang::ensym(x)
  y = rlang::ensym(y)
  return(
    list(geom_line(data = data, aes(x = !!x, y = !!y), col='red'),
         geom_point(data = data, aes(x = !!x, y = !!y),col='red',shape=15))
  )
}

# VPC con percentiles predichos
g_percs <- dfr_percs1 %>%
  ggplot(aes(x = TIME)) +
  geom_ribbon(aes(ymin = ME_li, ymax = ME_ls), alpha = 0.5, fill = 'gray70') +
  geom_ribbon(aes(ymin = LI_li, ymax = LI_ls), alpha = 0.5, fill = 'gray30') +
  geom_ribbon(aes(ymin = LS_li, ymax = LS_ls), alpha = 0.5, fill = 'gray30') +
  geom_line(aes(y=ME_me), lty='dashed') +
  geom_line(aes(y=LI_me), lty='dashed') +
  geom_line(aes(y=LS_me), lty='dashed') +
  linedots(data_OBS1, TIME, ME) +
  linedots(data_OBS1, TIME, LI) +
  linedots(data_OBS1, TIME, LS) +
  theme_bw() +
  coord_cartesian(ylim = c(0, 100)) +
  theme(panel.border = element_rect(fill = NULL, colour = 'black')) +
  xlab('TAD, tiempo tras dosis (h)') + 
  ylab('Concentración plasmática FEP (mg/L)')

# Almacenamiento en formato PDF
ggsave(filename = './FIGURAS/VPC_percentil.pdf', g_percs, device = 'pdf', 
       width = 6, height = 4)















##########################################################################-
# Simulación con pcVPC ----------------------------------------------------
##########################################################################-
# El objetivo es simular un set de datos con el diseño de dosis y 
# covariables original. Este set de datos se comparte en todas las simulaciones 
# ya realizadas por que son predicciones poblacionales. Los parámetros de 
# variabilidad se ajustan a cero para reflejar el valor esperado de acuerdo 
# al diseño y el tiempo seleccionado.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Ajustar los parámetros de variabilidad a cero.
##  2 Realizar la simulación teniendo en cuenta el archivo de proyecto con 
##  *simulx*. Este consiste de una sola tabla que se compartirá en todas 
##  las tablas de data_list, y se trata de simulaciones poblacionales (PRED).
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

##########################################################################-
# Cálculo del valor medio de PRED -----------------------------------------
##########################################################################-
# El cálculo del valor medio de PRED es importante para hacer la correción 
# en los valores de simulación.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Tomar *data_pred*
##  2 Con el df data_pred se deben generar 100 grupos (bins) de acuerdo al 
##  ordenamiento de tiempo.
##  3 Agrupar el df de acuerdo al grupo.
##  4 Resumir cada bin por el valor medio de tiempo y la mediana de las 
##  concentraciones (en este caso predicciones poblacionales).
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data_pred_1 <- data_pred %>%
  mutate(gr = ntile(time, 100)) %>%
  group_by(gr) %>% 
  summarise(time_bin = mean(time),
            ME_bin = quantile(x = y_1, probs = 0.50))

##########################################################################-
# Adicionar los valores medios de PRED en el bin al data.frame con las 
# predicciones poblacionales de acuerdo al diseño.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar data_pred
##  2 Calcular los bins, teniendo en cuenta que deben haber 100 de ellos.
##  3 Adicionar los valores medios de PRED a los data_pred en cada grupo 
##  correspondiente.
##  4 Redondear los tiempos en el df a 4 cifras decimales
##  5 Convertir la tabla en tibble.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data_pred = data_pred %>%
  mutate(gr = ntile(time, 100)) %>%
  left_join(., data_pred_1, by = 'gr') %>% 
  mutate(time = round(time, 4)) %>% 
  as_tibble(.) 
  
##########################################################################-
# Adicionar los valores de PRED y PRED_bin_med
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar el archivo de lista con data.frames de cada simulación teórica. 
##  2 En cada df, se debe ajustar el tiempo a 4 cifras significativas.
##  3 En cada df, se debe crear una variable con 100 bins.
##  4 En cada df, se deben adicionar las columnas relacionadas con PRED, 
##  PRED_bin_med
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
data_list <-  data_list %>%
    map( ~ mutate(.x, time = round(time, 4))) %>% 
    map(~ mutate(.x, gr = ntile(time, 100))) %>%
    map(~ left_join(.x, data_pred, 
                    by = c('id', 'gr', 'time')))

##########################################################################-
# Modificación de los archivos de datos -----------------------------------
##########################################################################-
# Calcular intervalos de predicción del 80% (P10, P50, P100) para cada una 
# de las simulaciones teóricas (dentro de cada posición de la lista). 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar la lista con los data frames *data_list*
##  2 Convertir la lista en un data.frame, este contiene una columna que 
##  referencia a la posición *sim_id*.
##  3 Calcular una columna con simulaciones corregidas por predicción PPRED
##  4 Agrupar el data.frame por la variable sim_id y grupo. 
##  5 Resumir la variable tiempo por su valor promedio, en cada grupo calcular 
##  la media e IP80% para la variable concentración (IP). 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

dfr_percs_pcVPC <- data_list %>% 
  map_dfr( ~ as.data.frame(.x), .id = 'sim_id') %>% 
  mutate(pcy_1 = y_1.x * ME_bin / y_1.y) %>% 
  group_by(sim_id, gr) %>% 
  summarise(TIME = mean(time),
            ME = quantile(x = pcy_1, probs = 0.50),
            LI = quantile(x = pcy_1, probs = 0.10),
            LS = quantile(x = pcy_1, probs = 0.90))

##########################################################################-
# Calcular intervalos de confianza del 95% para cada intervalo de predicción
# Se debe realizar un resumen del valor de IP dentro de cada bin, entre todas 
# las simulaciones, se obtienen medianas y percentiles P5% y P95% 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Tomar el data frame con los intervalos de predicción calculados
##  2 Agrupar por bins (hay 100 x 1000)
##  3 Resumir en cada bin el valor promedio del bin e IC95% para el límite 
##  inferior y el límite superior del IP
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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

##########################################################################-
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

##########################################################################-
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

##########################################################################-
# Adicionar los datos de ME_PRED al archivo que contiene los PRED
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar el archivo con ME_PRED
##  2 Agregar una columna con 6 bins, de acuerdo a la ordenación de la variable 
##  time, estos grupos contienen 15 puntos (pues ahí 15 individuos).
##  3 Adicionar el archivo con ME_PRED al archivo con PRED.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data_OBS_PRED <- data_OBS_PRED %>%
  mutate(gr = ntile(time, 6)) %>%
  left_join(., data_OBS_PRED_tile, by = 'gr') %>% 
  mutate(time = round(time, 4)) 

##########################################################################-
# Adicionar los datos de PRED y ME_PRED al archivo de observaciones originales
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar data_OBS
##  2 Unir el archivo de datos con PRED y ME_PRED al archivo con las 
##  observaciones originales.
##  3 Calcular las predicciones corregidas por PRED para las observaciones
##  4 Agrupar por gr (bins)
##  5 Resumir por tiempo, mediana y percentiles empíricos del 10% y 90% para 
##  las observaciones.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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


# VPC con percentiles predichos corregidos por PRED
g_percs1 <-
  dfr_percs1_pcVPC %>%
  ggplot(aes(x = TIME)) +
  geom_ribbon(aes(ymin = ME_li, ymax = ME_ls), alpha = 0.5, fill = 'gray70') +
  geom_ribbon(aes(ymin = LI_li, ymax = LI_ls), alpha = 0.5, fill = 'gray30') +
  geom_ribbon(aes(ymin = LS_li, ymax = LS_ls), alpha = 0.5, fill = 'gray30') +
  geom_line(aes(y=ME_me), lty='dashed') +
  geom_line(aes(y=LI_me), lty='dashed') +
  geom_line(aes(y=LS_me), lty='dashed') +
  linedots(data_OBS_PRED_sum, TIME, ME) +
  linedots(data_OBS_PRED_sum, TIME, LI) +
  linedots(data_OBS_PRED_sum, TIME, LS) +
  theme_bw() +
  theme(panel.border = element_rect(fill = NULL, colour = 'black')) +
  coord_cartesian(ylim = c(0, 100)) +
  xlab('TAD, tiempo tras dosis (h)') + 
  ylab('Concentración plasmática FEP \n Corregida por concentración (mg/L)')

g_percs1
# Almacenamiento del archivo PDF
ggsave(filename = './FIGURAS/pcVPC_percentil.pdf', g_percs1, 
       device = 'pdf', width = 6, height = 5)

# Eliminación del objeto data_list que es muy grande como para ser almacenado 
# y transferido
rm(data_list)


















##########################################################################-
# Chequeo Predictivo Numérico (NPC) ---------------------------------------
##########################################################################-
# Parámetros poblacionales
param <-
  read_csv("populationParameters.txt", col_types = cols()) %>%
  select(-se_sa, -rse_sa) %>%
  column_to_rownames(var = "parameter")

param <- setNames(pull(param), as.character(rownames(param)))

# Prealocación de lista
npc_list <- vector(mode = "list", length = 500)


list_individual <- function(i) {
  select_tto <- function(var, i) {
    data_TAD$treatment[data_TAD$treatment[, 'id'] == i, var]
  }
  
  inter1 = list(
    time = select_tto('time', i),
    amount = select_tto('amount', i),
    rate = select_tto('rate', i)
  )
  
  select_obs <- function(var, i) {
    data_TAD$y[data_TAD$y[, 'id'] == i, var]
  }
  
  inter2 = list(name = 'y_1',
                time = select_obs('time', i),
                y_1 = select_obs('y', i))
  
  return(
    list(treatment = inter1, output = inter2, size=1, level='longitudinal')
    )
}
list_indivi_tot <- list()

for (i in 1:15) {
  list_indivi_tot[[i]] <- list_individual(i)
}

model1 <- inlineModel("
[INDIVIDUAL]
input = {Cl_pop, Q_pop, omega_Q, V1_pop, omega_V1, V2_pop, omega_V2, omega_Cl}

DEFINITION:
Cl = {distribution=logNormal, typical=Cl_pop, sd=omega_Cl}
Q = {distribution=logNormal, typical=Q_pop, sd=omega_Q}
V1 = {distribution=logNormal, typical=V1_pop, sd=omega_V1}
V2 = {distribution=logNormal, typical=V2_pop, sd=omega_V2}

[LONGITUDINAL]
input = {a, b}

file = 'lib:infusion_2cpt_ClV1QV2.txt'

DEFINITION:
y_1 = {distribution=normal, prediction=Cc, errorModel=combined2(a, b)}
")

for (i in 1:500) {
  npc_list[i] <- simulx(model = model1,
                        parameter = param,
                        group = list_indivi_tot
                        )['y_1']
  print(paste('Lista la interacción N.º: ', i))
}


sim_vec <- npc_list %>%
  map_dfr( ~ as.data.frame(.x), .id = 'sim_id') %>%
  as_tibble() %>%
  magrittr::use_series("y_1")


obs_vec <- data_TAD$y$y


npc <- function (obs, sim,PI){
  nobs <- length(obs)
  nrep <- length(sim)[1]/nobs
  percentiles <- c(50-PI/2,(50+PI/2))
  
  matsimfull <- matrix(sim, ncol = nrep)
  matsimper <- t(apply(matsimfull,1,quantile,percentiles/100))
  matsimper_lower <- matsimper[,1:length(PI)]
  matsimper_upper <- matsimper[,(length(PI)+1):length(percentiles)]

  sim_outlower <- c()
  sim_outupper <- c()
  
  for (i in 1:length(PI)) {
    matsim_outlower <- ifelse(matsimfull < matsimper_lower[, i], 1, 0)
    
    sim_outlower <-
      cbind(sim_outlower, quantile(apply(matsim_outlower, 2, mean), c(0.025, 0.5, 0.975)))
    
    matsim_outupper <-
      ifelse(matsimfull > matsimper_upper[, i], 1, 0)
    
    sim_outupper <-
      cbind(sim_outupper, quantile(apply(matsim_outupper, 2, mean), c(0.025, 0.5, 0.975)))
  }
  
  yobs_lower <-
    yobs_upper <- matrix(rep(obs, length(PI)), ncol = length(PI))
  yobs_lower <- ifelse(yobs_lower < matsimper_lower, 1, 0)
  yobs_upper <- ifelse(yobs_upper > matsimper_upper, 1, 0)
  obsper_outlower <- apply(yobs_lower, 2, mean)
  obsper_outupper <- apply(yobs_upper, 2, mean)
  lower_cover <- obsper_outlower / sim_outlower[2, ]
  CI_lower_cover <- t(t(sim_outlower[c(1, 3), ]) / sim_outlower[2, ])
  upper_cover <- obsper_outupper / sim_outupper[2, ]
  CI_upper_cover <- t(t(sim_outupper[c(1, 3), ]) / sim_outupper[2, ])
  
  lower_out <-
    (lower_cover < CI_lower_cover[1, ]) + (lower_cover > CI_lower_cover[2, ])
  
  upper_out <-
    (upper_cover < CI_upper_cover[1, ]) + (upper_cover > CI_upper_cover[2, ])
  
  cover <- as.data.frame(cbind(
      rep(PI, 2),
      c(lower_cover, upper_cover),
      c(lower_out, upper_out),
      rbind(t(CI_lower_cover), t(CI_upper_cover))
    ), row.names = F)
  
  cover$Type <-
    rep(c("Limite IP inferior", "Limite IP superior"), each = length(PI))
  colnames(cover) <- c("PI","Ratio","Outliers","Lwr","Sup","Tipo")
  return(cover)}


df <-npc(obs_vec, sim_vec, PI = seq(0, 90, by = 5))  


NPC_plot <- df %>% 
  ggplot(aes(x = PI, y = Ratio, 
             group = Tipo, col = Tipo, fill = Tipo)) + 
  geom_ribbon(aes(ymin = Lwr, ymax = Sup), alpha = 0.8) +
  geom_line(col = 'black') + 
  geom_point(col = 'black') + 
  facet_grid(Tipo ~ .) + 
  geom_hline(yintercept = 1, lty = 'dotted') +
  ylab('Ratio O/E') + xlab('IP') + 
  theme_bw() +
  theme(panel.border = element_rect(fill = NULL, colour = 'black'),
        panel.grid = element_line(colour = 'NA'),
        legend.position = "none") +
  scale_fill_manual(values = c('#85B8FF', '#FF886B')) +
  scale_color_manual(values = c('#85B8FF', '#FF886B'))

NPC_plot
  

# Almacenamiento del archivo PDF
ggsave(filename = './FIGURAS/NPC_plot.pdf', NPC_plot, device = 'pdf', 
       width = 4/1.2, height = 6/1.2)














