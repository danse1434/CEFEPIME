## --------------------------- -
## Nombre del Script: Transformación de Datos de formato NONMEM a formato Monolix
##
## Propósito del Script: cambia un archivo plano con información farmacocinética 
## en formato NONMEM a un archivo plano en formato Monolix
##
## Autor: Daniel S. Parra González
##
## Fecha de creación: 2020-01-24
##
## Copyright (c) Daniel S. Parra, 2020
## Email: dsparrag@unal.edu.co
## 
## Carga de paquetes ------------------------------------------------------
require('tidyverse')
require('lubridate')
require('rlang')
require('grid')
require('data.table')
require('lixoftConnectors')

##########################################################################-
# Directorio Principal ----------------------------------------------------
##########################################################################-
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '03. Censura', 'MONOLIX'))

##########################################################################-
# Carga y alistamiento de archivo de datos --------------------------------
##########################################################################-
data <- read_csv(file = 'DATA/NONMEM_Cefepime.csv', 
                 col_names = T, 
                 na = c('.'),
                 col_types = cols(EVID = col_factor(),
                                  MDV = col_factor(),
                                  SS = col_factor(),
                                  SEXF = col_factor(), 
                                  ANTU = col_factor(),
                                  LLP = col_factor(),
                                  LMP = col_factor()))

##########################################################################-
# Modificación del archivo de acuerdo a requerimientos establecidos previamente
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## 1 Especificación de formato fecha en columna correspondiente
## 2 Aplicación fórmula Dubois & Dubois
## 3 Creación de columna con fecha y hora acoplada
## 4 Creación de columna con tipo de diagnóstico agregado
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
data = data %>% 
  mutate(DATE = dmy(DATE)) %>% 
  mutate(SCM2 = 0.007184*(WTKG^0.425)*(HCM^0.725)) %>% 
  mutate(DATEHOUR = dmy_hms(paste(data$DATE, data$TIME))) %>% 
  mutate(Dx = case_when(LLP == 1 ~ 'LLA',
                        LMP == 1 ~ 'LMA/LMC',
                        (LLP==0 & LMP==0) ~ 'Otros'))

##########################################################################-
# Selección de tipo de tiempo que se va a modelar
## Almacenar como expresión
t = "TSFD" %>%
  parse_expr(.) 
## Número de individuos en el estudio
nobs = unique(data$ID) %>% length()

##########################################################################-
# Especificación de Límite inferior de cuantificación ---------------------
##########################################################################-
LLOQ = 50/2^5
## Reemplazar esta variable en lugares correspondientes
# data = data %>%
#   mutate(CENSORING = if_else(DV <= LLOQ, 1, 0)) %>% 
#   mutate(DV = if_else(DV <= LLOQ, LLOQ, DV))



##########################################################################-
# Modificación a formato Monolix ------------------------------------------
##########################################################################-
# Modificación a formato Monolix
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Eliminar columnas TIME y DATE
##  2 Renombrar varias columnas para cumplir con nombres de Monolix
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data1 <-  data %>%
  select(-TIME, -DATE, -CMT) %>%
  rename(TIME = !!t) %>% 
  mutate('YTYPE' = ifelse(EVID == 0, 1, NA_real_)) %>%
  mutate('ADM' = ifelse(EVID == 1, 1, NA_real_)) %>% 
  select(-TAD, -DATEHOUR)
  # rename(
  #   'EVENT ID' = EVID,
  #   'IGNORED OBSERVATION' = MDV,
  #   'OBSERVATION' = DV,
  #   'AMOUNT' = AMT,
  #   'INFUSION RATE' = RATE,
  #   'ADDITIONAL DOSES' = ADDL,
  #   'INTERDOSE INTERVAL' = II,
  #   'STEADY STATE' = SS
  # ) %>%
  
##########################################################################-
# Escritura de archivos de datos
write_delim(x = data1, path = './RESULTADOS/Normal.csv', 
            delim = ';', na = '.')

# Escritura de archivos de datos
write_delim(x = data1, path = './RESULTADOS/Normal.txt', 
            delim = '\t', na = '.')


##########################################################################-
# Archivo de datos sin datos censurados -----------------------------------
##########################################################################-
# Eliminar los eventos de observación que sean menores a LLOQ
data2 <- data1 %>% 
  filter(!(DV <= LLOQ & (EVID == 0)))
##########################################################################-
# Escritura de archivos de datos
write_delim(x = data2, path = './RESULTADOS/Eliminado.csv', 
            delim = ';', na = '.')
write_delim(x = data2, path = './RESULTADOS/Eliminado.txt', 
            delim = '\t', na = '.')


##########################################################################-
# Archivo de datos con interpretación censurada por izquierda -------------
##########################################################################-
data3 <- data1 %>%
  mutate(CENS = if_else(!(DV <= LLOQ & (EVID == 0)), 0, 1),
         DV = if_else(!(DV <= LLOQ & (EVID == 0)), DV, LLOQ))

##########################################################################-
# Escritura de archivos de datos
write_delim(x = data3, path = './RESULTADOS/left_censored.csv', 
            delim = ';', na = '.')

write_delim(x = data3, path = './RESULTADOS/left_censored.txt', 
            delim = '\t', na = '.')

##########################################################################-
# Archivo de datos con interpretación censurada por intervalo -------------
##########################################################################-
data4 <- 
  data1 %>%
  mutate(
    LIMIT = if_else(!(DV <= LLOQ & (EVID == 0)), '.', '0'),
    CENS = if_else(!(DV <= LLOQ & (EVID == 0)), 0, 1),
    DV = if_else(!(DV <= LLOQ & (EVID == 0)), DV, LLOQ)
  )

##########################################################################-
# Escritura de archivos de datos
write_delim(x = data4, path = './RESULTADOS/interv_censored.csv', 
            delim = ';', na = '.')

write_delim(x = data4, path = './RESULTADOS/interv_censored.txt', 
            delim = '\t', na = '.')












