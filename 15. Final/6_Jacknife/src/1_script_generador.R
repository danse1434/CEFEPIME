##========================================================================#
## Nombre del Script: análisis de individuos atípicos mediante jacknife
##  construcción de set de datos ------------------------------------------
##  
## Proposito del Script: preparación de carpetas para análisis por Jacknife 
##  con el fin de evaluar la influencia de cada paciente en la estimación de 
##  parámetros poblacionales.
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion: 08-03-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##========================================================================#
# Carga de paquetes
require(rlang)
require(tidyverse)



populationParameters <-
  read_csv(
    "../4_Minimizacion/BASE_MODEL/populationParameters.txt"
  )

 # Creación de directorio de evaluación
dir.create(file.path(getwd(), 'assessment'))

##########################################################################-
# Lectura de modelo base --------------------
##########################################################################-
# Apertura de archivo de datos del modelo base
data <- read_delim("data/1_data_TSFD.csv",";",  na = c('.'),
                   escape_double = FALSE, trim_ws = TRUE)

id_vector <- unique(data$ID)


for (i in 1:length(id_vector)) {
  dir.create(file.path(getwd(), 'assessment', paste0('J', i)))
  
  jack_data <- data %>%
    filter(ID != id_vector[i])

  write_delim(x = jack_data, 
              path = file.path(getwd(), 'assessment', paste0('J', i), 'data.csv'),
              delim = ";", na = '.')
  
  file.copy(
    from = file.path('data/3_M_Error_1.mlxtran'),
    to = file.path(getwd(), 'assessment', paste0('J', i)))
}
