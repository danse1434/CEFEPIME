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
# Definición de directorio principal
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '08. Outliers', 'Jacknife'))


# Creación de directorio de evaluación
dir.create(file.path('assessment'))

##########################################################################-
# Lectura de modelo base --------------------
##########################################################################-
# Apertura de archivo de datos del modelo base
data <- read_delim("interv_censored.csv",
                   ";", 
                   na = c('.'),
                   escape_double = FALSE,
                   trim_ws = TRUE)


id_vector <- unique(data$ID)


for (i in 1:length(id_vector)) {
  dir.create(file.path('assessment', paste0('J', i)))
  
  jack_data <- data %>%
    filter(ID != id_vector[i])

  write_delim(x = jack_data, 
              path = file.path('assessment', paste0('J', i), 'data.csv'),
              delim = ";", na = '.')
  
  file.copy(
    from = file.path('Control.mlxtran'),
    to = file.path('assessment', paste0('J', i)))
}
