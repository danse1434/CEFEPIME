##------------------------------------------------------------------------#
## Nombre del Script: Creación de data sets con instrucciones de bootstrap 
## y almacenamiento en carpetas correspondientes
##  
## Proposito del Script: creación de muestreos sin reemplazo para obtener 
## estimación de intervalos de confianza de parámetros poblacionales.  
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion:  12-02-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
##########################################################################-
# Introducción ------------------------------------------------------------
##########################################################################-
# Carga de paquetes
require(tidyverse)
require(rlang)
require(grid)
require(gridExtra)
require(readxl)
##########################################################################-
# Selección de directorio de trabajo
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '05. Base_Boot', 'SCRIPT'))

# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

##########################################################################-
set.seed(12356)

data <- read_delim("interv_censored.csv", delim = ";")

data.ls <- data %>%
  group_by(ID) %>%
  group_split(.)
  

table.ls = list()


for (i in 1:1000) {
  table.ls[[i]] <- data.ls %>%
    sample(., 15, replace = T) %>%
    map_dfr( ~ as.data.frame(.x), .id = 'new_id')
}


for (i in 1:1000) {
  dir.create(file.path('..', 'BOOT', paste0('Data', i)))
  
  write_delim(
    x = table.ls[[i]],
    path = file.path('..', 'BOOT', paste0('Data', i), paste0('data_', i, '.csv')),
    delim = ';',
    na = '.'
  )
}


for (i in 1:1000) {
  file.copy(
    from = file.path('..', 'BOOT', 'RES_M1.properties'),
    to = file.path('..', 'BOOT', paste0('Data', i))
  )
}

##########################################################################-
# Apertura del archivo de control
fileName <- 'BASE_MODEL.mlxtran'
A = readChar(fileName, file.info(fileName)$size)

for (i in 1:1000) {
  B <- A %>%
    str_replace("(?<=data\\_)i(?=\\.csv)", paste(i))
  
  write_lines(B, file.path('..',
                           'BOOT',
                           paste0('Data', i),
                           'BASE_NEW.mlxtran'), sep = '\n')
  
}
