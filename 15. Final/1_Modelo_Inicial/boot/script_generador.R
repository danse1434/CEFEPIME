##------------------------------------------------------------------------------#
## Nombre del Script: Creación de data sets con instrucciones de bootstrap ------
## y almacenamiento en carpetas correspondientes
##  
## Propósito del Script: creación de muestreos sin reemplazo para obtener 
## estimación de intervalos de confianza de parámetros poblacionales.  
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  12-02-2020
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
require(glue)

#-------------------------------------------------------------------------------#
# Generación de subdirectorios de trabajos
set.seed(12356)

data <- read_delim('DATA/1_data_TSFD.csv', delim = ";", na = '.')

data.ls <- data %>%
  group_by(ID) %>%
  nest()

table.ls = list()

for (i in 1:1000) {
  table.ls[[i]] <- data.ls %>% 
    ungroup() %>% 
    sample_n(15, replace = TRUE) %>% 
    add_column(new_ID = 1:15) %>% 
    unnest(cols = c(data))
}

for (i in 1:1000) {
  dir.create(glue('boot/B{i}/'))
  
  write_delim(x = table.ls[[i]],
    path = glue('boot/B{i}/data_{i}.csv'), delim = ';', na = '.')
}

for (i in 1:1000) {
  file.copy(glue('boot/final_model.properties'),
            glue('boot/B{i}/'))
}

#-------------------------------------------------------------------------------#
# Apertura del archivo de control
fileName <- file.path('boot/final_model.mlxtran')
A = readChar(fileName, file.info(fileName)$size)

# Modificación de efectos aleatorios para mejorar la estimación SAEM
# Con esto se modifican los valores iniciales se dejan en 1.00 y esto permite una 
# mejor búsqueda en el espacio muestral. 

A1 <- A %>% 
  str_replace('(?<=CENS)\\}', ', new_ID}') %>% 
  str_replace('ID(?=\\s\\=\\s)', 'new_ID') %>% 
  str_replace_all('(?<=a\\s\\=\\s\\{value\\=)(\\-0|.)\\.\\d+(?=\\,)','1.00') %>%
  str_replace_all('(?<=b\\s\\=\\s\\{value\\=)(\\-0|.)\\.\\d+(?=\\,)','0.3') %>%
  str_replace_all('(?<=omega_Cl\\s\\=\\s\\{value\\=)(\\-0|.)\\.\\d+(?=\\,)','1.00') %>%
  str_replace_all('(?<=omega_Q\\s\\=\\s\\{value\\=)(\\-0|.)\\.\\d+(?=\\,)','1.00') %>%
  str_replace_all('(?<=omega_V1\\s\\=\\s\\{value\\=)(\\-0|.)\\.\\d+(?=\\,)','1.00') %>%
  str_replace_all('(?<=omega_V2\\s\\=\\s\\{value\\=)(\\-0|.)\\.\\d+(?=\\,)','1.00') %>% 
  str_replace('model_TBS', 'final_model') 


for (i in 1:1000) {
  B <- A1 %>%
    str_replace("(?<=data\\_)i(?=\\.csv)", glue('{i}'))
  
  write_lines(B,
              glue('boot/B{i}/final_model.mlxtran'), sep = '\n')
}
