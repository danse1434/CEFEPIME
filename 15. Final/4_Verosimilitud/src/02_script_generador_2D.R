##------------------------------------------------------------------------------#
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
##------------------------------------------------------------------------------#
# Carga de paquetes
require(rlang)
require(tidyverse)

#-------------------------------------------------------------------------------#
# Especificación de variación en condiciones iniciales --------------------
#-------------------------------------------------------------------------------#
# Apertura de archivo de datos de parámetros de modelo base
populationParameters <-
  read_csv(
    "../../04. Residual/Modelo_8/RES_M1/populationParameters.txt"
  )
#-------------------------------------------------------------------------------#
# Cálculo de valor mínimo (50%) y valor máximo (150%) respecto al valor 
# nominal estimado en el modelo base.
pop_par <- 
  populationParameters %>% 
  filter(parameter %in% c('Cl_pop', 'V1_pop')) %>% 
  select(value) %>% 
  magrittr::use_series(value)

pop_vector_1 <- pop_par[1] * seq(0.5, 1.5, length.out = 30)
pop_vector_2 <- pop_par[2] * seq(0.5, 1.5, length.out = 30)

# Grilla expandida con valor de Cl (13.5), y V1_pop (23.9)
pop_vec_df <- expand.grid(pop_vector_1, pop_vector_2)
  
#-------------------------------------------------------------------------------#
# Selección de parámetros a evaluar
par_eval = 'Cl_pop-Q_pop'
# Preparación de carpetas
dir.create(file.path(par_eval))
dir.create(file.path(par_eval, 'assessment'))

#-------------------------------------------------------------------------------#
# Creación de archivos individuales ---------------------------------------
#-------------------------------------------------------------------------------#
# Apertura del archivo de control
#................................................................................
#  1 Abrir archivo de texto de control de Monolix.
#  2 Modificar el archivo para ser leído como una string de R.
#  3 Reemplazar el nombre de archivo a leer por el control de Monolix por
#  el correspondiente en cada carpeta.
#  Se cambian los valores de parámetro inicial, se conserva la misma semilla
#  de simulación para cada iteración.
#  4 Crear carpetas para contener los archivos.
#  5 Almacenar el archivo en el directorio correspondiente con el nombre
#  de "BASE_NEW.mlxtran", con formato de texto.
#................................................................................

fileName <- 'MODELO_BASE.mlxtran'
Z = readChar(fileName, file.info(fileName)$size)

for (i in 1:dim(pop_vec_df)[1]) {
  dir.create(file.path(par_eval, 'assessment', paste0('A', i)))

  Y <- Z %>%
    str_replace("(?<=Cl\\_pop\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)",
                paste(pop_vec_df[i,1])) %>% 
    str_replace("(?<=V1\\_pop\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)",
                paste(pop_vec_df[i,2]))
  
  write_lines(Y,
              file.path(par_eval, 'assessment',
                        paste0('A', i),
                        'MODELO_BASE.mlxtran'),
              sep = '\n')

  file.copy(
    from = file.path('data.csv'),
    to = file.path(par_eval, 'assessment', paste0('A', i)))
}

