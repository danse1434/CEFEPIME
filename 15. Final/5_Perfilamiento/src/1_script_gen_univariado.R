##------------------------------------------------------------------------------#
## Nombre del Script: Perfilamiento de la función de verosimilitud (-2LL) -------
## 1 Preparación de archivos y carpetas de evaluación univariada
## 
## Propósito del Script: este script tiene como fin realizar una preparación 
##  de archivos y carpetas para la evaluación de mapeo de verosimilitud 
##  univariada.
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  07-03-2020 (12-06-2020)
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
    "../4_Minimizacion/BASE_MODEL/populationParameters.txt"
  )

#-------------------------------------------------------------------------------#
# Selección de parámetro a evaluar
par_eval = 'Cl_pop'
# Esta se ejecuto teniendo en cuenta los parámetros V1, Cl, V2, beta_Cl, para 
# beta_Cl_tSCRMGDL tener en cuenta un menos que debe agregarse en REGEX en línea 87.

#-------------------------------------------------------------------------------#
# Cálculo de valor mínimo (50%) y valor máximo (150%) respecto al valor 
# nominal estimado en el modelo base.
pop_par <- 
  populationParameters %>% 
  filter(parameter == par_eval) %>% 
  pull(value)

pop_vector <- pop_par * seq(0.5, 1.5, length.out = 100)
  
#-------------------------------------------------------------------------------#
# Preparación de carpetas
aux_dir <- file.path('1_univariado', par_eval)
# Creación de carpeta resultados parámetro
dir.create(aux_dir)

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
fileName <- '../4_Minimizacion/1_M_Error_1.mlxtran'

Z = readChar(fileName, file.info(fileName)$size)

# Quitar tareas innecesarias, sólo estimación de parámetros
Z1 <- Z %>% 
  str_replace_all('\\r\\nindividualParameters\\(.+\\}\\)', 
                  '\r\nindividualParameters(method = conditionalMode)') %>%
  str_replace_all('\\r\\nfim\\(.+\\)', 
                  '\r\nfim(method = StochasticApproximation)') %>% 
  # str_replace_all('\\r\\nlogLikelihood\\(.+\\)', '') %>% 
  str_replace_all('\\r\\nplotResult\\(.+\\}\\)', 
                  '\r\nplotResult(method = {saemresults})') %>% 
  str_replace_all('data/(?=1_data)', '../')  %>% 
  # Lleva a cero todas las iteraciones del algoritmo SAEM
  str_replace('(?<=burniniterations\\s\\=\\s)\\d{1,6}', "0") %>% 
  str_replace('(?<=smoothingiterations\\s\\=\\s)\\d{1,6}', "0") %>% 
  str_replace('(?<=exploratoryiterations\\s\\=\\s)\\d{1,6}', "0") %>% 
  str_replace('(?<=simulatedannealingiterations\\s\\=\\s)\\d{1,6}', "0") %>% 
  str_replace('(?<=exploratoryinterval\\s\\=\\s)\\d{1,6}', "0") %>% 
  str_replace('(?<=smoothinginterval\\s\\=\\s)\\d{1,6}', "0")
Z1
str <-
  paste0('(?<=',
         par_eval,
         '\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)')

for (i in 1:length(pop_vector)) {
  dir.create(file.path(aux_dir, paste0('A', i)))
  
  Y <- Z1 %>% 
    str_replace(str, paste(pop_vector[i]) )

  write_lines(Y,file.path(aux_dir, paste0('A', i), 'MODELO_FINAL.mlxtran'),
              sep = '\n')
}


file.copy(
  from = file.path('data/1_data_TSFD.csv'),
  to = file.path(aux_dir) )
