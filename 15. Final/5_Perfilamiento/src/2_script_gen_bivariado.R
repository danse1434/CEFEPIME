##========================================================================#
## Nombre del Script: Perfilamiento de la función de verosimilitud (-2LL) #1 Preparación de archivos y carpetas de evaluación
##  
## Proposito del Script: este script tiene como fin realizar una preparación 
##  de archivos y carpetas para la evaluación de mapeo de verosimilitud 
##  univariada.
##
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion: 07-03-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##========================================================================#
# Carga de paquetes
require(rlang)
require(tidyverse)

##########################################################################-
# Especificación de variación en condiciones iniciales --------------------
##########################################################################-
# Apertura de archivo de datos de parámetros de modelo base
populationParameters <-
  read_csv(
    "../4_Minimizacion/BASE_MODEL/populationParameters.txt"
  )
##########################################################################-
# Selección de parámetros a evaluar
par_eval1 = 'V1_pop'
par_eval2 = 'beta_Cl_tSCRMGDL'

##########################################################################-
# Cálculo de valor mínimo (50%) y valor máximo (150%) respecto al valor 
# nominal estimado en el modelo base.
pop_par <- 
  populationParameters %>% 
  filter(parameter %in% c(par_eval1, par_eval2)) %>% 
  pull(value) 

pop_vector_1 <- pop_par[1] * seq(0.5, 1.5, length.out = 30)
pop_vector_2 <- pop_par[2] * seq(0.5, 1.5, length.out = 30)

# Grilla expandida con valor de Cl (13.5), y V1_pop (23.9)
pop_vec_df <- expand.grid(pop_vector_1, pop_vector_2)


aux_dir <- file.path('2_bivariado', paste0(par_eval1,'-',par_eval2))

# Preparación de carpetas
dir.create(file.path(aux_dir))

##########################################################################-
# Creación de archivos individuales ---------------------------------------
##########################################################################-
# Apertura del archivo de control
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Abrir archivo de texto de control de Monolix.
##  2 Modificar el archivo para ser leído como una string de R.
##  3 Reemplazar el nombre de archivo a leer por el control de Monolix por
##  el correspondiente en cada carpeta.
##  Se cambian los valores de parámetro inicial, se conserva la misma semilla
##  de simulación para cada iteración.
##  4 Crear carpetas para contener los archivos.
##  5 Almacenar el archivo en el directorio correspondiente con el nombre
##  de "BASE_NEW.mlxtran", con formato de texto.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
fileName <- '../4_Minimizacion/1_M_Error_1.mlxtran'

Z = readChar(fileName, file.info(fileName)$size)

# Quitar tareas innecesarias, sólo estimación de parámetros
Z1 <- Z %>% 
  str_replace_all('\\r\\nindividualParameters\\(.+\\}\\)', '') %>% 
  str_replace_all('\\r\\nfim\\(.+\\)', '') %>% 
  # str_replace_all('\\r\\nlogLikelihood\\(.+\\)', '') %>% 
  str_replace_all('\\r\\nplotResult\\(.+\\}\\)', '') %>% 
  # Sólo hay un archivo d edatos
  str_replace_all('data/(?=1_data)', '../')  %>% 
  # Lleva a cero todas las iteraciones del algoritmo SAEM
  str_replace('(?<=burniniterations\\s\\=\\s)\\d{1,6}', "0") %>% 
  str_replace('(?<=smoothingiterations\\s\\=\\s)\\d{1,6}', "0") %>% 
  str_replace('(?<=exploratoryiterations\\s\\=\\s)\\d{1,6}', "0") %>% 
  str_replace('(?<=simulatedannealingiterations\\s\\=\\s)\\d{1,6}', "0") %>% 
  str_replace('(?<=exploratoryinterval\\s\\=\\s)\\d{1,6}', "0") %>% 
  str_replace('(?<=smoothinginterval\\s\\=\\s)\\d{1,6}', "0")

str1 <-
  paste0('(?<=',
         par_eval1,
         '\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)')
str2 <- 
  paste0('(?<=',
         par_eval2,
         '\\s\\=\\s\\{value\\=)-\\d+\\.\\d+(?=\\,\\smethod)')
  
for (i in 1:dim(pop_vec_df)[1]) {
  dir.create(file.path(aux_dir,paste0('A', i)))
  
  Y <- Z1 %>%
    str_replace(str1,
                paste(pop_vec_df[i,1]) ) %>% 
    str_replace(str2,
                paste(pop_vec_df[i,2]))
  
  write_lines(Y, file.path(aux_dir,paste0('A', i), 'MODELO_FINAL.mlxtran'),
              sep = '\n')
  
}

file.copy(from = file.path('data/1_data_TSFD.csv'),
          to = file.path(aux_dir) )
