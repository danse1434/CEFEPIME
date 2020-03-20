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
# Definición de directorio principal
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '07. Minimizacion', '3_Mapeo'))

##########################################################################-
# Especificación de variación en condiciones iniciales --------------------
##########################################################################-
# Apertura de archivo de datos de parámetros de modelo base
populationParameters <-
  read_csv(
    "../../04. Residual/Modelo_8/RES_M1/populationParameters.txt"
  )

##########################################################################-
# Selección de parámetro a evaluar
par_eval = 'Q_pop'


##########################################################################-
# Cálculo de valor mínimo (50%) y valor máximo (150%) respecto al valor 
# nominal estimado en el modelo base.
pop_par <- 
  populationParameters %>% 
  filter(parameter == par_eval) %>% 
  select(value) %>% 
  magrittr::use_series(value)

pop_vector <- pop_par * seq(0.5, 1.5, length.out = 100)
  
##########################################################################-
# Preparación de carpetas
dir.create(file.path(par_eval))
dir.create(file.path(par_eval, 'assessment'))

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
fileName <- 'MODELO_BASE.mlxtran'
Z = readChar(fileName, file.info(fileName)$size)

for (i in 1:length(pop_vector)) {
  dir.create(file.path(par_eval, 'assessment', paste0('A', i)))

  Y <- Z %>%
    str_replace("(?<=Cl\\_pop\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)",
                paste(pop_vector[i]))
  
  write_lines(Y,
              file.path(par_eval, 'assessment', paste0('A', i),
                        'MODELO_BASE.mlxtran'),
              sep = '\n')

  file.copy(
    from = file.path('data.csv'),
    to = file.path(par_eval, 'assessment', paste0('A', i)))
}


# ##########################################################################-
# # Lectura de archivo verosimilitud ----------------------------------------
# ##########################################################################-
# # Creación de lista vacía pre-alocada
# LL_ls <- vector(mode = 'list', length = 100)
# 
# # Apertura de tablas en cada uno de los directorios de manera iterativa
# for (i in 1:100) {
#   LL_ls[[i]] <-
#     read_csv(
#       file.path('assessment', paste0('A', i), 'MODELO_BASE', 
#                 'LogLikelihood', 'logLikelihood.txt'), col_types = cols())
# }
# 
# A <- data.frame(Run = as.character(1:100),
#                 Parameter = pop_vector)
# 
# # Selección de tema
# theme_set(theme_classic() +
#             theme(panel.border = element_rect(fill = NA, colour = 'black')))
# 
# 
# LL_df <- LL_ls %>%
#   map_dfr(~ as.data.frame(.x), .id = 'Run') %>%
#   as_tibble(.) %>% 
#   inner_join(A, by = c("Run")) 
# 
# min_LL_df <- LL_df %>%
#   filter(criteria == "-2LL") %>%
#   slice(which.min(importanceSampling))
# 
# 
# LL_df %>% 
#   filter(criteria == "-2LL") %>% 
#   ggplot(aes(x = Parameter, y = importanceSampling)) +
#   geom_line() +
#   geom_vline(xintercept = pop_par, lty = 'dashed', col = 'blue4') +
#   geom_point(data = min_LL_df, col = 'red3', shape = 8) +
#   xlab('Cl') + ylab('-2LL')
#   
