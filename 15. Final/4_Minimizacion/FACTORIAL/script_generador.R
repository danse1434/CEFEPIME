##------------------------------------------------------------------------------#
## Nombre del Script: Generación de set de datos para evaluación de -------------
## convergencia mediante un diseño factorial completo 
##  
## Propósito del Script: generar set de datos con diferentes conjuntos de set de 
## datos, definiendo valores medio (M), bajo (L) y alto (H).
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  12-06-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
# Carga de paquetes
require(rlang)
require(tidyverse)

# Carga de directorio principal
setwd(file.path('C:', 'Users', 'Daniel', 'OneDrive', 'Documents', 
                '(Proyecto)_Estudio_PKPD', 'CEFEPIME', '15. Final', 
                '4_Minimizacion'))

#-------------------------------------------------------------------------------#
# Introducción -----------------------------------------------------
#-------------------------------------------------------------------------------#

# Carga de archivo de datos con parámetros estimados
populationParameters <- read_csv("BASE_MODEL/populationParameters.txt")

# Calcular valores bajos (L) y altos (H), definidos como 50% y 100%
pop_par <- populationParameters %>%
  rename(MV = value) %>%
  mutate(LV = 0.5 * MV, HV = 1.5 * MV)

A <- pop_par %>%
  filter(parameter %in% c('Cl_pop', 'V1_pop', 'V2_pop')) %>%
  select(parameter, LV, MV, HV) %>%
  pivot_longer(cols = c(LV, MV, HV),
               names_to = 'tipo',
               values_to = 'valores') %>%
  pivot_wider(names_from = parameter, values_from = valores)

A

B <- expand.grid(A$Cl_pop, A$V1_pop, A$V2_pop)
colnames(B) <- colnames(A)[-1]

# Matriz con re-arreglos
B

#-------------------------------------------------------------------------------#
# Creación de archivos de inicio de evaluación ----------------------------------
#-------------------------------------------------------------------------------#
fileName <- '1_M_Error_1.mlxtran'
Z = readChar(fileName, file.info(fileName)$size)


# Quitar las tareas que no sean necesarias
Z1 <- Z %>%
  str_replace_all('\\r\\nindividualParameters\\(.+\\}\\)', '') %>% 
  str_replace_all('\\r\\nfim\\(.+\\)', '') %>% 
  str_replace_all('\\r\\nlogLikelihood\\(.+\\)', '') %>% 
  str_replace_all('\\r\\nplotResult\\(.+\\}\\)', '') %>% 
  str_replace_all('data/(?=1_data)', '')

for (i in 1:dim(B)[1]) {
  
  Y <- Z1 %>%
    str_replace("(?<=Cl\\_pop\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)",
                paste(B[i, 1])) %>%
    str_replace("(?<=V1\\_pop\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)",
                paste(B[i, 2])) %>%
    str_replace("(?<=V2\\_pop\\s\\=\\s\\{value\\=)\\d+\\.\\d+(?=\\,\\smethod)",
                paste(B[i, 3])) %>% 
    str_replace_all("(?<=exportpath\\s\\=\\s\\')BASE_MODEL", paste0('A_', i))
  
  write_lines(Y, file.path('FACTORIAL',
                           paste0('A', i, '.mlxtran')), sep = '\n')
  
  file.copy(
    from = file.path('data/1_data_TSFD.csv'),
    to = file.path('FACTORIAL', '1_data_TSFD.csv'))
}
