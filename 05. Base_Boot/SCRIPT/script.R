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
# Modificación del archivo de datos ---------------------------------------
##########################################################################-
# Especificar la semilla de simulación
set.seed(12356)
# Lectura de archivo de datos originales
data <- read_delim("interv_censored.csv", delim = ";")
##########################################################################-
# Conversión de archivo de datos originales a una lista fragmentada con los 
# datos originales
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Agrupar la tabla *data* por ID
##  2 Separación de la tabla por el grupo como elementos de una lista
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
data.ls <- data %>%
  group_by(ID) %>%
  group_split(.)

##########################################################################-
# Creación de lista con data.frames con los remuestreos correspondientes
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Creación de un objeto tipo lista vacío
##  2 Hacer un muestreo con reemplazo de 15 individuos de la lista data.ls, 
##  generar una lista de estos muestreos.
##  3 Convertir las listas remuestradas en un data frame con una columna 
##  adicional new_id que actúa como un nuevo índice de individuos alternativos 
##  a id (esto permite que no individuos repetidos se comporten como individuos 
##  nuevos). Almacenar estas tablas dentro de cada una de los 1000 subelementos 
##  de la lista.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
table.ls = list()

for (i in 1:1000) {
  table.ls[[i]] <- data.ls %>%
    sample(., 15, replace = T) %>%
    map_dfr( ~ as.data.frame(.x), .id = 'new_id')
}

##########################################################################-
# Preparación de directorios con las tablas remuestreadas -----------------
##########################################################################-
# Creación de directorios y asignación de los datos de escritura
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Crear de 1000 directorios dentro de la ruta "../BOOT/", con el 
##  código Data{i}.
##  2 Escribir cada uno de los elementos de la lista *table.ls* con la clave 
##  "data_{i}.csv" y ubicarlo dentro de cada uno de las carpetas creadas.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

for (i in 1:1000) {
  dir.create(file.path('..', 'BOOT', paste0('Data', i)))
  
  write_delim(
    x = table.ls[[i]],
    path = file.path('..', 'BOOT', paste0('Data', i), paste0('data_', i, '.csv')),
    delim = ';',
    na = '.'
  )
}

##########################################################################-
# Copiar los archivos RES_M1.properties en cada uno de los directorios creados

for (i in 1:1000) {
  file.copy(
    from = file.path('..', 'BOOT', 'RES_M1.properties'),
    to = file.path('..', 'BOOT', paste0('Data', i))
  )
}

##########################################################################-
# Modificación de archivos de control de Monolix --------------------------
##########################################################################-
# Apertura del archivo de control
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Abrir archivo de texto de control de Monolix.
##  2 Modificar el archivo para ser leído como una string de R.
##  2 Reemplazar el nombre de archivo a leer por el control de Monolix por 
##  el correspondiente en cada carpeta.
##  4 Almacenar el archivo en el directorio correspondiente con el nombre 
##  de "BASE_NEW.mlxtran", con formato de texto.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

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
