##------------------------------------------------------------------------#
## Nombre del Script: Revisión de trayectoria de algoritmos de evaluación 
##  covariables  ----------------------------------------------------
##  
## Propósito del Script: extraccion de datos sobre criterios de convergencia
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 19-05-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
# Selección de directorio de trabajo
setwd(file.path('C:', 'Users', 'Daniel', 'OneDrive', 'Documents', 
                '(Proyecto)_Estudio_PKPD', 'CEFEPIME', '14. Covariables_2', 
                '5_Modelo_SCM_LRT'))


#-------------------------------------------------------------------------------#
# Carga de paquetes
require(tidyverse)
require(ggplot2)
require(gridExtra)
require(gganimate)
require(patchwork)

#-------------------------------------------------------------------------------#
# Introducción ------------------------------------------------------------
#-------------------------------------------------------------------------------#
# Lectura de archivo de datos
# El archivo de datos original, se llama "arguments.dat" y tiene una 
# codificación desconocida, por lo cual se pasó a una hoja de cálculo, y se
# realizó un preprocesamiento con separación de las líneas por "|" y 
# conversión a un archivo "*.csv".

data <- read_csv("BASE_MODEL/ModelSelection/.Internals/Argumentos_1.csv",
                 skip = 12)

# Se eliminan las primeras 12 filas del archivo

#-------------------------------------------------------------------------------#
# Modificación de archivo de datos ----------------------------------------
#-------------------------------------------------------------------------------#

# Eliminación de espacios blancos en nombres de la tabla
dcol <- colnames(data) %>% 
  str_replace_all(" ", "") # Remoción espacios

# Desplazamiento una celda a la derecha de los nombres
for (i in 2:29) {
  colnames(data)[1] <- 'Intro' 
  colnames(data)[i] <- dcol[i-1]
}

#-------------------------------------------------------------------------------#
# Creación de objeto *data1*
#................................................................................
#  1 Detección de filas con letras en la columna X1
#  2 Renombrar a X1 como la variable Parametro
#  3 Remoción de columna Intro
#  4 Eliminación de espacios blancos en la variable parámetro
#  5 Seleccionar casos en "Parámetro" que correspondan a un param. PK
#................................................................................

data1 <- data %>% 
  filter(str_detect(X1, "\\w")) %>% 
  rename(Parametro = X1) %>% 
  select(-Intro) %>% 
  mutate(Parametro = str_replace_all(Parametro, "\\s", "")) %>% 
  filter(str_detect(Parametro, "Cl|V1|Q|V2")) 

#-------------------------------------------------------------------------------#
# Modificación a objeto *data2*
#................................................................................
#  1 Añadir columna con valores de iteración repetiendo 4 veces (por cada 
#  parámetro), la secuencia que indica la iteración.
#  2 Seleccionar covariables evaluadas en el assesment
#  3 Cambiar por 0 aquellos lugares donde exista NA en la tabla (covariable 
#  ausente para el parámetro), y 1 donde no (covariable presente para el 
#  parámetro).
#  4 Colapsar la tabla por covariables y cumplimiento de criterio, se tienen 
#  ahora dos columnas que indican el para cov-param involucrado por 
#  iteración.
#  5 Concatenar las columnas Parámetro y Covariable en el *par* involucrado
#................................................................................

data2 <- data1 %>% 
  add_column(Iteracion = rep(seq(1, dim(data1)[1]/4, 1), each = 4), 
             .before = "Parametro") %>% 
  select(matches("logt|ANTU|LLP|LMP|SEXF|Iteracion|Parametro|tSCRMGDL")) %>% 
  mutate_at(vars(matches("logt|ANTU|LLP|LMP|SEXF|tSCRMGDL")), 
            ~ ifelse(is.na(.), 0, 1)) %>% 
  pivot_longer(cols = matches("logt|ANTU|LLP|LMP|SEXF|tSCRMGDL"), 
               names_to = "Covariable", 
               values_to = "Criterio") %>% 
  unite(Parametro, Covariable, col = "Par", sep = '_')

#-------------------------------------------------------------------------------#
# Gráfico con pares estudiados ----------------------------------------------
#-------------------------------------------------------------------------------#
# Modificación de objeto a *data3*
#................................................................................
#  1 Modificación de *Par* a factor.
#  2 Creación de *Par1* que anonimiza los factores en enteros.
#  3 Modificación de *Criterio* a factor.
#  4 Separación de Par a dos columnas teniendo en cuenta el separador "_"
#................................................................................

data3 <- data2 %>%
  mutate(Par = factor(Par),
         Par1 = as.integer(fct_anon(Par)),
         Criterio = factor(Criterio)) %>% 
  separate(Par, into = c('Parametro', 'Covariable'), sep = "\\_", 
           remove = FALSE)

#-------------------------------------------------------------------------------#
# Gráfico con trayectorias de convergencia --------------------------------
#-------------------------------------------------------------------------------#
# Modificación de objeto a *data4*
#................................................................................
#  1 Selección de la primera columna
#  2 Filtrar filas que contengan las expresiones "LL" o "BI", estas 
#  tienen los criterios de convergencia alcanzados en cada paso.
#  3 Eliminar espacios blancos en columna
#  4 Separar la variable en *Parametro* y *Valor*
#  5 Convertir la variable *valor* en número
#  6 Adicionar una columna que indique la iteración, como se eliminó la 
#  primera iteración con el comando skip en la lectura se arranca desde 2.
#  7 Convertir los valores `-2*LL` en LRT en la variable *Parametro*
#................................................................................

data4 <- data %>% 
  select(Intro) %>% 
  filter(str_detect(Intro, "LL|BI")) %>% 
  mutate(Intro = str_replace_all(Intro, "\\s", "")) %>% 
  separate(Intro, into = c('Parametro', 'Valor'), "\\:") %>% 
  mutate(Valor = as.double(Valor)) %>% 
  add_column(Iteracion = rep(seq(2, (dim(.)[1] / 2)+1, 1), each = 2), 
             .before = "Parametro") %>% 
  mutate(Parametro = if_else(Parametro == "BICc", 'BICc', 'LRT'))

#-------------------------------------------------------------------------------#
# Creación de un archivo único con parámetros y valor de criterios --------
#-------------------------------------------------------------------------------#
# Creación de objeto *data6* 
#................................................................................
#  1 Colocación de valor de cada de criterio en su columna
#  2 Unión de tabla de criterios con tabla padre de parámetros por iteración
#................................................................................

data5 <- data4 %>% 
  pivot_wider(names_from = Parametro, values_from = Valor)

data6 <- data3 %>% 
  left_join(data5, by = 'Iteracion')


# Producto almacenado
write_csv(data6, "resumen_convergencia.csv")
