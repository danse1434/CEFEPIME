## --------------------------- -
## Nombre del Script: Transformación de Datos de formato NONMEM a formato ADAPT V
##
## Propósito del Script: cambia un archivo plano con información farmacocinética 
## en formato NONMEM a un archivo plano en formato ADAPT - V
##
## Autor: Daniel S. Parra González
##
## Fecha de creación: 2020-01-04
##
## Copyright (c) Daniel S. Parra, 2020
## Email: dsparrag@unal.edu.co
## 
## Carga de paquetes ------------------------------------------------------
require('tidyverse')
require('lubridate')
require('rlang')
require('grid')
require('data.table')

##########################################################################-
# Directorio Principal ----------------------------------------------------
##########################################################################-
setwd(dir = "F:/Documentos/(Proyecto)_Estudio_PKPD/CEFEPIME/02. Base/NONMEM")

##########################################################################-
# Carga y alistamiento de archivo de datos --------------------------------
##########################################################################-
data <- read_csv(file = 'DATA/NONMEM_Cefepime.csv', 
                 col_names = T, 
                 na = c('.'),
                 col_types = cols(EVID = col_factor(),
                                  MDV = col_factor(),
                                  SS = col_factor(),
                                  SEXF = col_factor(), 
                                  ANTU = col_factor(),
                                  LLP = col_factor(),
                                  LMP = col_factor()))
##########################################################################-
# Modificación del archivo de acuerdo a requerimientos establecidos previamente
## Especificación de formato fecha en columna correspondiente
## Aplicación fórmula Dubois & Dubois
## Creación de columna con fecha y hora acoplada
## Creación de columna con tipo de diagnóstico agregado

data = data %>% 
  mutate(DATE = dmy(DATE)) %>% 
  mutate(SCM2 = 0.007184*(WTKG^0.425)*(HCM^0.725)) %>% 
  mutate(DATEHOUR = dmy_hms(paste(data$DATE, data$TIME))) %>% 
  mutate(Dx = case_when(LLP == 1 ~ 'LLA',
                        LMP == 1 ~ 'LMA/LMC',
                        (LLP==0 & LMP==0) ~ 'Otros'))

##########################################################################-
# Selección de tipo de tiempo que se va a modelar
## Almacenar como expresión
t = "TSFD" %>%
  parse_expr(.) 
## Número de observaciones
nobs = unique(data$ID) %>% length()

##########################################################################-
# Componentes de Eventos de Administración --------------------------------
##########################################################################-
## Crea las filas de terminación de infusiones
### Selección de eventos de administración
### Creación de duración de infusión
### Selección de columnas para especificar los eventos de terminación de infusión
### Agrupar por ID, y luego tiempos
### Seleccionar filas en cada grupo por posición
### Desagrupar el data frame
### Recalcular TSFD para que corresponda con terminación de t. de inf.
### Resetear la tasa de infusión a cero por ser terminación
data_A = data %>% 
  filter(EVID ==1) %>% 
  mutate(DUR = AMT/RATE) %>% 
  select(ID, RATE, !!t, DUR) %>% 
  group_by(ID, !!t) %>% 
  filter(row_number()==1 | row_number()==n()) %>% 
  ungroup() %>% 
  mutate(!!t := !!t + DUR) %>% 
  mutate(RATE = 0)
##########################################################################-
## Creación de Tabla de Eventos de Administración EVID == 1
### Selección de eventos de administración
### Creación de duración de infusión
### Selección de columnas para especificar los eventos de terminación de infusión
### Insertar los tiempos de terminación de infusión en *data_A* en la tabla
### Ordenamiento de la tabla por las variables ID y luego TSFD
### Seleccionar sólo el ID, TSFD, RATE
data_ADM = data %>%
  filter(EVID == 1) %>%
  mutate(DUR = AMT / RATE) %>%
  select(ID, RATE, !!t, DUR) %>%
  bind_rows(., data_A) %>%
  arrange(., ID, !!t) %>%
  select(ID, !!t, RATE)
##########################################################################-
## Creación de lista con tablas de administración 
### Tabla con todos los eventos de administración
### Conversión a objeto Lista con índice por ID
### Eliminación de nombre de columnas en cada tibble
### Conversión a archivo tipo data.frame
ADM_List_ID = data_ADM %>%
  split(x = select(., -c('ID')), f = as.factor(.$ID), drop = T) %>%
  map(~ setNames(.x, NULL)) %>%
  map(~ as.data.frame(.x))
##########################################################################-
## Creación de parámetros de administración de tipo ADAPT
### Tabla con todos los eventos de administración
### Agrupar por ID
### Resumir por número de eventos por individuo, crear dos columnas dummys 1, y 0
### Conversión a objeto Lista con índice por ID
### Conversión a objeto tipo data.frame
### Eliminación de nombre de columnas en cada data.frame
### Transponer la matriz de administración
ADM_List_ID0 = data_ADM %>%
  group_by(ID) %>%
  summarise(C1 = 1, C2 = 0, C3 = n())  %>%
  split(x = select(.,-c('ID')),f = as.factor(.$ID),drop = T) %>%
  map( ~ as.data.frame(.x)) %>%
  map( ~ setNames(.x, NULL)) %>%
  map( ~ t(.x))

##########################################################################-
# Componentes de Eventos de Observación -----------------------------------
##########################################################################-
## Creación de Tabla de Eventos de Observación EVID == 0
### Se filtran sólo los datos de observaciones
### Se seleccionan las variables de interés ID, !!t, DV
data_OBS = data %>%
  filter(EVID == 0) %>%
  select(ID, !!t, DV)
##########################################################################-
## Creación de lista con tablas de observaciones
### Se crea una lista con data.frames de observaciones (Tiempo y DV), con índice 
### ID
### Conversión a objeto tipo data.frame
### Eliminación de nombre de columnas en cada data.frame
OBS_List = data_OBS %>%
  split(x = .[, c(expr_text(t), 'DV')],f = as.factor(.$ID), drop = T) %>%
  map(~ as.data.frame(.x)) %>%
  map(~ setNames(.x, NULL))

##########################################################################-
## Creación de parámetros de observación de tipo ADAPT
### Tabla con todos los eventos de observación
### Agrupar por ID
### Resumir por número de eventos por individuo, crear una columna dummy 1.
### Conversión a objeto Lista con índice por ID
### Eliminación de nombre de columnas en cada data.frame
### Transponer las matrices de observaciones
OBS_List0 = data_OBS %>%
  group_by(ID) %>%
  summarise(n1 = 1, n2 = n()) %>% 
  split(x = .[, c('n1', 'n2')], f = as.factor(.$ID), drop = T) %>%
  map( ~ setNames(.x, NULL)) %>%
  map( ~ t(.x))

##########################################################################-
# Componentes de Eventos de Identificación --------------------------------
##########################################################################-
### Tabla
### Agrupar por ID
### Selecionar las primeros filas por grupo
### Crear una nueva variable que coloque el prefijo 'Sujeto_' al ID
### Selecionar sólo la variable ID_NAME
### Se crea una lista con los nombres, con índice ID
### Conversión a objeto tipo data.frame
### Eliminación de nombre de columnas en cada data.frame
### Transponer las matrices de observaciones
### Eliminación de variable ID en cada tabla
data_ID_names = data %>%
  group_by(ID) %>%
  filter(row_number() == 1) %>%
  mutate(ID_NAME = str_c('Sujeto_', ID)) %>%
  select(ID_NAME) %>%
  split(x = select(.,-ID),f = as.factor(.$ID),drop = T) %>%
  map( ~ as.data.frame(.x)) %>%
  map( ~ setNames(.x, NULL)) %>%
  map( ~ t(.x)) %>%
  map( ~ .x[-1])

##########################################################################-
# Creación de lista maestra -  --------------------------------------------
##########################################################################-
### Creación de lista
List_Total = list()
### Formación de una lista total con los cinco componentes
### Escribir vectores con líneas, los vectores simples se quedan así
### los data.frames se colapasan a vectores tipo CSV
for (i in seq(1,nobs*5,5)) {
  j = ((i-1)/5)+1
  List_Total[i+0] = data_ID_names[j]
  List_Total[i+1] = ADM_List_ID0[j]
  List_Total[i+2] = ADM_List_ID[j] %>% 
    map(~apply(.x, 1, function(y){paste(y,collapse=",")})) %>% 
    map( ~ t(.x))
  List_Total[i+3] = OBS_List0[j]
  List_Total[i+4] = OBS_List[j] %>% 
    map(~apply(.x, 1, function(y){paste(y,collapse=",")})) %>% 
    map( ~ t(.x))
}

##########################################################################-
# Creación de archivo *.CSV -----------------------------------------------
##########################################################################-
# Abrir la conexión con el archivo tipo *.csv receptor de la información 
fileConn = file("output.csv")
### Seleccionar la lista
### Eliminar estructura de lista a vector
### Escribir las líneas 
List_Total %>% 
  unlist(.) %>% 
  writeLines(., sep = "\n", con = fileConn)
### Cerrar la conexión
close(fileConn)


















