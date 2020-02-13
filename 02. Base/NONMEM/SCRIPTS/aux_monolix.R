## --------------------------- -
## Nombre del Script: Transformación de Datos de formato NONMEM a formato Monolix
##
## Propósito del Script: cambia un archivo plano con información farmacocinética 
## en formato NONMEM a un archivo plano en formato Monolix
##
## Autor: Daniel S. Parra González
##
## Fecha de creación: 2020-01-24
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
require('lixoftConnectors')

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
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
## 1 Especificación de formato fecha en columna correspondiente
## 2 Aplicación fórmula Dubois & Dubois
## 3 Creación de columna con fecha y hora acoplada
## 4 Creación de columna con tipo de diagnóstico agregado
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
## Número de individuos en el estudio
nobs = unique(data$ID) %>% length()

##########################################################################-
# Especificación de Límite inferior de cuantificación ---------------------
##########################################################################-
LLOQ = 50/2^5
## Reemplazar esta variable en lugares correspondientes
# data = data %>%
#   mutate(CENSORING = if_else(DV <= LLOQ, 1, 0)) %>% 
#   mutate(DV = if_else(DV <= LLOQ, LLOQ, DV))



##########################################################################-
# Modificación a formato Monolix ------------------------------------------
##########################################################################-
# Modificación a formato Monolix
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Eliminar columnas TIME y DATE
##  2 Renombrar varias columnas para cumplir con nombres de Monolix
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data1 <-  data %>%
  select(-TIME, -DATE, -CMT) %>%
  rename(TIME = !!t) %>% 
  mutate('YTYPE' = ifelse(EVID == 0, 1, NA_real_)) %>%
  mutate('ADM' = ifelse(EVID == 1, 1, NA_real_)) %>% 
  select(-TAD, -DATEHOUR)
  # rename(
  #   'EVENT ID' = EVID,
  #   'IGNORED OBSERVATION' = MDV,
  #   'OBSERVATION' = DV,
  #   'AMOUNT' = AMT,
  #   'INFUSION RATE' = RATE,
  #   'ADDITIONAL DOSES' = ADDL,
  #   'INTERDOSE INTERVAL' = II,
  #   'STEADY STATE' = SS
  # ) %>%
  
##########################################################################-
# Escritura de archivos de datos
write_delim(x = data1, path = './RESULTADOS/Monolix_data.csv', 
            delim = ';', na = '.')

# Escritura de archivos de datos
write_delim(x = data1, path = './RESULTADOS/Monolix_data.txt', 
            delim = '\t', na = '.')

##########################################################################-
# Correr en forma de lote -------------------------------------------------
##########################################################################-
initializeLixoftConnectors(software = 'pkanalix')

head1 <- c('ID', 'observation', 'mdv',	'evid', 'amount', 'rate', 'addl', 
           'ii', 'ss', 'time',  rep('ignore', 18), 'OBSERVATION ID', 
           'admId')

head2 <- c('ID', 'observation', 'mdv',	'evid', 'amount', 'rate', 'addl', 
           'ii', 'ss', 'time',  'catcov', 'cov', 'cov', 'cov', 'cov', 
           'cov', 'cov', 'cov', 'cov', 'cov', 'cov', 'cov', 'cov', 'catcov', 
           'catcov', 'catcov', 'catcov', 'catcov', 'observation id', 
           'ADMINISTRATION ID')

path1 <- file.path(getwd(),'RESULTADOS','Monolix_data.txt')

# Creación de un proyecto de PKanalix
newProject(data = list(dataFile = path1, 
                       headerTypes = head1,
                       observationTypes = 'continuous'))


# Selección de tarea NCA
# undebug(setNCASettings)

##########################################################################-
# Colocar las configuraciones de NCA
setNCASettings(administrationType = list("1" = "Intravenous"), 
               integralMethod = "LinLogTrapLinLogInterp",
               blqMethodBeforeTmax = "missing",
               lambdaRule = "adjustedR2")

##########################################################################-
# Almacenar el proyecto en archivo
saveProject(projectFile = './RESULTADOS/proyecto_1.pkx')


##########################################################################-
# Análisis por NCA --------------------------------------------------------
##########################################################################-
# Calcular NCA
runNCAEstimation()

##########################################################################-
# Seleccion parámetros NCA específicos
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Obtener parámetros de NCA individuales
##  2 Obtención de encabezado con datos
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
indivParams <- getNCAIndividualParameters('AUCINF_pred', 'Cmax', 'Cmin', 
                                          'Rsq')$parameters
indivParams %>% head(.) %>% print(.)

##########################################################################-
# Post-procesamiento de los datos
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Abrir paquete 'table1'
##  2 Configurar nombre de AUCINF
##  3 Configurar nombre de Cmax
##  4 Configurar nombre de Cmin
##  5 Configurar nombre de Rsq
##  6 Crear una tabla HTML configurado de manera específico en idioma espa-
##  ñol. 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
require('table1')
table1::label(indivParams$AUCINF_pred) <- 'AUC a Infinito'
table1::label(indivParams$Cmax) <- 'Cmax'
table1::label(indivParams$Cmin) <- 'Cmin'
table1::label(indivParams$Rsq) <- 'Rsq'

table1::table1( ~ AUCINF_pred + Cmax + Cmin + Rsq, data = indivParams, 
                overall = 'Total', 
                render.continuous = c("Media [CV%]" = "Mean [CV%]",
                                     "Mediana [Min, Max]" = "Median [Min, Max]"))


##########################################################################-
# Análisis Compartimental -------------------------------------------------
##########################################################################-
# Configurar modelo estructural para Análisis Compartimental
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar el modelo de infusión de dos compartimentos con parámetros
##  ClV1QV2
##  2 Correr estimación por análisis compartimental
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
setStructuralModel('lib:infusion_2cpt_ClV1QV2.txt')

runCAEstimation()

# Obtener tabla de resumen estadístico de parámetros individuales
getCAIndividualParameters() %>% 
  magrittr::use_series('statistics') %>%  
  print()


##########################################################################-
# Creación de archivo de datos con la DV transformado en logarítmo --------
##########################################################################-
data2 <- data1 %>% 
  mutate(DV = ifelse(EVID == 0, log(DV), NA))

##########################################################################-
# Escritura de archivos de datos
write_delim(x = data2, path = './RESULTADOS/Monolix_data_LOG.csv', 
            delim = ';', na = '.')











##########################################################################-
# Creación de tabla de administración corregida para 72 horas
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  Especificación de número de días de esquema a simular (ndias)
##  Especificación de intervalo entre dosis (ii)
##  Especificación de número de dosis a simular (addl)
##  Especificación de vector de tiempos de administración (tvec)
##  Especificación de vector de covariables a incluir (cov_vec)
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ndias = 3
ii = 8
addl = (24/ii) * ndias
tvec = seq(0,addl*ii,ii)
cov_vec = list("ID", 'DV', 'MDV', 'CMT', 'EVID', 'AMT', 'RATE', 'ADDL', 'II', 
               'SS', 'DUR', "SEXF", "AGEA", "WTKG", "HCM", "IMCKGM2", "SCM2", 
               "SCRMGDL", "CLCRMLMIN","PROGDL", "ALBGDL", "DIND", "DNFD", 
               "RAL", "RAN", 'ANTU', 'LLP', 'LMP', 'Dx')
# Creación de tabla con covariables fijas para cada paciente
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  Tabla original
##  Filtrar los eventos de administración de la tabla
##  Seleccionar el evento donde TAD es cero, última dosis
##  Calcular la duración de la infusión (DUR)
##  Seleccionar a ID y las covariables a incluir
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
data0 = data %>%
  filter(EVID == 1) %>%
  filter(TAD == 0) %>%
  mutate(DUR = AMT / RATE) %>%
  select(!!!cov_vec)
##########################################################################-
# Creación de una tabla con los eventos simulados de acuerdo a ii y ndias
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  Tabla original
##  Filtrar los eventos de administración de la tabla
##  Calcular la duración de la infusión (DUR)
##  Filtrar eventos reales para reemplazarlos, sólo queda la dosis inicial
##  Recalcular TAD teniendo en cuenta tvec
##  Añadir filas con tiempos simulados para cada ID 
##  Unir con la tabla de covariables 
##  Mezclar las covariables separadas en una sóla columna
##  Borras las columnas innecesarias
##  Eliminar aquellas columnas que tengan todavía valores de TSFD
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data1 =  data %>%
  filter(EVID == 1) %>%
  mutate(DUR = AMT / RATE) %>%
  filter(TAD >= 0) %>%
  mutate(TAD = TAD + addl * ii) %>%
  add_row(ID = rep(seq(1, nobs, 1), each = length(tvec)),
          TAD = rep(tvec, nobs)) %>%
  left_join(., data0, by = "ID", copy = T)

for (i in 2:length(x)) {
  data1 = data1 %>%
    mutate(!!parse_expr(x[[i]]) := coalesce(!!parse_expr(paste0(x[[i]], ".x")),
                                            !!parse_expr(paste0(x[[i]], ".y"))))
}

data1 = data1 %>%
  select(-matches(".+\\.[xy]")) %>%
  filter(is.na(TSFD)) 

##########################################################################-
# Creación de una tabla de observaciones corregida
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  Filtrar tabla de observaciones
##  Calcular la duración de infusiones 
##  Recalcular TAD teniendo en cuenta tvec
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
data2 = data %>%
  filter(EVID == 0) %>%
  mutate(DUR = AMT / RATE) %>%
  mutate(TAD = TAD + addl * ii)
##########################################################################-
# Creación de la tabla completa
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  Seleccionar tabla con observaciones construida
##  Modificar la variable TAD para desplazar por DUR
##  Cambiar la cantidad en estos registros a cero (terminación infusión)
##  Unir con datos de administración de inicio de infusión
##  Unir con datos de observaciones
##  Ordenar por ID y luego por tiempos
##  Eliminar cualquier dato duplicado por ID y TAD
##  Filtrar datos de TAD por si existen errores
##  Asignar objeto data a objeto data1
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
data1 = data1 %>%
  mutate(TAD = TAD + DUR) %>%
  mutate(RATE = 0) %>%
  bind_rows(., data1) %>%
  bind_rows(., data2) %>%
  arrange(., ID, !!t) %>%
  distinct(ID, TAD, .keep_all = TRUE) %>%
  filter(!(TAD == 'inf'))
##########################################################################-
data = data1
##########################################################################-
# Componentes de Eventos de Administración --------------------------------
##########################################################################-
## Crea las filas de terminación de infusiones
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Selección de eventos de administración
### Creación de duración de infusión
### Selección de columnas para especificar los eventos de terminación de infusión
### Agrupar por ID, y luego tiempos
### Seleccionar filas en cada grupo por posición
### Desagrupar el data frame
### Recalcular TSFD para que corresponda con terminación de t. de inf.
### Resetear la tasa de infusión a cero por ser terminación
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Selección de eventos de administración
### Creación de duración de infusión
### Selección de columnas para especificar los eventos de terminación de infusión
### Insertar los tiempos de terminación de infusión en *data_A* en la tabla
### Ordenamiento de la tabla por las variables ID y luego TSFD
### Seleccionar sólo el ID, TSFD, RATE
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
data_ADM = data %>%
  filter(EVID == 1) %>%
  mutate(DUR = AMT / RATE) %>%
  select(ID, RATE, !!t, DUR) %>%
  # bind_rows(., data_A) %>%
  arrange(., ID, !!t) %>%
  select(ID, !!t, RATE)
##########################################################################-
## Creación de lista con tablas de administración 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Tabla con todos los eventos de administración
### Conversión a objeto Lista con índice por ID
### Eliminación de nombre de columnas en cada tibble
### Conversión a archivo tipo data.frame
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
ADM_List_ID = data_ADM %>%
  split(x = select(., -c('ID')), f = as.factor(.$ID), drop = T) %>%
  map(~ setNames(.x, NULL)) %>%
  map(~ as.data.frame(.x))
##########################################################################-
## Creación de parámetros de administración de tipo ADAPT
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Tabla con todos los eventos de administración
### Agrupar por ID
### Resumir por número de eventos por individuo, crear dos columnas dummys 1, y 0
### Conversión a objeto Lista con índice por ID
### Conversión a objeto tipo data.frame
### Eliminación de nombre de columnas en cada data.frame
### Transponer la matriz de administración
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Se filtran sólo los datos de observaciones
### Se seleccionan las variables de interés ID, !!t, DV
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
data_OBS = data %>%
  filter(EVID == 0) %>%
  select(ID, !!t, DV)
##########################################################################-
## Creación de lista con tablas de observaciones
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Se crea una lista con data.frames de observaciones (Tiempo y DV), con índice 
### ID
### Conversión a objeto tipo data.frame
### Eliminación de nombre de columnas en cada data.frame
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
OBS_List = data_OBS %>%
  split(x = .[, c(expr_text(t), 'DV')],f = as.factor(.$ID), drop = T) %>%
  map(~ as.data.frame(.x)) %>%
  map(~ setNames(.x, NULL))

##########################################################################-
## Creación de parámetros de observación de tipo ADAPT
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Tabla con todos los eventos de observación
### Agrupar por ID
### Resumir por número de eventos por individuo, crear una columna dummy 1.
### Conversión a objeto Lista con índice por ID
### Eliminación de nombre de columnas en cada data.frame
### Transponer las matrices de observaciones
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Agrupar por ID
### Selecionar las primeros filas por grupo
### Crear una nueva variable que coloque el prefijo 'Sujeto_' al ID
### Selecionar sólo la variable ID_NAME
### Se crea una lista con los nombres, con índice ID
### Conversión a objeto tipo data.frame
### Eliminación de nombre de columnas en cada data.frame
### Transponer las matrices de observaciones
### Eliminación de variable ID en cada tabla
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
fileConn = file("aux3.csv")
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
### Seleccionar la lista
### Eliminar estructura de lista a vector
### Escribir las líneas 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
List_Total %>% 
  unlist(.) %>% 
  writeLines(., sep = "\n", con = fileConn)
### Cerrar la conexión
close(fileConn)


















