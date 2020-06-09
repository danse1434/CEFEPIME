##------------------------------------------------------------------------------#
## Nombre del Script: Análisis de mapeo de la función objetivo
##  
## Propósito del Script: El mapeo de OFV permite conocer si se ha alcanzado 
## un mínimo global del modelo, además permite conocer IC no asintóticos 
## alrededor de los parámetros del modelo. En este script se analizan los 
## resultados obtenidos mediante la evaluación de verosimilitud con la suite 
## de Monolix. En este script se leen datos y se produce un gráfico con 
## funciones de optimización.

## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion: 07-03-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#

#-------------------------------------------------------------------------------#
# Introducción ------------------------------------------------------------
#-------------------------------------------------------------------------------#
# Carga de paquetes
require(rlang)
require(tidyverse)

# Definición de directorio principal
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '07. Minimizacion', '3_Mapeo'))

# Apertura de archivo de datos de parámetros de modelo base final
populationParameters <-
  read_csv(
    "../../04. Residual/Modelo_8/RES_M1/populationParameters.txt"
  )

#-------------------------------------------------------------------------------#
# Lectura de archivos de verosimilitud ------------------------------------
#-------------------------------------------------------------------------------#
# Función de extracción de datos de archivos de verosimilitud (*extractor*)
#' @param popVal Parámetro a extraer de los archivos de datos
#' @param data Tabla de datos con información de parámetros poblacionales
#' 
#' @return DataFrame con 
#................................................................................
#  1 Convertir *popVal* en quosure
#  2 Prealocar al vector *X* como una lista de 100 items
#  3 Crear *pop_par* - valor de parámetro de interés 
#    a Seleccionar el archivo *data*
#    b Filtrar el parámetro de interés *popVal*
#    c Seleccionar la columna valor
#    d Tomar el valor
#  4 Crear *pop_vector* - vector con perturbación de parámetros; multiplicar 
#  por factores de corrección en el intervalo (0.5, 1.5), 100 veces.
#  5 Crear *A*, data.frame con la columna Run que permite identificar la 
#  corrida, y la columna Parameter es *pop_vector* y corresponde a la 
#  perturbación aplicada en la corrida.
#  6 Agregar a cada elemento una tabla leída desde un directorio especificado 
#  por el parámetro, de manera iterativa.
#  7 Creación de la tabla *Y* que contiene las columnas IS, y el valor 
#  perturbado.
#    a Seleccionar X
#    b Convertir a *X* en un data.frame
#    c Convertir en tibble
#    d Seleccionar el criterio -2LL
#    e Renombrar a la variable Parameter en el nombre del parámetro
#    f Seleccionar columnas importanceSampling y la que contiene al parámetro
#................................................................................

extractor <- function(popVal, data = populationParameters) {
  popVal <- rlang::enquo(popVal)
  
  X = vector(mode = 'list', length = 100)
  
  pop_par <- data %>%
    filter(parameter == !!popVal) %>%
    select(value) %>%
    magrittr::use_series(value)
  
  pop_vector <- pop_par * seq(0.5, 1.5, length.out = 100)
  
  A <- data.frame(Run = as.character(1:100),
                  Parameter = pop_vector)
  
  for (i in 1:100) {
    X[[i]] = read_csv(file.path(paste(quo_name(popVal)), 'GRUPO',
                                paste0('LL', i, '.txt')), col_types = cols())
  }
  
  Y <- X %>%
    map_dfr(~ as.data.frame(.x), .id = 'Run') %>%
    as_tibble(.) %>%
    filter(criteria == "-2LL") %>%
    inner_join(A, by = c("Run")) %>%
    rename(!!quo_name(popVal) := Parameter) %>%
    select(importanceSampling, !!popVal)
    
  return(Y)
}

#-------------------------------------------------------------------------------#
# Creación de *data* que es un data.frame que contiene a los valores de 
# parámetros perturbados, y sus valores de verosimilitud alcanzados.
# 
# Se unen las columnas de manera lateral
data = bind_cols(extractor("Cl_pop"), extractor("V1_pop"), 
                 extractor("V2_pop"), extractor("Q_pop"),
                 extractor("omega_Cl"), extractor("omega_V1"),
                 extractor("omega_Q"), extractor("omega_V2"),
                 extractor("a"), extractor("b") )

# Modificación de nombres de columnas de data
#................................................................................
#  1 Crear vector C
#  2 Colocar en C una versión modificada de los nombres de columna, si en 
#  la posición contiene el nombre "importance" se toma la posición siguiente 
#  y se le adiciona '_LL', de resto se utiliza el valor en la posición.
#  3 Asignar los nombres de columna de data como C
#................................................................................

C <- vector()

for (i in 1:dim(data)[2]) {
  if (str_detect(colnames(data)[i], "importance") == TRUE) {
    C[i] <- paste0(colnames(data)[i + 1], '_LL')
    # C[i] <- paste0('LL')
  } else {
    C[i] <- colnames(data)[i]
    
  }
}

C -> colnames(data)

#-------------------------------------------------------------------------------#
# Modificación de valores leídos ------------------------------------------
#-------------------------------------------------------------------------------#
# Convertir a la tabla data en una tabla colapsada pero estructurada.
#................................................................................
#  1 Seleccionar la tabla *data*
#  2 Colocar los nombres de la tabla como columna, esto expresa la pos 
#  relativa dentro del intervalo de exploración univariada.
#  3 Conversión a datos colapsados a todas las columnas excepto la 
#  correspondiente a la columna ID.
#  4 Adicionar una columna con NA (AUXILIAR) que se llama *blanco*.
#  5 Separar la columna por parámetro por la barra al piso, si esta le es 
#  seguida de 'LL', si no tiene nada a la derecha se le adiciona el NA de 
#  la columna *blanco* a la derecha. Se separa en Parametro IZQ y Verosimil DER
#  6 Eliminar la columna AUXILIAR blanco
#  7 Reemplazar los valores NA en la columna Verosimil por la palabra Valor
#  8 Convertir en datos cartesianos a la columna Verosimil, y colocar resultados 
#  como valor, esto crea 3 columnas una que corresponde a la perturbación, y 
#  otra a la verosimilitud.
#  9 Transformar la variable Parametro en un factor ordenado. 
#................................................................................

data1 <- data %>%
  rownames_to_column() %>%
  gather(matches("\\_pop|omega\\_|^a|^b"),
         key = "Parametro",
         value = 'Valor') %>%
  add_column(Blanco = NA_character_) %>%
  separate(
    Parametro,
    sep = "\\_(?=LL)",
    into = c('Parametro', 'Verosimil'),
    fill = "right"
  ) %>%
  select(-Blanco) %>%
  mutate(Verosimil = replace_na(Verosimil, 'Valor')) %>%
  spread(Verosimil, Valor) %>% 
  mutate(Parametro = 
           factor(Parametro, 
                  levels = c('Cl_pop', 'V1_pop', 'Q_pop', 'V2_pop',
                             'omega_Cl', 'omega_V1', 'omega_Q', 'omega_V2',
                             'a', 'b')))

#-------------------------------------------------------------------------------#
# Creación de gráfico -----------------------------------------------------
#-------------------------------------------------------------------------------#

# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))
#-------------------------------------------------------------------------------#
# Calcular una columna *LL1* corregida para tomar como 0 al valor mínimo del 
# perfil de convergencia.
data1 <- data1 %>% 
  group_by(Parametro) %>% 
  mutate(LL1 = LL - min(LL))

#-------------------------------------------------------------------------------#
# Filtrar tabla *data1* con filas que contienen el valor mínimo de la 
# función de verosimilitud.
min_data1 <- data1 %>%
  group_by(Parametro) %>% 
  slice(which.min(LL1)) %>% 
  filter(!str_detect(Parametro, "omega\\_|^a|^b"))

##########################################################################-
# Filtar la tabla *populationParameters* con los parámetros de interés
pop_par = populationParameters %>%
  rename(Parametro = parameter) %>%
  filter(!str_detect(Parametro, "omega\\_|^a|^b"))

##########################################################################-
# Crear gráfico G1, con perfiles de verosimilitud, y líneas auxiliares
G1 <- data1 %>%
  filter(!str_detect(Parametro, "omega\\_|^a|^b")) %>%
  ggplot(aes(x = Valor, y = LL1)) +
  geom_line() +
  geom_point(data = min_data1, col = 'red3', shape = 8) +
  geom_vline(data = pop_par,
             aes(xintercept = value),
             lty = 'dashed',
             col = 'blue4') +
  facet_wrap( ~ Parametro, ncol = 2, scales = "free") +
  ylab('LL - min(LL)') +
  geom_hline(yintercept = 3.84, lty = 'dotted') 

G1  # Vista preliminar

#-------------------------------------------------------------------------------#
#' Función que toma el valor 3,84 que corresponde a chi^2 con df = 1 y 
#' alpha de 0.05
#' @param data un tabla de datos específica
#'
#' @return Retornar los dos puntos donde se corta la recta
#................................................................................
#  1 Asignar v0 como 3.84
#  2 Crear una función aproximada con el valor de perturbación como x y la 
#  verosimilitud corregida como y.
#  3 Crear el vector A
#  4 Extraer el punto medio del perfil como el punto mínimo de la función 
#  de verosimilitud.
#    a Seleccionar fila con el menor valor de LL1
#    b Desagrupar 
#    c Selecionar la columna valor
#    d Extraer el valor de eje x como punto medio
#  5 Extraer los puntos donde se corta el intersecto optimizando las funciones 
#  y definiendo dos intervalos que limitan por izq y der al valor medio.
#................................................................................

xval.func <- function(data) {
  v0 <- 3.84
  f1 <- approxfun(x = data$Valor, y = data$LL1)
  
  A <- vector()
  
  m <- data %>% 
    slice(which.min(LL1)) %>% 
    ungroup() %>% 
    select(Valor) %>% 
    magrittr::extract2("Valor")
  
  A[[1]] <- optimize(function(t0)
    abs(f1(t0) - v0), interval = range(data1$Valor), upper = m)$minimum
  
  A[[2]] <- optimize(function(t0)
    abs(f1(t0) - v0), interval = range(data1$Valor), lower = m)$minimum
  
  return(A)
}

limit <- function(data = data1, param, L) {
  param <- rlang::enquo(param)
  
  data %>%
    filter(Parametro == !!param) %>%
    xval.func(.) %>%
    magrittr::extract2(L)
}

limit = Vectorize(limit, vectorize.args = c('param'))


approxfunf <- function(param, x = Valor, y = LL1, data = data1) {
  param <- rlang::enquo(param)
  
  df = data %>%
    filter(Parametro == !!param)
  
  x = df %>% magrittr::use_series(Valor)
  y = df %>% magrittr::use_series(LL1)
  z = approxfun(x = x, y = y)
  return(z)
}

f_Cl <- approxfunf(param = 'Cl_pop')
f_Q <- approxfunf(param = 'Q_pop')
f_V1 <- approxfunf(param = 'V1_pop')
f_V2 <- approxfunf(param = 'V2_pop')

min_data2 <- min_data1 %>% 
  mutate(x1 = limit(param = Parametro, L = 1)) %>% 
  mutate(x2 = limit(param = Parametro, L = 2)) %>% 
  add_column(y = 3.84)
  

min_data2[min_data2["Parametro"] == "Cl_pop", "x2"] <- 13.0
min_data2[min_data2["Parametro"] == "V2_pop", "x2"] <- 20.0

min_data2[min_data2["Parametro"] == "V2_pop", "x1"] <- 14.3

G1 <- G1 +
  geom_point(data = min_data2,
             mapping = aes(x = x1, y = y),col = 'green1') + 
  geom_point(data = min_data2,
             mapping = aes(x = x2, y = y), col = 'green1')

##########################################################################-
# Almacenamiento de gráfico
ggsave(G1, filename = 'perfiles_verosimilitud.pdf', device = 'pdf', 
       width = 5, height = 4)



  
  

#-------------------------------------------------------------------------------#
# Gráficos bidimensionales Cl vs V1 ---------------------------------------
#-------------------------------------------------------------------------------#
# Selección de parámetros desde el valor poblacional estimado en el modelo
# base final.
#................................................................................
#  1 Seleccionar de *populationParameters*
#  2 Filtrarlosparámetros requeridos
#  3 Extraer la columna "value"
#................................................................................

#-------------------------------------------------------------------------------#
# Conversión a un vector con perturbación, se escogen 30 puntos en el intervalo 
# para generar hasta 900 puntos en la malla 3D.
# Grilla expandida con valor de Cl (13.5), y V1_pop (23.9)

pop_par_func <- function(var){
  
  if (length(var) != 2) {
    stop("Error: la longitud del vector de variables debe ser 2!")
  }
  
  pop_par <- populationParameters %>% 
    filter(parameter %in% c({{ var }})) %>% 
    pull(value)
  
  pop_vector <- list()
  
  pop_vector[[1]] <- pop_par[1] * seq(0.5, 1.5, length.out = 30)
  pop_vector[[2]] <- pop_par[2] * seq(0.5, 1.5, length.out = 30)
  
  pop_vec_df <- expand.grid(pop_vector[[1]], pop_vector[[2]]) %>%
    rownames_to_column(., var = "Run") %>% 
    as_tibble(.)
  
  colnames(pop_vec_df) <- c('Run', var)
  
  return(list(var = var,
              pop_par = pop_par,
              pop_vector = pop_vector,
              pop_vec_df = pop_vec_df))
}

#-------------------------------------------------------------------------------#
# Lectura de archivo de datos con grilla
reader_func <- function(var, design_df){
  
  if (length(var) != 2) {
    stop("Error: la longitud del vector de variables debe ser 2!")
  }

  v <- rlang::enquos(var)
  
  X = vector(mode = 'list', length = 900)

  for (i in 1:900) {
    X[[i]] = read_csv(file.path(paste0(var[1], '-', var[2]), 'GRUPO',
                                paste0('LL', i, '.txt')), col_types = cols())
  }

  Y <- X %>%
    map_dfr( ~ as.data.frame(.x), .id = 'Run') %>%
    as_tibble(.) %>%
    filter(criteria == "-2LL") %>%
    inner_join(design_df, by = c("Run")) %>%
    rename(LL = importanceSampling) %>%
    mutate(LL1 = LL - min(LL))

  Ymin <- Y %>%
    slice(which.min(LL1))
  
  return(list(v,
    gridvalues <- Y,
    minvalues <- Ymin))
  }

#-------------------------------------------------------------------------------#
# Cl_pop vs V1_pop
X <- pop_par_func(var = c('Cl_pop', 'V1_pop'))
Y <- reader_func(var = c('Cl_pop', 'V1_pop'),
                 design_df = X$pop_vec_df)

G2 <- 
  ggplot2::ggplot(Y[[2]], aes(Cl_pop, V1_pop, z = LL1)) + 
  ggplot2::geom_contour_filled(aes(fill = stat(level)), 
                               breaks = c(seq(0, 2400, by = 200))) +
  ggplot2::geom_point(data = Y[[3]], shape = 8, col = 'red') +
  ggplot2::geom_point(x = X$pop_par[1], y = X$pop_par[2], shape = 4, col = 'green1') +
  ggplot2::guides(fill = guide_colorsteps(barheight = unit(4, "cm"))) + 
  scale_fill_viridis_d(name = 'Nivel') +
  xlab(expression(V[1])) + ylab("Cl")

ggsave(filename = 'ClvsV1.pdf', plot = G2, device = 'pdf', 
       width = 6/1.5, height = 4/1.5)

ggsave(filename = 'ClvsV1.png', plot = G2, device = 'png', 
       width = 6/1.5, height = 4/1.5, dpi = 300)


##########################################################################-
# Cl_pop-Q_pop
X <- pop_par_func(var = c('Cl_pop', 'Q_pop'))
Y <- reader_func(var = c('Cl_pop', 'Q_pop'),
                 design_df = X$pop_vec_df)

G3 <- ggplot2::ggplot(Y[[2]], aes(Cl_pop, Q_pop, z = LL1)) + 
  ggplot2::geom_contour_filled(aes(fill = stat(level)), 
                               breaks = c(seq(0, 2400, by = 200))) +
  ggplot2::geom_point(data = Y[[3]], shape = 8, col = 'red') +
  ggplot2::geom_point(x = X$pop_par[1], y = X$pop_par[2], shape = 4, col = 'green1') +
  ggplot2::guides(fill = guide_colorsteps(barheight = unit(4, "cm"))) + 
  scale_fill_viridis_d(name = 'Nivel') +
  xlab("Cl") + ylab("Q")

ggsave(filename = 'ClvsQ.pdf', plot = G3, device = 'pdf', 
       width = 6/1.5, height = 4/1.5)

##########################################################################-
# V1_pop-V2_pop
X <- pop_par_func(var = c('V1_pop', 'V2_pop'))
Y <- reader_func(var = c('V1_pop', 'V2_pop'),
                 design_df = X$pop_vec_df)

G4 <- 
  ggplot2::ggplot(Y[[2]], aes(V1_pop, V2_pop, z = LL1)) + 
  ggplot2::geom_contour_filled(aes(fill = stat(level))) +
  ggplot2::geom_point(data = Y[[3]], shape = 8, col = 'red') +
  ggplot2::geom_point(x = X$pop_par[1], y = X$pop_par[2], shape = 4, col = 'green1') +
  ggplot2::guides(fill = guide_colorsteps(barheight = unit(4, "cm"))) + 
  scale_fill_viridis_d(name = 'Nivel') +
  xlab(expression(V[1])) + ylab(expression(V[2]))

ggsave(filename = 'V1vsV2.pdf', plot = G4, device = 'pdf', 
       width = 6/1.5, height = 4/1.5)


##########################################################################-
# V2_pop-Q_pop
X <- pop_par_func(var = c('Q_pop', 'V2_pop'))
Y <- reader_func(var = c('V2_pop', 'Q_pop'),
                 design_df = X$pop_vec_df)

G5 <- 
  ggplot2::ggplot(Y[[2]], aes(V2_pop, Q_pop, z = LL1)) + 
  ggplot2::geom_contour_filled(aes(fill = stat(level))) +
  ggplot2::geom_point(data = Y[[3]], shape = 8, col = 'red') +
  ggplot2::geom_point(x = X$pop_par[1], y = X$pop_par[2], shape = 4, col = 'green1') +
  ggplot2::guides(fill = guide_colorsteps(barheight = unit(4, "cm"))) + 
  scale_fill_viridis_d(name = 'Nivel') +
  ylab(expression(V[2])) + xlab("Q")

G5

ggsave(filename = 'V2vsQ.pdf', plot = G5, device = 'pdf', 
       width = 6/1.5, height = 4/1.5)





##########################################################################-
# V1_pop-Q_pop
X <- pop_par_func(var = c('V1_pop', 'Q_pop'))
Y <- reader_func(var = c('V1_pop', 'Q_pop'),
                 design_df = X$pop_vec_df)


G6 <-
  ggplot2::ggplot(Y[[2]], aes(V1_pop, Q_pop, z = LL1)) + 
  ggplot2::geom_contour_filled(aes(fill = stat(level))) +
  ggplot2::geom_point(data = Y[[3]], shape = 8, col = 'red') +
  ggplot2::geom_point(x = X$pop_par[1], y = X$pop_par[2], shape = 4, col = 'green1') +
  ggplot2::guides(fill = guide_colorsteps(barheight = unit(4, "cm"))) + 
  scale_fill_viridis_d(name = 'Nivel') +
  xlab(expression(V[1])) + ylab("Q")

ggsave(filename = 'V1vsQ.pdf', plot = G6, device = 'pdf', 
       width = 6 / 1.5, height = 4 / 1.5)









