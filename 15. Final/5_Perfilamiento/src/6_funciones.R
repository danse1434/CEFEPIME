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

extractor <- function(popVal, n = 100) {
  X = vector('list', length = n); 
  Y = vector('list', length = n); 
  # Z = vector('list', length = n);

  for (i in 1:n) {
    dir_m = file.path(aux_dir, popVal, paste0('A',i), 'BASE_MODEL')
    
    X[[i]] = read_csv(file.path(dir_m, 'LogLikelihood', 'logLikelihood.txt'), 
                      col_types = cols())
    
    Y[[i]] = read_csv(file.path(dir_m, 'populationParameters.txt'), 
                      col_types = cols())
    
    # Z[[i]] = read_csv(file.path(dir_m, 'ChartsData', 'Saem', 'CvParam.txt'),
    #                   col_types = cols())
  }
  
  X1 <- X %>% 
    map_dfr(~ as_tibble(.x), .id = 'Run') %>%
    pivot_wider(id_cols = 'Run', names_from = 'criteria', values_from = 'importanceSampling')
  
  Y1 <- Y %>% 
    map_dfr(~ as_tibble(.x), .id = 'Run') %>% 
    pivot_wider(id_cols = 'Run', names_from = 'parameter', values_from = 'value')
  
  # Z1 <- Z %>%
  #   map_dfr(~ as_tibble(.x), .id = 'Run') %>%
  #   rename('convI' = 'convergenceIndicator') %>%
  #   select('Run', 'convI')

  
  W = left_join(X1, Y1, by = 'Run') %>% 
    # left_join(Z1, by = 'Run') %>%
    mutate(LL1 = `-2LL` - min(`-2LL`))
  
  return(W)
}


#-------------------------------------------------------------------------------#
#' Función que toma el valor 3,84 que corresponde a chi^2 con df = 1 y 
#' alpha de 0.05
#'
#' @param data una tabla de datos específica
#' @param x parámetro poblacional variante 
#' @param y 
#' @param v0 
#' @param xmin 
#' @param xmax 
#'
#' @return Retornar los dos puntos donde se corta la recta

xval.func <- function(data, x, y, v0 = 3.84, 
                      xmin = NA_real_, xmax = NA_real_) {
  # Conversión a quosures
  x_quo = rlang::ensym(x);
  y_quo = rlang::ensym(y);
  # Extracción de vectores **x** y **y**
  x_range = pull(data,!!x_quo);
  y_range = pull(data,!!y_quo);
  
  # Selección de valores mínimos o máximo por defecto
  x_min = max(xmin, min(x_range), na.rm = TRUE);
  x_max = min(xmax, max(x_range), na.rm = TRUE);
  
  # Valor de x con mínimo de LL
  m <- slice(data, which.min(!!y_quo)) %>%
    pull(!!x_quo) 
  
  # Optimización de límite de v0 para puntos de 
  # Crear una función aproximada con el valor de perturbación como x y la 
  # verosimilitud corregida como y.
  f1 <- approxfun(x = x_range, y = y_range)
  
  # Extraer los puntos donde se corta el intersecto optimizando las funciones 
  # y definiendo dos intervalos que limitan por izq y der al valor medio.
  A <- vector(length = 2)
  
  A[[1]] <- optimize(
    function(t0) abs(f1(t0) - v0),
    interval = c(x_min, x_max), upper = m)$minimum
  
  A[[2]] <- optimize(
    function(t0) abs(f1(t0) - v0),
    interval = c(x_min, x_max), lower = m)$minimum
  
  Minimo  = tibble(!!x_quo := m, !!y_quo := f1(m))
  p_corte = tibble(!!x_quo := A, !!y_quo := f1(A))
  
  Result = list('Minimo' = Minimo, 'Punto_Corte' = p_corte)
  return(Result)
}

