#-------------------------------------------------------------------------------#
#' Función de conversión a matriz
#' 
#' @param df datos en formato compacto columna de id, tiempo, y Cpred
#' @param id_col (string) nombre de columna identificadora
#' @param names (string) variable que se convierte en columnas
#' @param values (string) variable que se convierte en celdas
#' 
#' @return matriz con individuos por columnas 
#' @export
#' 
conv_matriz <- function(df, id_col = time, names = id, values = ufCc) {
  id_col_q <- rlang::enquo(id_col)
  names_q  <- rlang::enquo(names)
  values_q <- rlang::enquo(values)
  
  df %>% 
    pivot_wider(data = ., 
                id_cols = !!id_col_q, 
                names_from = !!names_q, 
                names_prefix = 'ID',
                values_from = !!values_q) %>% 
    select(-!!id_col_q) %>% 
    as.matrix()
}


#-------------------------------------------------------------------------------#
#' Función de conversión a matriz con identificación del tiempo
#' 
#' @param df datos en formato compacto columna de id, tiempo, y Cpred
#' @param id_col (string) nombre de columna identificadora
#' @param names (string) variable que se convierte en columnas
#' @param values (string) variable que se convierte en celdas
#' 
#' @return matriz con individuos por columnas 
#' @export
#' 
conv_tibble <- function(df, id_col = time, names = id, values = ufCc) {
  id_col_q <- rlang::enquo(id_col)
  names_q  <- rlang::enquo(names)
  values_q <- rlang::enquo(values)
  
  df %>% 
    pivot_wider(data = ., 
                id_cols = !!id_col_q, 
                names_from = !!names_q, 
                names_prefix = 'ID',
                values_from = !!values_q)
}
