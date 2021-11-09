#-------------------------------------------------------------------------------#
#' Crea una lista con historial de dosificación
#'
#' @param DD (decimal) dosis diaria en mg/d.
#' @param ii (entero) intervalo entre dosis, en horas.
#' @param tinf (decimal) tiempo de infusión, en horas.
#' @param n (entero) número de dosis a simular
#' @param tfd (decimal) tiempo inicial a simular
#'
#' @return (lista) lista con nombres:
#'         - amount: vector con dosis individuales
#'         - time: vector con tiempos de dosificación
#'         - tinf: tiempo de infusión
#' 
#' @export
#'
#' @examples
#' listaTratamiento(200, 2, 4, 10)
#' listaTratamiento(30, 6, 4, 10, 2)
#' listaTratamiento(200, 24, 4, 10)
#' listaTratamiento(15*70, 8, 4, 10)
#' 
listaTratamiento <- function(DD, ii, tinf, n, tfd = 0) {
  
  ls_d <- rep(DD * ii/24, n)
  ls_t <- seq(tfd, tfd + (ii*(n - 1)), by = ii)
  
  res <- list(amount = ls_d, time = ls_t, tinf = tinf)
  return(res)
}



