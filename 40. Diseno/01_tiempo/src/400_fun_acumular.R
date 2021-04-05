#' Acumular Data Frame
#'
#' @param dat data.frame original
#' @param var variable a acumular
#'
#' @return
#' @export
#'
#' @examples
#' 
accumulate_by <- function(dat, var) {
  var <- lazyeval::f_eval(var, dat)
  lvls <- plotly:::getLevels(var)
  dats <- lapply(seq_along(lvls), function(x) {
    cbind(dat[var %in% lvls[seq(1, x)], ], frame = lvls[[x]])
  })
  dplyr::bind_rows(dats)
}