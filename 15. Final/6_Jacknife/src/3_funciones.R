#-------------------------------------------------------------------------------#
#' Función de gráfico PC tipo biplot
#'
#' @param data matriz con resultados PCA, cada columna corresponde a un PC
#' @param x Componente en el eje X
#' @param y Componente en el eje Y
#'
#' @return Lista con tabla W, y gráfico G1 (biplot PCA)
#' @export
#'
#' @examples
#' pcout(pca2, PC1, PC2)
#' 
pcout <- function(data, x, y) {
  # Quosures
  xq <- rlang::enquo(x); yq <- rlang::enquo(y)
  
  # Tibble con columna identificadora
  W <- as_tibble(data) %>% 
    rownames_to_column(var = 'ID')
  
  # Gráfico
  G1 <- ggplot(W, aes(!!xq, !!yq, colour = is.out)) +
    geom_point() +
    # Líneas intersección (0,0)
    geom_hline(yintercept = 0, lty = 'dashed', col = "grey50") +
    geom_vline(xintercept = 0, lty = 'dashed', col = "grey50") +
    # Puntos de color negro/azul (adentro/fuera)
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none")
  
  return(list(data = W, graph = G1))
}

#-------------------------------------------------------------------------------#
#' Función de elipse
#'
#' @param data matriz con resultados PCA, cada columna corresponde a un PC
#' @param which vector de componentes PC a seleccionar
#' @param alpha nivel de significancia
#'
#' @return Lista con geom de tipo polígono para adición
#' @export
#'
#' @examples
#' elipsogk(pca2, c(1, 2), 0.10)
#' 
elipsogk <- function(data = pca2, which, alpha){
  # Matriz de centro, covarianza, y radio
  Ycent <- bigutilsr::covrob_ogk(data)$center[which]
  Ycov  <- bigutilsr::covrob_ogk(data)$cov[which, which]
  Yrad  <- sqrt(qchisq(1 - (alpha / 2), 
                       df = ncol(data[, which])))
  # Data.Frame con elipse
  ellipse_ogk <- data.frame(
    car::ellipse(center = Ycent, shape = Ycov, radius = Yrad,
                 segments = 100, draw = FALSE))
  
  # Nombre de columna
  colnames(ellipse_ogk) <- colnames(data[,which])
  
  # Lista con polígono para adicionar a gráfico
  list(geom_polygon(data = ellipse_ogk, color = "red", fill = "red", 
                    alpha = 0.01, lty = 'solid'))
}
