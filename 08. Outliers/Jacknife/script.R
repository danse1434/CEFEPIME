##========================================================================#
## Nombre del Script: análisis de individuos atípicos mediante jacknife
##  análisis de datos ----------------------------------------------------#
##  
## Proposito del Script: análisis de influencia de individuos
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion: 08-03-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##========================================================================#
# Carga de paquetes
require(rlang)
require(tidyverse)
require(bigutilsr)

# Definición de directorio principal
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '08. Outliers', 'Jacknife'))

##########################################################################-
# Lectura de archivos de parámetros poblacionales --------------------
##########################################################################-
# Creación de una lista con los parámetros poblacionales de los estimados 
# de Jacknife.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Prealocar la lista *data_vector*
##  2 En cada elemento de la lista, colocar una tabla con los datos de 
##  parámetros poblacionales en cada carpeta.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data_vector = vector(mode = 'list', length = 15L)

for (i in 1:15) {
  data_vector[[i]] <- read_csv(
    file.path('assessment', paste0('J', i) , 'BASE_MODEL',
              'populationParameters.txt'), col_types = cols()
    )
}

##########################################################################-
# Conversión de la lista en una tabla con componentes adicionales 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Abrir el objeto data_vector
##  2 Convertir en una tabla con identificación de cada elemento como ID
##  3 Convertir en un tibble
##  4 Agrupar por parámetro
##  5 Convertir la variable sujeto en número
##  6 Calcular el valor estandarizado para cada parámetro
##  7 Desagrupar
##  8 Modificar la variable parameter como un factor ordenado
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data_df <- data_vector %>% 
  map_dfr( ~ as.data.frame(.x), .id = 'Sujeto') %>%
  as_tibble(.) %>% 
  group_by(parameter) %>% 
  mutate(Sujeto = as.double(Sujeto)) %>% 
  mutate(std_value = (value - mean(value)) / sd(value)) %>% 
  ungroup() %>% 
  mutate(parameter = 
           factor(parameter, 
                  levels = c('Cl_pop', 'V1_pop', 'Q_pop', 'V2_pop',
                             'omega_Cl', 'omega_V1', 'omega_Q', 'omega_V2',
                             'a', 'b'))) #%>% 
  # filter(!(Sujeto %in% c(11)))
 
##########################################################################-
# Gráficos iniciales -----------------------------------------------------
##########################################################################-
# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

# Gráfico de cambio en cada parámetro en forma univariada 
G1 <- data_df %>% 
  ggplot(aes(x = Sujeto, y = std_value, 
             fill = ifelse(abs(std_value) > 1.960, 'out', 'in'))) +
  geom_hline(yintercept = 0) +
  geom_text(data = filter(data_df, abs(std_value) > 1.960),
            aes(y = std_value * 1.1, label = paste0('ID ', Sujeto)), size = 2.4) +
  geom_bar(stat = "identity", col = 'black', size = 0.04) + 
  geom_hline(yintercept = c(-1.960, 1.960), lty = 'dashed') +
  facet_wrap( ~ parameter, ncol = 4) + 
  scale_fill_manual(values = c('green1', 'red')) +
  ylab("Valor Estandarizado") + xlab("Índice") +
  coord_cartesian(ylim = c(-3.3,3.3)) +
  theme(legend.position = "none")

ggsave(filename = "Figura/valores_estandarizados.pdf", plot = G1, 
       device = 'pdf', width = 5, height = 4)

##########################################################################-
# Análisis de componentes principales -------------------------------------
##########################################################################-
# Carga de paquetes adicionales
# require(ggfortify)

##########################################################################-
# Conversión de tabla original en una tabla con cada parámetro en una columna
# (forma expandida / cartesiana).
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Selección del objeto *data_df*
##  2 Selección de las columnas Sujeto, parameter, y std_values (valores 
##  estandarizados).
##  3 Expansión de los datos tomando como parameter el valor de las columnas 
##  y std_value como las respuestas.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data_dfs <- data_df %>% 
  select(Sujeto, parameter, std_value) %>% 
  spread(parameter, std_value) 

##########################################################################-
# Realización de análisis de componentes principales
# 
# Se utilizó la función prcomp de la librería stats, se elimina la variable 
# sujeto ya que sólo se aceptan las variables que componen las características 
# del sujeto. Se almacena todo en el objeto *pca1*.

pca1 <- data_dfs %>%
  select(-Sujeto, -Q_pop) %>%
  stats::prcomp(.)

# Selección de la matriz con las variables transformadas en los componentes 
# principales. *pca2*

(pca1$sdev)^2 / sum((pca1$sdev)^2)

require(ggfortify)
pca1 %>% autoplot(x = 1, y = 9)

pca2 <- pca1$x

##########################################################################-
# Aplicación de criterios de identificación de outliers -------------------
##########################################################################-
# Aplicación de criterio de distancia de al menos 6 desviaciones estándar 
# desde la media.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar la matriz
##  2 A cada columna (2), aplicar una función que identifique cuales valores 
##  dentro de la matriz cumplen la condición que el valor abs de la 
##  diferencia con la media son mayores a 6 veces sd.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

pca2 %>%
  apply(., 2, function(x)
    which(abs(x - mean(x)) > (6 * sd(x))))

pca2 %>%
  apply(., 2, function(x)
    which((abs(x - median(x)) / 
             mad(x)) > 6)) 

# critdf <- data.frame(x = seq(1, 6, 0.01))
# for (i in 1:dim(critdf)[1]) {
#   critdf[i, 'y'] <- pca2 %>%
#     apply(., 2, function(x)
#       which(abs(x - mean(x)) > (critdf[i, 'x'] * sd(x)))) %>%
#     length(.)
# }
# critdf %>%
#   ggplot(aes(x,y)) +
#   geom_point()

# El criterio no es capaz de identificar algún outlier variando el umbral de 
# desviación estándar, debido a que no considera ninguno de los datos outlier 
# o los identifica a todos como tal.

dist <- pca2 %>%
  apply(., 2, function(x)
    (abs(x - median(x)) / mad(x) > 6)) %>%
  apply(1, max)

##########################################################################-
# Aplicación de la distancia robusta de Mahalanobis
# 
# Se aplica la función dist_ogk del paquete "bigutilsr" que permite calcular 
# la distancia de Mahalanobis para cada uno de los sujetos. Esta medida 
# permite calcular la distancia entre un Punto P y una distribución D. 
# 
# Es una generalización multidimensional de la idea de medir cuantas desvest 
# esta lejos de la media de P. 

dist2 <- bigutilsr::dist_ogk(pca2)
qplot(dist, sqrt(dist2))

pval <- pchisq(dist2, df = 9, lower.tail = FALSE)
hist(pval)

# Corrección de Bonferroni
is.out <- (pval < (0.05 / length(dist2)))  


qplot(pca2[, 1], pca2[, 2], color = is.out, size = I(3)) + coord_equal()


##########################################################################-
# Aplicación de Factor de Outlier Local (LOF)
dist3 <- LOF(pca2, seq_k = c(1,4, 14))

qplot(dist2, dist3)

qplot(pca2[, 1], pca2[, 2], color = dist3, size = I(3)) + coord_equal() +
  scale_color_viridis_c(breaks = NULL)




##########################################################################-
# Gráficos finales PCA ----------------------------------------------------
##########################################################################-
# Gráfico que acepta las dimensiones en cada PCA
pcout <- function(data, x, y) {
  x <- rlang::enquo(x)
  y <- rlang::enquo(y)
  W <- data %>%
    as_tibble(.) %>%
    rownames_to_column(var = 'ID')
  
  G1 <- W %>%
    ggplot(aes(!!x, !!y, colour = is.out)) +
    geom_point() +
    geom_hline(yintercept = 0,
               lty = 'dashed',
               col = "grey50") +
    geom_vline(xintercept = 0,
               lty = 'dashed',
               col = "grey50") +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none")
  
  return(list(data = W, graph = G1))
}

##########################################################################-
# 
elipsogk <- function(data = pca2, which, alpha){
  
  Y_center <- bigutilsr::covrob_ogk(data)$center[which]
  Y_cov <- bigutilsr::covrob_ogk(data)$cov[which, which]
  Y_radius <- sqrt(qchisq(1 - (alpha / 2), df = ncol(data[, which])))
  
  ellipse_ogk <- data.frame(
  car::ellipse(center = Y_center, shape = Y_cov, radius = Y_radius,
               segments = 100, draw = FALSE))
  
  colnames(ellipse_ogk) <- colnames(pca2[,which])
  
  list(geom_polygon(data = ellipse_ogk, color = "red", fill = "red", 
                    alpha = 0.01, lty = 'solid'))
  }



gm1 <- pcout(pca2, PC1, PC2)$graph +
  elipsogk(pca2, c(1, 2), 0.10) +
  elipsogk(pca2, c(1, 2), 0.01) +
  elipsogk(pca2, c(1, 2), 0.05) +
  geom_point(aes(x = bigutilsr::covrob_ogk(pca2)$center[1], 
                 y = bigutilsr::covrob_ogk(pca2)$center[2]), 
             color = "red", size = 3, shape = 8) +
  coord_cartesian(xlim = c(-6.0, 6.0), ylim = c(-6.0, 6.0)) +
  geom_text(
    pcout(pca2, PC1, PC2)$data %>%
      add_column(is.out) %>%
      filter(is.out == TRUE),
    mapping = aes(label = ID),
    nudge_x = 0.5
  )

gm2 <- pcout(pca2, PC1, PC3)$graph +
  elipsogk(pca2, c(1, 3), 0.10) +
  elipsogk(pca2, c(1, 3), 0.01) +
  elipsogk(pca2, c(1, 3), 0.05) +
  geom_point(aes(x = bigutilsr::covrob_ogk(pca2)$center[1], 
                 y = bigutilsr::covrob_ogk(pca2)$center[2]), 
             color = "red", size = 3, shape = 8) +
  coord_cartesian(xlim = c(-6.0, 6.0), ylim = c(-6.0, 6.0)) +
  geom_text(
    pcout(pca2, PC1, PC3)$data %>%
      add_column(is.out) %>%
      filter(is.out == TRUE),
    mapping = aes(label = ID),
    nudge_x = 0.5
  )


gm3 <- pcout(pca2, PC2, PC3)$graph +
  elipsogk(pca2, c(2, 3), 0.10) +
  elipsogk(pca2, c(2, 3), 0.01) +
  elipsogk(pca2, c(2, 3), 0.05) +
  geom_point(aes(x = bigutilsr::covrob_ogk(pca2)$center[1], 
                 y = bigutilsr::covrob_ogk(pca2)$center[2]), 
             color = "red", size = 3, shape = 8) +
  coord_cartesian(xlim = c(-6.0, 6.0), ylim = c(-6.0, 6.0)) +
  geom_text(
    pcout(pca2, PC2, PC3)$data %>%
      add_column(is.out) %>%
      filter(is.out == TRUE),
    mapping = aes(label = ID),
    nudge_x = 0.5
  )

ggsave('Figura/ComponentePC1vsPC2.pdf', gm1, device = 'pdf', 
       width = 5, height = 4)

ggsave('Figura/ComponentePC1vsPC3.pdf', gm2, device = 'pdf', 
       width = 5, height = 4)

ggsave('Figura/ComponentePC2vsPC3.pdf', gm3, device = 'pdf', 
       width = 5, height = 4)

##########################################################################-
# Gráfico 3D --------------------------------------------------------------
##########################################################################-

W <- pcout(pca2, PC1, PC3)$data

rot <- seq(-90, 100, length.out = 20)
rot_index <- seq(1, 20, length.out = 20)

for (i in rot_index) {
  
png(filename = paste("Figura/grafico_3D/PC_", rot_index[i], ".png", sep = " "), 
    width = 465, height = 225, units = 'mm', res = 300)

plot3D::scatter3D(
  x = W$PC1, y = W$PC2, z = W$PC3,
  colvar = NULL, col = "blue", pch = 19, cex = 0.5,
  bty = 'b2', type = 'h', theta = rot[i], phi = 20,
  xlab = 'PC1', ylab = 'PC2', zlab = 'PC3', ticktype = "detailed",
  xlim = c(-6.0, 6.0), ylim = c(-6.0, 6.0), zlim = c(-6.0, 6.0)
)

plot3D::text3D(x = W$PC1 + 0.05, 
               y = W$PC2 + 0.05, z = W$PC3 + 0.5,  
       labels = W$ID, add = TRUE, colkey = FALSE, cex = 1)
dev.off()
}

i = 5



pdf(width = 5*1.8, height = 4*1.8)

plot3D::scatter3D(
  x = W$PC1, y = W$PC2, z = W$PC3,
  colvar = NULL, col = "blue", pch = 19, cex = 0.5,
  bty = 'b2', type = 'h', theta = rot[i], phi = 20,
  xlab = 'PC1', ylab = 'PC2', zlab = 'PC3', ticktype = "detailed",
  xlim = c(-6.0, 6.0), ylim = c(-6.0, 6.0), zlim = c(-6.0, 6.0)
)

plot3D::text3D(x = W$PC1 + 0.05, 
               y = W$PC2 + 0.05, z = W$PC3 + 0.5,  
               labels = W$ID, add = TRUE, colkey = FALSE, cex = 1)

par(xpd=TRUE)

dev.off()

##########################################################################-
# Crear una animación con fines ilustrativos
# library(animation)
# 
# saveGIF({
#   for (i in rot_index) {
#     plot3D::scatter3D(
#       x = W$PC1, y = W$PC2, z = W$PC3,
#       colvar = NULL, col = "blue", pch = 19, cex = 0.5,
#       bty = 'b2', type = 'h', theta = rot[i], phi = 20,
#       xlab = 'PC1', ylab = 'PC2', zlab = 'PC3', ticktype = "detailed",
#       xlim = c(-6.0, 6.0), ylim = c(-6.0, 6.0), zlim = c(-6.0, 6.0)
#     )
#     
#     plot3D::text3D(x = W$PC1 + 0.05, 
#                    y = W$PC2 + 0.05, z = W$PC3 + 0.5,  
#                    labels = W$ID, add = TRUE, colkey = FALSE, cex = 1)
#     # dev.off()
#   }
# })




##########################################################################-
# Animación en el paquete RGL
# library(rgl)
# 
# x <- W$PC1
# y <- W$PC2
# z <- W$PC3
# 
# Y_center <- bigutilsr::covrob_ogk(pca2)$center[1:3]
# Y_cov <- bigutilsr::covrob_ogk(pca2)$cov[1:3, 1:3]
# 
# rgl_add_axes <- function(x, y, z, axis.col = "black",
#                          xlab = "", ylab="", zlab="", show.plane = FALSE, 
#                          show.bbox = FALSE, bbox.col = c("#b469ff","black"))
# { 
#   
#   lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
#   # Add axes
#   xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
#   rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
#   rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
#   rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
#   
#   # Add a point at the end of each axes to specify the direction
#   axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
#                 c(0, 0, zlim[2]))
#   rgl.points(axes, color = axis.col, size = 3)
#   
#   # Add axis labels
#   rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
#             adj = c(0.5, -0.8), size = 2)
#   
#   # Add plane
#   if(show.plane){ 
#     xlim <- xlim/1.1; zlim <- zlim /1.1
#     rgl.quads( x = rep(xlim, each = 2), y = c(0, 0, 0, 0),
#                z = c(zlim[1], zlim[2], zlim[2], zlim[1]))
#   }
#   # Add bounding box decoration
#   if(show.bbox){
#     rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
#              emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
#              xlen = 3, ylen = 3, zlen = 3) 
#   }
# }
# rgl.open()
# rgl.bg(color = "white") 
# rgl_add_axes(x, y, z, show.bbox = TRUE, show.plane = FALSE,
#              xlab = 'PCA1', ylab = 'PCA2', zlab = 'PCA3')
# rgl.spheres(x, y, z, r = 0.2, color = "black")
# 
# ellips <- ellipse3d(Y_cov, 
#                     centre= Y_center, level = 0.95)
# 
# shade3d(ellips, col = "#D95F02", alpha = 0.1, lit = FALSE)
# wire3d(ellips, col = "#D95F02",  lit = FALSE)
# # plot3d(ellips, col = "red", alpha = 0.2, add = TRUE, box = FALSE)
# aspect3d(1,1,1)
# 
# 
# # rgl.snapshot(filename = "plot.png")
# rgl.postscript("plot.pdf",fmt="pdf")
# 
# movie3d(spin3d(axis = c(1, 0, 0)), duration = 3,
#         dir = getwd())
# 
# 
# 





