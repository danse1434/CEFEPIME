##------------------------------------------------------------------------------#
## Nombre del Script: análisis de outliers mediante jacknife --------------------
##  
## Propósito del Script: análisis de influencia de individuos
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 08-03-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
# Carga de paquetes
require(rlang)
require(tidyverse)
require(patchwork)
require(bigutilsr)

# Lectura de archivo de funciones
source('src/3_funciones.R', encoding = 'UTF-8')

# Órden de parámetros
par_lev = c('Cl_pop', 'V1_pop', 'V2_pop', 'beta_Cl_tSCRMGDL',
            'omega_Cl', 'omega_V1', 'omega_Q', 'omega_V2', 'a')

#-------------------------------------------------------------------------------#
# Lectura de archivos de parámetros poblacionales --------------------
#-------------------------------------------------------------------------------#
# Creación de una lista con parámetros POP estimados con Jacknife.
#................................................................................
#  1 Prealocar la lista *data_vector*
#  2 En cada elemento de la lista, colocar una tabla con los datos de 
#  parámetros poblacionales en cada carpeta.
#................................................................................
data_vector = vector('list', 15L)

for (i in 1:15) {
  data_vector[[i]] <- read_csv(
    file.path('assessment', paste0('J', i) , 'BASE_MODEL',
              'populationParameters.txt'), col_types = cols()
    )
}

#-------------------------------------------------------------------------------#
# Conversión de la lista en una tabla con componentes adicionales 
#................................................................................
#  1 Abrir el objeto data_vector
#  2 Convertir en una tabla con identificación de cada elemento como ID
#  3 Convertir en un tibble
#  4 Agrupar por parámetro
#  5 Convertir la variable sujeto en número
#  6 Calcular el valor estandarizado para cada parámetro
#  7 Desagrupar
#  8 Modificar la variable parameter como un factor ordenado
#................................................................................

data_df <- data_vector %>% 
  map_dfr( ~ as_tibble(.x), .id = 'Sujeto') %>%
  group_by(parameter) %>% 
  mutate(Sujeto = as.integer(Sujeto),
         std_value = (value - mean(value)) / sd(value) ) %>% 
  ungroup() %>% 
  mutate(parameter = factor(parameter, levels = par_lev) ) %>% 
  filter(!is.na(parameter))
 
#-------------------------------------------------------------------------------#
# Gráficos iniciales -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Selección de tema
theme_set(theme_bw() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

# Gráfico de cambio en cada parámetro en forma univariada 
G1 <- data_df %>% 
  ggplot(aes(x = Sujeto, y = std_value, fill = parameter)) +
  geom_hline(yintercept = 0) +
  geom_bar(stat = "identity", col = 'black', size = 0.04) + 
  geom_hline(yintercept = c(-1.960, 1.960), lty = 'dashed') +
  facet_wrap(~parameter, ncol = 4) + 
  scale_fill_viridis_d() +
  ylab("Valor Estandarizado") + xlab("Índice") +
  coord_cartesian(ylim = c(-3.1,3.1)) +
  theme(legend.position = "none")
G1
ggsave("1_valores_estandarizados.pdf", G1, 'pdf', 'figures', 1, 5, 4)

#-------------------------------------------------------------------------------#
# Análisis de componentes principales -------------------------------------
#-------------------------------------------------------------------------------#
# Carga de paquetes adicionales
# require(ggfortify)

#-------------------------------------------------------------------------------#
# Conversión de tabla original en una tabla con cada parámetro en una columna
# (forma expandida / cartesiana).
#................................................................................
#  1 Selección del objeto *data_df*
#  2 Selección de las columnas Sujeto, parameter, y std_values (valores 
#  estandarizados).
#  3 Expansión de los datos tomando como parameter el valor de las columnas 
#  y std_value como las respuestas.
#................................................................................

data_dfs <- data_df %>%
  select(Sujeto, parameter, std_value) %>%
  pivot_wider(id_cols = Sujeto,
              names_from = 'parameter',
              values_from = 'std_value') 

#-------------------------------------------------------------------------------#
# Realización de análisis de componentes principales (PCA)
# 
# Se utilizó la función "prcomp" de la librería stats, se elimina la variable 
# sujeto ya que sólo se aceptan las variables que componen las características 
# del sujeto. Se almacena todo en el objeto *pca1*.

pca1 <- data_dfs %>%
  select(-Sujeto) %>%
  stats::prcomp(.)

# Selección de la matriz con las variables transformadas en los componentes 
# principales. *pca2*

(pca1$sdev)^2 / sum((pca1$sdev)^2)

require(ggfortify)
pca1 %>% autoplot(x = 3, y =4)
# Los cuatro primeros componentes explican 31.04%, 28.23%, 14.27%, y 8.51% de la 
# variabilidad del modelo (en total 82.05%).
pca2 <- pca1$x

#-------------------------------------------------------------------------------#
# Aplicación de criterios de identificación de outliers -------------------
#-------------------------------------------------------------------------------#
# Criterio de distancia de al menos 6 desviaciones estándar desde la media.
#................................................................................
#  1 Seleccionar la matriz
#  2 A cada columna (2), aplicar una función que identifique cuales valores 
#  dentro de la matriz cumplen la condición que el valor abs de la 
#  diferencia con la media son mayores a 6 veces sd.
#................................................................................

apply(pca2, 2, function(x)
  which(abs(x - mean(x)) > (6 * sd(x))))

apply(pca2, 2, function(x)
  which((abs(x - median(x)) / mad(x)) > 6)) 

# Ninguno de los individuos se puede descartar con estos criterios.

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

dist <- apply(pca2, 2, function(x)
    (abs(x - median(x)) / mad(x) > 6)) %>%
  apply(1, max)

#-------------------------------------------------------------------------------#
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

#-------------------------------------------------------------------------------#
# Aplicación de Factor de Outlier Local (LOF)
dist3 <- LOF(pca2, seq_k = c(1,4, 14))

qplot(dist2, dist3)

qplot(pca2[, 1], pca2[, 2], color = dist3, size = I(3)) + coord_equal() +
  scale_color_viridis_c(breaks = NULL)

#-------------------------------------------------------------------------------#
# Gráficos finales PCA ----------------------------------------------------
#-------------------------------------------------------------------------------#
# Gráficos con PCA para cada par (PC1,PC2)-(PC2,PC3)-(PC1,PC3)
# 
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
                 y = bigutilsr::covrob_ogk(pca2)$center[3]), 
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
  geom_point(aes(x = bigutilsr::covrob_ogk(pca2)$center[2], 
                 y = bigutilsr::covrob_ogk(pca2)$center[3]), 
             color = "red", size = 3, shape = 8) +
  coord_cartesian(xlim = c(-6.0, 6.0), ylim = c(-6.0, 6.0)) +
  geom_text(
    pcout(pca2, PC2, PC3)$data %>%
      add_column(is.out) %>%
      filter(is.out == TRUE),
    mapping = aes(label = ID),
    nudge_x = 0.5
  )

gm4 <- pcout(pca2, PC1, PC4)$graph +
  elipsogk(pca2, c(1, 4), 0.10) +
  elipsogk(pca2, c(1, 4), 0.01) +
  elipsogk(pca2, c(1, 4), 0.05) +
  geom_point(aes(x = bigutilsr::covrob_ogk(pca2)$center[1], 
                 y = bigutilsr::covrob_ogk(pca2)$center[4]), 
             color = "red", size = 3, shape = 8) +
  coord_cartesian(xlim = c(-6.0, 6.0), ylim = c(-6.0, 6.0)) +
  geom_text(
    pcout(pca2, PC1, PC4)$data %>%
      add_column(is.out) %>%
      filter(is.out == TRUE),
    mapping = aes(label = ID),
    nudge_x = 0.5
  )

# Almacenamiento de gráficos
ggsave('2_CompPC1vsPC2.pdf', gm1, 'pdf', 'figures', 1, 5, 4)
ggsave('2_CompPC1vsPC3.pdf', gm2, 'pdf', 'figures', 1, 5, 4)
ggsave('2_CompPC2vsPC2.pdf', gm3, 'pdf', 'figures', 1, 5, 4)

#-------------------------------------------------------------------------------#
# Gráfico 3D --------------------------------------------------------------
#-------------------------------------------------------------------------------#
W         <- pcout(pca2, PC1, PC3)$data
rot       <- seq(-90, 100, length.out = 20)
rot_index <- seq(1, 20, length.out = 20)

for (i in rot_index) {
  
png(filename = paste("figures/3D/PC_", rot_index[i], ".png", sep = " "), 
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

# Ángulo de rotación seleccionado
i = 10


pdf(width = 5*1.8, height = 4*1.8)

a <- {plot3D::scatter3D(
  x = W$PC1, y = W$PC2, z = W$PC3,
  colvar = NULL, col = "blue", pch = 19, cex = 0.5,
  bty = 'b2', type = 'h', theta = rot[i], phi = 20,
  xlab = 'PC1', ylab = 'PC2', zlab = 'PC3', ticktype = "detailed",
  xlim = c(-6.0, 6.0), ylim = c(-6.0, 6.0), zlim = c(-6.0, 6.0)
)

plot3D::text3D(x = W$PC1 + 0.05, 
               y = W$PC2 + 0.05, z = W$PC3 + 0.5,  
               labels = W$ID, add = TRUE, colkey = FALSE, cex = 1)

par(xpd=TRUE)}

dev.off()



#-------------------------------------------------------------------------------#
# Creación de gráfico compuesto

gmc <- ((gm1 + gm2) / (gm3 + gm4)) + plot_annotation(tag_levels = 'A')

ggsave('4_CompuestoPC.pdf', gmc, 'pdf', 'figures', 1, 5*2, 4*2)
ggsave('4_CompuestoPC.png', gmc, 'png', 'figures', 2, dpi = 'print')
