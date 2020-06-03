##------------------------------------------------------------------------#
## Nombre del Script: Revisión de trayectoria de algoritmos de evaluación 
##  covariables  ----------------------------------------------------
##  
## Propósito del Script: extraccion de datos sobre criterios de convergencia
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 19-05-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
# Selección de directorio de trabajo
setwd(file.path('C:', 'Users', 'Daniel', 'OneDrive', 'Documents', 
                '(Proyecto)_Estudio_PKPD', 'CEFEPIME', '14. Covariables_2'))


# Carga de paquetes
require(tidyverse)
require(ggplot2)
require(gridExtra)
require(gganimate)
require(patchwork)
require(ggrepel)

#-------------------------------------------------------------------------------#
# Carga de archivos de manejo

# Covariables en espacio muestral
samp_space <- c("Cl_tSCRMGDL", "Cl_logtWTKG", "Cl_logtIMCKGM2", "Cl_logtSCM2", 
                "Cl_logtRAN", "Cl_LLP", "Cl_LMP", "V1_logtAGEA", "V1_tSCRMGDL", 
                "V1_logtPROGDL", "V1_LMP", "Q_logtAGEA", "Q_logtIMCKGM2", 
                "Q_logtCLCRMLMIN", "Q_logtRAL", "Q_ANTU", "Q_LLP", "V2_SEXF", 
                "V2_logtAGEA", "V2_logtWTKG", "V2_tSCRMGDL", "V2_logtPROGDL", 
                "V2_logtALBGDL", "V2_logtRAL", "V2_ANTU", "V2_LLP", "V2_LMP")

#-------------------------------------------------------------------------------#
# Lectura de archivos de modelos
# Se abren los archivos y se filtran los pares especificados como parte del espacio 
# muestral

dfCOSSACBIC <- read_csv("2_Modelo_COSSAC_BICc/resumen_convergencia.csv") %>%
  filter(Par %in% samp_space)

dfCOSSACLRT <- read_csv("3_Modelo_COSSAC_LRT/resumen_convergencia.csv") %>%
  filter(Par %in% samp_space)

dfSCMBIC    <- read_csv("4_Modelo_SCM_BICc/resumen_convergencia.csv") %>%
  filter(Par %in% samp_space)

dfSCMLRT    <- read_csv("5_Modelo_SCM_LRT/resumen_convergencia.csv") %>%
  filter(Par %in% samp_space)

#-------------------------------------------------------------------------------#
# Gráficos de tilas --------------------------------------------------------------
#-------------------------------------------------------------------------------#
# Color asignado de acuerdo a criterio incluido (1) o no (0), para mostrar cada 
# tila en un color de acuerdo a la presencia o no del par covariable-parámetro 
# en la iteración.

Criterio <- c("1" = 'red', "0" = 'gray99') # Cuando no se incluye el par no se muestra

# Configurar el tema de gráficos sin leyendas
theme_set(theme_bw() + 
            theme(panel.border = element_rect(fill = NA, colour = 'black'),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  legend.position = "none" ))


# Gráfico 1 - Muestra las covariables incluidas en el modelo
#' Gráfico de tilas
#'
#' @param data tabla de datos con información de convergencia 
#'
#' @return objeto ggplot con gráfico de tilas
#'
#' @examples tile.plot(data_COSSAC_BIC)
#' 
tile.plot <- function(data) {
  data %>%
  # Limpiar los nombres de covariables en el set de datos, y creación de variable 
  # con los pares en cada iteración.
  mutate(Covariable = str_replace(Covariable, 'logt|^t', ''),
         Covariable = factor(Covariable),
         Par = paste0(Parametro, '_', Covariable),
         Criterio = factor(Criterio), 
         Par = fct_rev(Par)) %>% 
  # Gr{afico}
  ggplot(aes(x = Iteracion, y = Par, fill = Criterio)) + 
  geom_tile(col = 'gray50', size = 0.1) +
  xlab("Iteración") + ylab('Parámetro - Covariable') +
  scale_x_continuous(position = 'top') +
  scale_fill_manual(values = Criterio) +
  theme(axis.text.x = element_text(angle = 90),
        axis.text.y = element_text(size = 6, hjust = 0), 
        axis.title.y = element_blank())
}

#-------------------------------------------------------------------------------#
# Creación de gráficos y asignación a variables
# Se nota en este momento que los algoritmos generan resultados similares sin 
# importar el criterio de seguimiento de GIC. 

Gtile1a <- tile.plot(dfCOSSACBIC)
Gtile1b <- tile.plot(dfCOSSACLRT)
Gtile2a <- tile.plot(dfSCMBIC)
Gtile2a <- tile.plot(dfSCMLRT)

# Crear gráfico compuesto
Gtile.comp <- (Gtile1a + Gtile2a) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(widths = c(1,2))

# Almacenamiento
ggsave("Res_modelos_tilas.pdf", Gtile.comp, device = 'pdf', 
       width = 8, height = 4)

#-------------------------------------------------------------------------------#
# Gráfico con seguimiento de convergencia ---------------------------------------
#-------------------------------------------------------------------------------#

#' Función para resumir los datos de convergencia
#' Permite obtener una tabla que sólo muestre el criterio de convergencia obtenido 
#' en cada iteración. 
#' @param data tabla
#'
#' @return data.frame con GIC de cada convergencia
#'
#' @examples sum.fun(dfCOSSACBIC)
#' 
sum.fun <- function(data) {
  group_by(data, Iteracion) %>%
    summarise(LRT  = mean(LRT), BICc = mean(BICc))
}  

#-------------------------------------------------------------------------------#
#' Función para obtener un gráfico de seguimiento de convergencia
#'
#' @param data  tabla de datos
#' @param var_y variable de criterio de información 
#' @param crit  criterio que indica superiorirdad de modelo
#' @param colour color de tema
#'
#' @return objeto gráfico de tipo dispersión
#' 
#' @examples sum.fun(dfCOSSACBIC)
#' 
plot.fun <- function(data, var_y = 'LRT', crit=420, colour = 'black', ylab='OFV') {
  var_y_quo <- rlang::ensym(var_y)
  
  G1 <- data %>% 
    ggplot(aes(x = Iteracion, y = !!var_y_quo)) + 
    geom_point(colour = ifelse(data[var_y] <= crit, 'red', colour)) +
    geom_text_repel(aes(label = ifelse(data[var_y] <= crit, Iteracion, "")),
                    box.padding = 0.5, max.overlaps = Inf) + 
    geom_line(col = colour) + 
    xlab('Iteración') + ylab(ylab)
  
  return(G1)
}

#-------------------------------------------------------------------------------#
# Creación de gráficos compuestos 

set.seed(123465)

# En la selección por el algoritmo COSSAC, la iteraciòn con el menor valor de 
# LRT es N. 24 con 474; y BICc es N. 19 con 531.
# En la selecciòn por el algoritmo SCM, la iteraciòn con el menor valor de LRT 
# es N. 96 con 399; y BICc es N. 96 con 446
# slice(dfSCMBIC, which.min(BICc))

G1_comp <- 
  sum.fun(dfCOSSACBIC) %>% plot.fun(crit = 474+3.84, colour="#0000EE") + 
  sum.fun(dfSCMBIC) %>% plot.fun(crit = 400+3.84) +
  plot_layout(widths = c(1,2))

G2_comp <- 
  sum.fun(dfCOSSACBIC) %>% plot.fun(var_y = 'BICc', crit = 530.9+3.84, colour="#0000EE", ylab = 'BICc') + 
  sum.fun(dfSCMBIC) %>% plot.fun(var_y = 'BICc', crit = 449.0+3.84, ylab = 'BICc')  +
  plot_layout(widths = c(1,2))


G_tot_comp <- 
  (G1_comp & coord_cartesian(ylim = c(400,570))) / 
  (G2_comp & coord_cartesian(ylim = c(420,630))) +
  plot_annotation(tag_levels = 'A')

# Almacenamiento
ggsave("Res_Converg_modelos.pdf", G_tot_comp, device = 'pdf', 
       width = 8, height = 6)

 


filtrar.data <- function(data, var_orden) {
    var_orden_quo = rlang::ensym(var_orden)
      
    df1 <- data %>%
    select(-Parametro, -Covariable, -Par1) %>% 
    pivot_wider(names_from = Par, values_from = Criterio) %>% 
    rowwise() %>% 
    mutate(N = sum(c_across(matches("\\_")))) %>% 
    ungroup() %>% 
    arrange(!!var_orden_quo) %>%
    select(Iteracion, LRT, BICc, N)
  
  df2 <- gt::gt(df1)
  
  return(list(df1, df2))
}

# df1 <- dfSCMBIC %>%
#   select(-Parametro, -Covariable, -Par1) %>% 
#   pivot_wider(names_from = Par, values_from = Criterio) %>% 
#   select(-Iteracion, -LRT, -BICc) %>% 
#   as.matrix() 
# 
# dim(df1)
# 
# df2 <- matrix(nrow = dim(df1)[1], ncol = dim(df1)[2])
# df3 <- vector(length = dim(df1)[1])
# 
# for (i in 1:dim(df1)[1]) {
#   for (j in 1:dim(df1)[2]) {
#     df2[i,j] <- paste(colnames(df1)[j], df1[i,j])
#   }
# }
# 
# for (i in 1:dim(df1)[1]) {
#   df3[i] <- paste(df2)
# }
# 
# apply(., c(1,2), function(x) paste(colnames(x),x,sep='-'))


filtrar.data(dfSCMBIC, 'LRT')
filtrar.data(dfSCMBIC, 'BICc')

filtrar.data(dfCOSSACBIC, 'LRT')[[1]]x
filtrar.data(dfCOSSACBIC, 'BICc')[[1]]









