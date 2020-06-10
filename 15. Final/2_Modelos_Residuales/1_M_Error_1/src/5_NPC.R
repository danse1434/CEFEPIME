##------------------------------------------------------------------------------#
## Nombre del Script: Cálculo de NPC modelo de covariables final 
##  
## Propósito del Script:  Realizar una chequeo predictivo numérico en el modelo 
## de covariables final
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 30-05-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
require(mlxR)
require(tidyverse)
require(patchwork)

#-------------------------------------------------------------------------------#
# Introducción  -----------------------------------------------------------------
#-------------------------------------------------------------------------------#
# 
project.file <- 'M_Error.mlxtran'

# Remuestreo de pacientes iniciales

# Números de simulaciones del diseño experimental
N <- 500*15 # 7500

# Realizar simulaciones desde el archivo
res1 <- simulx(project = project.file, 
               group = list(size = N))

# Crear un data.frame con simulaciones
res2 <- res1$y_1 %>% 
  as_tibble() %>% 
  left_join(res1$originalId, by = c('id' = 'newId')) %>% 
  add_column(sim_id = rep(1:(N/15), each=15*6)) %>% 
  group_by(oriId) %>% 
  mutate(TAD = time - min(time))

res2

# ggplot(res2, aes(x=TAD, y=y_1, col=oriId)) +
#   geom_point() + 
#   facet_wrap(.~id)

# Leer archivo TAD con observaciones
data_TAD <- read_csv("data/2_data_TAD.csv", na = '.')
data_TAD_1 <- data_TAD %>% filter(EVID == 0)

#-------------------------------------------------------------------------------#
# Chequeo Predictivo Numérico (NPC)   -------------------------------------------
#-------------------------------------------------------------------------------#
#' Función de procedimiento NPC
#'
#' @param obs vector con concentraciones observados (Cobs) 
#' @param sim vector con concentraciones simuladas (Csim)
#' @param PI  vector de intervalos de predicción (IP)
#'
#' @return
#' Los vectores Cobs y Csim tienen que estar alineados de manera que se repite 
#' la simulación varios veces. 
#' @export 
#' data.frame con ratio de OE, intervalo inferior, y superior, y tipo de límite 
#' de intervalo de predicción (IP).
#'
#' @examples npc(obs = data_TAD_1$DV, sim = res2$y_1, PI = seq(5, 90, by = 5))
#' 
npc <- function (obs, sim,PI){
  nobs <- length(obs)
  nrep <- length(sim)[1]/nobs
  # Cálculo de límites inferior y superior de cada intervalo de predicción (IP) 
  # en forma de percentiles
  percentiles <- c(50-PI/2, (50+PI/2))
  # Conversión de vector Csim en matriz con una columna por cada repetición
  matsimfull <- matrix(sim, ncol = nrep)
  # Cálculo de los percentiles en cada simulación
  matsimper <- t(apply(matsimfull,1,quantile,percentiles/100))
  matsimper_lower <- matsimper[,1:length(PI)]
  matsimper_upper <- matsimper[,(length(PI)+1):length(percentiles)]
  # Verificación de cobertura de datos simulados
  sim_outlower <- c()
  sim_outupper <- c()
  
  for (i in 1:length(PI)) {
    # Verificación si datos simulados son menores a percentiles
    matsim_outlower <- ifelse(matsimfull < matsimper_lower[, i], 1, 0)
    
    sim_outlower <-
      cbind(sim_outlower, quantile(apply(matsim_outlower, 2, mean), c(0.025, 0.5, 0.975)))
    
    # Verificación si datos simulados son mayores a percentiles
    matsim_outupper <-
      ifelse(matsimfull > matsimper_upper[, i], 1, 0)
    
    sim_outupper <-
      cbind(sim_outupper, quantile(apply(matsim_outupper, 2, mean), c(0.025, 0.5, 0.975)))
  }
  
  # Repetición de vector de observaciones con columnas para cada intervalo de predicción
  yobs_lower      <- matrix(rep(obs, length(PI)), ncol = length(PI))
  yobs_upper      <- matrix(rep(obs, length(PI)), ncol = length(PI))
  # Verificación si observación es menor a limite inferior para cada IP
  yobs_lower      <- ifelse(yobs_lower < matsimper_lower, 1, 0)
  # Verificación si observación es menor a limite superior para cada IP
  yobs_upper      <- ifelse(yobs_upper > matsimper_upper, 1, 0)
  # Cálculo de proporción de OUTLIERS como valor medio
  obsper_outlower <- apply(yobs_lower, 2, mean)
  obsper_outupper <- apply(yobs_upper, 2, mean)
  # Cálculo de cobertura como razón OE - proporción de valores observados (Cobs) 
  # sobre valores esperados (Csim)  
  lower_cover     <- obsper_outlower / sim_outlower[2,]
  CI_lower_cover  <- t(t(sim_outlower[c(1, 3),]) / sim_outlower[2,])
  upper_cover     <- obsper_outupper / sim_outupper[2,]
  CI_upper_cover  <- t(t(sim_outupper[c(1, 3),]) / sim_outupper[2,])
  # Número de datos OUTLIERS 
  lower_out <-
    (lower_cover < CI_lower_cover[1, ]) + (lower_cover > CI_lower_cover[2, ])
  
  upper_out <-
    (upper_cover < CI_upper_cover[1, ]) + (upper_cover > CI_upper_cover[2, ])
  
  # Creación de datos de cobertura
  cover <- as.data.frame(cbind(
    rep(PI, 2),
    c(lower_cover, upper_cover),
    c(lower_out, upper_out),
    rbind(t(CI_lower_cover), t(CI_upper_cover))
  ), row.names = F)
  
  cover$Type <-
    rep(c("Limite IP inferior", "Limite IP superior"), each = length(PI))
  colnames(cover) <- c("PI","Ratio","Outliers","Lwr","Sup","Tipo")
  return(cover)}

#-------------------------------------------------------------------------------#
# Tabla con resultados de NPC

df <- npc(obs = data_TAD_1$DV, sim = res2$y_1, 
          PI = seq(5, 90, by = 5))

# Gráfico de NPC

NPC_plot <- df %>% 
  ggplot(aes(x = PI, y = Ratio, 
             group = Tipo, col = Tipo, fill = Tipo)) + 
  geom_ribbon(aes(ymin = Lwr, ymax = Sup), alpha = 0.8) +
  geom_line(col = 'black') + 
  geom_point(col = 'black') + 
  facet_grid(Tipo ~ .) + 
  geom_hline(yintercept = 1, lty = 'dotted') +
  ylab('Ratio O/E') + xlab('IP') + 
  theme_bw() +
  theme(panel.border = element_rect(fill = NULL, colour = 'black'),
        legend.position = "none") +
  scale_fill_manual(values = c('#85B8FF', '#FF886B')) +
  scale_color_manual(values = c('#85B8FF', '#FF886B'))

NPC_plot

# Apertura de gráfico de VPC almacenado 
VPC_plot <- readRDS('figures/RDS/pcVPC_percentil.rds')

# Formación de gráfico compuesto
NPC_conplot <- (VPC_plot + theme_bw()) + 
                NPC_plot + plot_annotation(tag_levels = 'A')

# Almacenamiento del archivo PDF
ggsave(filename = 'figures/NPC_plot.pdf', NPC_plot, device = 'pdf', 
       width = 4/1.2, height = 6/1.2)

ggsave(filename = 'figures/VPC_NPC_plot.pdf', NPC_conplot, device = 'pdf', 
       width = 8, height = 5)

# Almacenamiento como RDS
saveRDS(NPC_plot, 'figures/RDS/NPC_percentil.RDS')

