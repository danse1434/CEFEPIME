##------------------------------------------------------------------------#
## Nombre del Script: Script de generación de VPC  ------------------------
##  
## Proposito del Script: crear gráficos de chequeo predictivo visual mediante  
##  simulación por medio del paquete MLXR de Monolix. 
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion: 18-02-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
##########################################################################-
# Introducción -----------------------------------------------------
##########################################################################-
# Carga de paquetes
require(tidyverse)
require(rlang)
require(mlxR)

##########################################################################-
# Selección de directorio principal
setwd(file.path('F:','Documentos','(Proyecto)_Estudio_PKPD','CEFEPIME',
                '04. Residual', 'Modelo_8', 'RES_M1'))

##########################################################################-
# Simulación de parámetros teóricos ---------------------------------------
##########################################################################-
# Especificación del modelo
model = inlineModel("
[LONGITUDINAL]
input = {Cl, V1, Q, V2, a, b}
EQUATION:
V = V1 
k = Cl/V1 
k12 = Q/V1 
k21 = Q/V2

Cc = pkmodel(V, k, k12, k21)

DEFINITION:
y = {distribution=normal, prediction=Cc, errorModel=combined2(a, b)}

;--------------------------------------------------------

[INDIVIDUAL]
input = {Cl_pop, omega_Cl, Q_pop, omega_Q, V1_pop, omega_V1, V2_pop, omega_V2}
DEFINITION:
Cl = {distribution=logNormal, typical=Cl_pop, sd=omega_Cl}
Q = {distribution=logNormal, typical=Q_pop, sd=omega_Q}
V1 = {distribution=logNormal, typical=V1_pop, sd=omega_V1}
V2 = {distribution=logNormal, typical=V2_pop, sd=omega_V2}
")

# Esquema de Administración
adm <- list(time = seq(0, 80, 8),
            amount = 2000,
            tinf = 15 / 60)
# Esquema de Observación
Cc1  <- list(name = 'y', time = seq(from = 72, to = 80, by = 0.01))
# Parámetros de modelo farmacocinético
param <- c(Cl_pop = 13.5, #L
           V1_pop = 21.8, #L/h
           Q_pop = 36.5, #L
           V2_pop = 16.4, #L/h
           omega_Cl = 0.258, 
           omega_V1 = 6.79E-8,
           omega_Q = 1.67,
           omega_V2 = 0.721,
           a = 1.55, b = 0.0572)
# Tamaño de Población
g1 <- list(size = 1E3,
           level = 'individual',
           parameter = param)
# Comando de simulación
res <- simulx(
  model = model,
  group = list(g1),
  output = list(Cc1),
  treatment = adm,
  settings = list(seed = 38918)
)

prctilemlx(res$y) +
  ylab('Concentración plasmática (mg/L)') + xlab('Tiempo (horas)')

##########################################################################-
# Apertura de archivo de datos --------------------------------------------
##########################################################################-
# Se debe abrir el archivo de datos con información farmacocinética en forma 
# de TAD para comparación con las simulaciones.
OBS_DATA <- read_delim("Monolix_data_TAD.csv", delim = ";",
                       escape_double = FALSE,
                       trim_ws = TRUE, na = ".") %>% 
  filter(!is.na(DV)) %>% 
  mutate(gr = ntile(TIME, 6))

##########################################################################-
# Intervalos de predicción ------------------------------------------------
##########################################################################-
PI <- 80 # Especificar IP en porcentaje. Es común el intervalo del 80% (10-90%)
CI <- 95 # Especificar el IC en porcentaje. Es común el 95%.

# Hacer un vector para percentil superior e inferior
perc_PI <- c(0+(1-PI/100)/2, 1-(1-PI/100)/2)
perc_CI <- c(0+(1-CI/100)/2, 1-(1-CI/100)/2)
           
##########################################################################-
# Especificar los tiempos de contenedores manualmente
bin_times <- c(72.00000, 72.29167, 72.54167, 73.25000, 74.50000, 76.25000, 80.0)
Number_of_bins <- length(bin_times)-1

res1 = res$y %>%
  as_tibble(.) %>%
  mutate(TAD = time - 72)

res2 = res1 %>% 
  mutate(gr = ntile(TAD, 100)) %>% 
  group_by(gr) %>% 
  summarise(TIME = mean(TAD), 
            ME = quantile(x = y, probs = 0.50), 
            LI = quantile(x = y, probs = 0.10), 
            LS = quantile(x = y, probs = 0.90),
            L95.1 = quantile(x = y, probs = 0.05),
            L95.2 = quantile(x = y, probs = 0.95),
            L99.1 = quantile(x = y, probs = 0.01),
            L99.2 = quantile(x = y, probs = 0.99))

vpc_colours <- c('P99' <- 'green1', 'P95' <- 'green1', 'P90' <- 'green1')
vpc_alphas <- c('P99' <- 1, 'P95' <- 0.2, 'P90' <- 0.05)

ribbon <- function(min, max, fialp) {
  min <- rlang::ensym(min)
  max <- rlang::ensym(max)
  list(geom_ribbon(mapping = aes(
    ymin = !!min, ymax = !!max, fill = fialp, alpha = fialp ))) %>%
    return(.)
}

ggplot(res2, aes(x = TIME)) +
  ribbon(LI, LS, 'P90') +
  ribbon(L95.1, L95.2, 'P95') +
  ribbon(L99.1, L99.2, 'P99') +
  geom_line(mapping = aes(y = ME), size = 1) +
  geom_vline(xintercept = 0, lty = 'dashed') +
  geom_point(
    OBS_DATA %>% filter(!is.na(DV)),
    mapping = aes(x = TIME, y = DV),
    inherit.aes = F, shape = 'o'
  ) + 
  theme_bw() +
  theme(panel.grid = element_line(colour = NULL, linetype = 'blank'),
        legend.position = c(0.8, 0.8)) +
  ylab('Concentración plasmática (mg/L)') + xlab('TAD') +
  scale_fill_manual(values = vpc_colours, name = 'Percentiles')+
  scale_alpha_manual(values = vpc_alphas, name = 'Percentiles')


##########################################################################-
# VPC de percentil simple -------------------------------------------------
##########################################################################-
OBS_DATA1 <- OBS_DATA %>%
  group_by(gr) %>%
  summarise(
    TAD = mean(TIME),
    ME = quantile(x = DV, probs = 0.50),
    LI = quantile(x = DV, probs = 0.10),
    LS = quantile(x = DV, probs = 0.90)
  )


vpc_colour <- c('P_obs' <- 'red', 'P_sim' <- 'black')
vpc_lty <- c('P_obs' <- 'dashed', 'P_sim' <- 'solid')


lines <- function(data, x, y, ltcol, type) {
  x = rlang::ensym(x)
  y = rlang::ensym(y)
  if(type == 1){
    return(
      list((geom_line(mapping = aes(y = !!y, colour = ltcol, lty = ltcol)))) 
    )
  } else {
    return(
      list(geom_line(data = data, mapping = aes(x = !!x, y = !!y, 
                                    colour = ltcol, lty = ltcol)),
            geom_point(data = data, mapping = aes(x = !!x, y = !!y, 
                                     colour = ltcol))) 
    )}
}

ggplot(res2, aes(x = TIME)) + 
  lines(y = LI, ltcol = 'P_sim', type = 1) +
  lines(y = LS, ltcol = 'P_sim', type = 1) +
  lines(y = ME, ltcol = 'P_sim', type = 1) +
  
  lines(OBS_DATA1, TAD, LI, 'P_obs', type = 2) +
  lines(OBS_DATA1, TAD, LS, 'P_obs', type = 2) +
  lines(OBS_DATA1, TAD, ME, 'P_obs', type = 2) +
  theme_bw() +
  theme(panel.grid = element_line(colour = NULL, linetype = 'blank'),
        legend.position = c(0.8, 0.8)) +
  ylab('Concentración plasmática (mg/L)') + xlab('TAD') +
  scale_color_manual(values = vpc_colour, name = 'Percentiles') + 
  scale_linetype_manual(values = vpc_lty, name = 'Percentiles')


##########################################################################-
# VPC de bandas de confianza ----------------------------------------------
##########################################################################-
y_1_obsVsPred <-
  read_csv("ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt") %>% 
  mutate(y_1 = round(y_1, 2))

OBS_DATA2 <- OBS_DATA %>%
  mutate(DV = round(DV, 2)) %>% 
  left_join(., y_1_obsVsPred %>% select(ID, y_1, popPred), 
            by = c('ID', 'DV' = 'y_1'))


OBS_DATA3 <- OBS_DATA2 %>% 
  group_by(gr) %>% 
  summarise(TAD = mean(TIME), 
            ME = quantile(x = popPred, probs = 0.50),
            LI = quantile(x = popPred, probs = 0.10),
            LS = quantile(x = popPred, probs = 0.90) )

OBS_DATA2 <- OBS_DATA2 %>% 
  left_join(., OBS_DATA3, by = 'gr')

OBS_DATA2 %>% 
  mutate(pcVPC = DV * ME/popPred) %>% 
  ggplot(.) +
  geom_point(aes(x = TIME, y = popPred))
  



res3 <- res1 %>% 
  mutate(gr = ntile(TAD, 100)) %>% 
  group_by(gr)

sim_PI <- NULL

for(i in 1:1000){
  
  sim_vpc_ci <- res3 %>%
    filter(id %in% i) %>% # Select an individual replicate
    group_by(gr) %>% # Calculate everything per bin
    summarize(C_median = median(y), 
              C_lower = quantile(y, perc_PI[1]), 
              C_upper = quantile(y, perc_PI[2])) %>% # Calculate prediction intervals
    mutate(replicate = i) # Include replicate number
  
  sim_PI <- rbind(sim_PI, sim_vpc_ci)
}

sim_PI

# Calculate confidence intervals around these prediction intervals calculated with each replicate
sim_CI <- sim_PI %>%
  group_by(gr) %>%
  summarize(
    C_median_CI_lwr = quantile(C_median, perc_CI[1]),
    C_median_CI_upr = quantile(C_median, perc_CI[2]),
    # Median
    
    C_low_lwr = quantile(C_lower, perc_CI[1]),
    C_low_upr = quantile(C_lower, perc_CI[2]),
    # Lower percentages
    
    C_up_lwr = quantile(C_upper, perc_CI[1]),
    C_up_upr = quantile(C_upper, perc_CI[2]) # High percentages
  )


res2 %>% 
  left_join(., sim_CI, by = 'gr') %>% 
  ggplot(mapping = aes(x = TIME)) +
  # geom_ribbon(aes(ymin = C_median_CI_lwr, ymax = C_median_CI_upr), alpha = 0.1) +
  geom_ribbon(aes(ymin = C_low_lwr, ymax = C_low_upr), alpha = 0.1) +
  geom_ribbon(aes(ymin = C_up_lwr, ymax = C_up_upr), alpha = 0.1)

##########################################################################-
# Adicionar el bin correspondiente
# for(i in Number_of_bins:1) {
#   res1$bin[res1$time <= bin_times[i + 1]] <- i
# }
# 
# 
# PRED_BIN <-
#   res1 %>%
#   group_by(bin) %>% # Calcular PRED por bin
#   summarize(predbin = median(y))
# 
# 
# prctilemlx(res$y) +
#   ylab('Concentración plasmática (mg/L)') + xlab('TAD') +
#   geom_point(
#     data = OBS_DATA %>% filter(!is.na(DV)),
#     mapping = aes(x = TIME + 72, y = DV),
#     inherit.aes = F
#   ) +
#   theme(panel.grid = element_line(colour = NULL, linetype = 'blank'))
# 
# 
# # [res1$replicate == 1, ] %>%
# #   group_by(BIN) %>% # Calculate the PRED per bin
# #   summarize(PREDBIN = median(PRED))
