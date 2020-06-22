##------------------------------------------------------------------------------#
## Nombre del Script: Simulación MC con modelo final ----------------------------
##  
## Propósito del Script: realizar una descripción del comportamiento de los 
## individuos en cuanto a PTA vs MIC. 
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 16-06-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
# Carga de paquetes
require(scales)
require(rlang)
require(tidyverse)
require(Rcpp)
require(RcppArmadillo)
require(mlxR)

# Carga de archivo Rcpp
Rcpp::sourceCpp('src/90_verificacion_PTA.cpp')
# Carga de archivo con funciones R
source('src/80_funciones.R', encoding = 'UTF-8')

# Carga de modelo final
model_ori <- 'data/Modelo_Final.mlxtran'

#-------------------------------------------------------------------------------#
# Definición de Parámetros Modelo Final  ----------------------------------------
#-------------------------------------------------------------------------------#
# Efectos fijos
p1 <- c("Cl_pop" = 20.6, "beta_Cl_tSCRMGDL" = -0.415, "V1_pop" = 23.8, 
        "V2_pop" = 13.3, "Q_pop" = 23.4)
# Efectos aleatorios
p2 <- c("omega_Cl" = 0.224, "omega_V1" = 0.299, "omega_Q" = 0.946, 
        "omega_V2" = 0.795, "a" = 1.86)
# Covariable
p3 <- c("SCRMGDL" = 0.54)
# Fracción libre
p4 <- c("f" = 0.2)

#-------------------------------------------------------------------------------#
# Esquemas de tratamiento
# Colocar los esquemas de tratamiento
# CEP 2000mg q8h inf. 30min (DD = 6g)
adm1 <- list(time = seq(0, 8 * 10, by = 8), amount = rep(2000, 11), 
            tinf = rep(0.5, 11)) 
# CEP 2000mg q8h inf. 2hr (DD = 6g)
adm2 <- list(time = seq(0, 8 * 10, by = 8), amount = rep(2000, 11), 
             tinf = rep(2, 11))
# CEP 2000mg q8h inf. 4hr (DD = 6g)
adm3 <- list(time = seq(0, 8 * 10, by = 8), amount = rep(2000, 11), 
             tinf = rep(4, 11))
# CEP 6000mg q24h inf. 24hr (DD = 6g)
adm4 <- list(time = seq(0, 24 * 3, by = 24), amount = rep(6000, 4), 
             tinf = rep(24, 4))


# Grupo de simulación
grp1 <- list(size = 5000, level = 'individual', treatment = adm1, 
             parameter = list(p1, p2, p3, p4))
grp2 <- list(size = 5000, level = 'individual', treatment = adm2, 
             parameter = list(p1, p2, p3, p4))
grp3 <- list(size = 5000, level = 'individual', treatment = adm3, 
             parameter = list(p1, p2, p3, p4))
grp4 <- list(size = 5000, level = 'individual', treatment = adm4, 
             parameter = list(p1, p2, p3, p4))


# out0 <- list(name = 'Cc', time = seq(0, 8, length.out = 100))
out1 <- list(name = 'ufCc', time = seq(0, 24, length.out = 100))

# Simulación
ptm <- proc.time()

res <- simulx(model  = model_ori,
              output = list(out1),
              group  = list(grp1, grp2, grp3, grp4))

print(proc.time() - ptm)


# Configuración de tema
theme_set( theme_bw() )

# Gráfico con perfiles
# ggplot(res$ufCc, aes(x = time, y = ufCc, group = id, col = group)) +
#   geom_line() +
#   # geom_line(data=res$ufCc, aes(x=time, y = ufCc, col = id)) +
#   # facet_wrap(.~id, ncol = 3) +
#   scale_color_viridis_d() +
#   theme(legend.position = "none")

#-------------------------------------------------------------------------------#
# Perfiles representativos -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Objeto con perfiles representativos, se crean IP
PRCT1 <- mlxR::prctilemlx(res$ufCc, plot = FALSE)

# Objeto con cálculo de parámetros de exposición
# EXPO1 <- exposure(
#         model  = model_ori,
#         output = out1,
#         group  = list(grp1, grp2, grp3, grp4)
#         )

EXPO1 <- res$ufCc %>% 
  group_by(group) %>% 
  nest() %>% 
  mutate(mat = map(data, ~conv_tibble(.x)), 
         expo = map(mat, ~UDF_exposure(.x))) %>% 
  mutate(tmax_cmax = map(expo, ~magrittr::use_series(.x, 'Pmaximos')),
         expo = map(expo, ~magrittr::use_series(.x, 'AUC_estm')))
  

#-------------------------------------------------------------------------------#
# Cumplimiento PTA -----------------------------------------------------
#-------------------------------------------------------------------------------#
MIC = c(1*(2^(seq(-10,10,1))))

a <- res$ufCc %>% 
  group_by(group) %>% 
  nest() %>% 
  mutate(mat = map(data, ~conv_matriz(.x))) %>% 
  mutate(PTA = map(mat, ~pta_verificador(.x, MIC, crit = 0.5))) %>% 
  mutate(
    fTmasMIC = map(PTA, ~ magrittr::use_series(.x, 'fTmasMIC')),
    PTA = map(PTA, ~ magrittr::use_series(.x, 'PTA'))
  )

# a %>% 
#   select(group, PTA) %>% 
#   unnest(PTA) %>% 
#   ggplot(aes(x = MIC, y = PTA, col = group)) +
#   geom_line() + 
#   scale_x_continuous(
#     trans = 'log2',
#     breaks = trans_breaks("log2", function(x) 2^x),
#     labels = trans_format("log2", math_format(2^.x)))
# 
# a %>% 
#   select(group, fTmasMIC) %>% 
#   unnest(fTmasMIC) %>% 
#   pivot_longer(cols = matches('^ID'),names_to = 'ID', values_to = 'ind_pkpd') %>% 
#   ggplot(aes(x = MIC, y = ind_pkpd, col = ID)) +
#   geom_line() + 
#   scale_x_continuous(
#     trans = 'log2',
#     breaks = trans_breaks("log2", function(x) 2^x),
#     labels = trans_format("log2", math_format(2^.x))) +
#   facet_wrap(.~group, ncol = 2) +
#   theme(legend.position = "none")

# microbenchmark::microbenchmark(pta_verificador(df1, c(1 * (2^(seq(-10, 10, 0.1)))), 0.5))
# Se verifica que la complejidad del algoritmo es O(n^2), este corre en tiempo polinómico
