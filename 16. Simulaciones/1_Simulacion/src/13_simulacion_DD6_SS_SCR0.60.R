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
p3 <- c("SCRMGDL" = 0.60)
# Fracción libre
p4 <- c("f" = 0.2)

#-------------------------------------------------------------------------------#
# Esquemas de tratamiento
# Colocar los esquemas de tratamiento
# CEP 2000mg q8h inf. 30min (DD 6g)
adm1 <- list(time   = seq(0, 8 * 21, by = 8), 
             amount = rep(2000, 22), 
             tinf    = rep(0.5, 22)) 
# CEP 2000mg q8h inf. 2hr (DD 6g)
adm2 <- list(time   = seq(0, 8 * 21, by = 8), 
             amount = rep(2000, 22), 
             tinf   = rep(2, 22))
# CEP 2000mg q8h inf. 4hr (DD 6g)
adm3 <- list(time   = seq(0, 8 * 21, by = 8), 
             amount = rep(2000, 22), 
             tinf   = rep(4, 22))
# CEP 6000mg q24h inf. 24hr (DD 6g)
adm4 <- list(time   = c(0.0, seq(0 + 0.5, 24 * 8, by = 24)), 
             amount = c(2000, rep(6000, 8)), 
             tinf   = c(0.5, rep(24, 8)))

# Grupo de simulación
grp1 <- list(size = 2500, level = 'individual', treatment = adm1, 
             parameter = list(p1, p2, p3, p4))
grp2 <- list(size = 2500, level = 'individual', treatment = adm2, 
             parameter = list(p1, p2, p3, p4))
grp3 <- list(size = 2500, level = 'individual', treatment = adm3, 
             parameter = list(p1, p2, p3, p4))
grp4 <- list(size = 2500, level = 'individual', treatment = adm4, 
             parameter = list(p1, p2, p3, p4))

# Definición de outputs
# out0 <- list(name = 'Cc', time = seq(0, 8, length.out = 100))
out1 <- list(name = 'ufCc', time = seq(120, 120+24, length.out = 250))

# Simulación
ptm <- proc.time()
res <- simulx(model  = model_ori,
              output = list(out1),
              group  = list(grp1, grp2, grp3, grp4))

print(proc.time() - ptm)

# Configuración de tema
theme_set( theme_bw() )

#-------------------------------------------------------------------------------#
# Perfiles representativos -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Objeto con perfiles representativos, se crean IP
# 
PRCT1 <- mlxR::prctilemlx(res$ufCc, plot = FALSE)

#-------------------------------------------------------------------------------#
# Parámetros de exposición de perfiles ------------------------------------------
#-------------------------------------------------------------------------------#
# Agrupar los datos y convertir en columna lista
# Convertir al data.frame en la columna lista en un archivo tibble con una 
# especificación dada.
# Calcular los parámetros de exposición mediante el código en C++
# Extraer las tabla de Pmaximos y AUC_estm en sus respectivas lista-columnas
# 
EXPO1 <- res$ufCc %>% 
  group_by(group) %>% nest() %>% 
  mutate(mat = map(data, ~conv_tibble(.x)), 
         expo = map(mat, ~UDF_exposure(.x))) %>% 
  mutate(
    tmax_cmax = map(expo, ~ magrittr::use_series(.x, 'Pmaximos')),
    expo      = map(expo, ~ magrittr::use_series(.x, 'AUC_estm'))
  )


#-------------------------------------------------------------------------------#
# Cumplimiento PTA -----------------------------------------------------
#-------------------------------------------------------------------------------#
MIC = c(1 * (2 ^ (seq(-10, 10, 0.05))))

RES_PTA <- res$ufCc %>% 
  group_by(group) %>% nest() %>% 
  mutate(mat = map(data, ~conv_matriz(.x))) %>% 
  mutate(PTA = map(mat, ~pta_verificador(.x, MIC, crit = 0.5, treshold = 4))) %>% 
  mutate(
    fTmasMIC = map(PTA, ~ magrittr::use_series(.x, 'fTmasMIC')),
    PTA      = map(PTA, ~ magrittr::use_series(.x, 'PTA'))
  )

#-------------------------------------------------------------------------------#
# Almacenamiento Carpetas -------------------------------------------------------
#-------------------------------------------------------------------------------#
CODIGO_CARP <- 'L013'
aux_dir <- file.path(getwd(), 'results', CODIGO_CARP)


dir.create(aux_dir)
saveRDS(res, file.path(aux_dir, 'RES.rds'))
saveRDS(PRCT1, file.path(aux_dir, 'PRCT1.rds'))
saveRDS(EXPO1, file.path(aux_dir, 'EXPO1.rds'))
saveRDS(RES_PTA, file.path(aux_dir, 'RES_PTA.rds'))

# rm(list = ls())

