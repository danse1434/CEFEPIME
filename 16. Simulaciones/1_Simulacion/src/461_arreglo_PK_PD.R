require(tidyverse)
require(mlxR)
require(Rcpp)
require(RcppArmadillo)

MIC = c(1 * (2 ^ (seq(-10, 10, 0.05))))
Rcpp::sourceCpp('src/90_verificacion_PTA.cpp')
source('src/80_funciones.R')

id_SRC <- c('L401','L402')

for (i in seq_along(id_SRC)) {
aux_dir <- file.path(getwd(), 'results', id_SRC[i])

res <- readRDS(file.path(aux_dir, 'RES.rds'))

ptm <- proc.time()
# Primer Objetivo PK-PD
RES_PTA_1 <- res$ufCc %>% 
  group_by(group) %>% nest() %>% 
  mutate(mat = map(data, ~conv_matriz(.x))) %>% 
  mutate(PTA = map(mat, ~pta_verificador(.x, MIC, crit = 0.5, treshold = 4))) %>% 
  mutate(
    fTmasMIC = map(PTA, ~ magrittr::use_series(.x, 'fTmasMIC')),
    PTA      = map(PTA, ~ magrittr::use_series(.x, 'PTA'))
  )
# 
print(proc.time() - ptm)

ptm <- proc.time()
# Segundo Objetivo PK-PD
RES_PTA_2 <- res$ufCc %>% 
  group_by(group) %>% nest() %>% 
  mutate(mat = map(data, ~conv_matriz(.x))) %>% 
  mutate(PTA = map(mat, ~pta_verificador(.x, MIC, crit = 0.95, treshold = 1))) %>% 
  mutate(
    fTmasMIC = map(PTA, ~ magrittr::use_series(.x, 'fTmasMIC')),
    PTA      = map(PTA, ~ magrittr::use_series(.x, 'PTA'))
  )
# 
print(proc.time() - ptm)

ptm <- proc.time()
# Segundo Objetivo PK-PD
RES_PTA_3 <- res$ufCc %>% 
  group_by(group) %>% nest() %>% 
  mutate(mat = map(data, ~conv_matriz(.x))) %>% 
  mutate(PTA = map(mat, ~pta_verificador(.x, MIC, crit = 0.6, treshold = 1))) %>% 
  mutate(
    fTmasMIC = map(PTA, ~ magrittr::use_series(.x, 'fTmasMIC')),
    PTA      = map(PTA, ~ magrittr::use_series(.x, 'PTA'))
  )
# 
print(proc.time() - ptm)

saveRDS(RES_PTA_1, file.path(aux_dir, 'RES_PTA_50tmas4MIC.rds'))
saveRDS(RES_PTA_2, file.path(aux_dir, 'RES_PTA_100tmasMIC.rds'))
saveRDS(RES_PTA_3, file.path(aux_dir, 'RES_PTA_60tmasMIC.rds'))
  
}
