require(tidyverse)
require(mlxR)
require(Rcpp)
require(RcppArmadillo)

MIC = c(1 * (2 ^ (seq(-10, 10, 0.05))))
Rcpp::sourceCpp('src/90_verificacion_PTA.cpp')
source('src/80_funciones.R')

vector_loc <- c('07', '10', '13', '16', '19', '40')

for (i in 1:6) {
  aux_dir <- file.path(getwd(), 'results', paste0('L0', vector_loc[[i]]))
  
  res <- readRDS(file.path(aux_dir, 'RES.rds'))
  
  ptm <- proc.time()
  # Segundo Objetivo PK-PD
  RES_PTA_2 <- res$ufCc %>% 
    group_by(group) %>% nest() %>% 
    mutate(mat = map(data, ~conv_matriz(.x))) %>% 
    mutate(PTA = map(mat, ~pta_verificador(.x, MIC, crit = 1.00, treshold = 1))) %>% 
    mutate(
      fTmasMIC = map(PTA, ~ magrittr::use_series(.x, 'fTmasMIC')),
      PTA      = map(PTA, ~ magrittr::use_series(.x, 'PTA'))
    )
  # 
  print(proc.time() - ptm)
  
  saveRDS(RES_PTA_2, file.path(aux_dir, 'RES_PTA_100tmasMIC.rds'))
  
}
