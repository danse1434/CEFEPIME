##------------------------------------------------------------------------------#
## Nombre del Script: Evaluación de consistencia asintótica de estimadores -------
##
## Propósito del Script: Este procedimiento se realiza para encontrar el número 
## mínimo de configuraciones para estimar de forma apropiada a estadísticos de 
## exposición AUC, y Cmax.
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  14-03-2021
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#

require(patchwork)
require(Rcpp)
require(RcppArmadillo)
require(progress)

mdir <- file.path('data', 'modeloFinal.txt')
# rdir <- file.path('model', 'reference_model.mlxtran')

# Carga de archivo Rcpp
Rcpp::sourceCpp('src/90_verificacion_PTA.cpp')
# Carga de archivo con funciones R
source(file.path('src', '80_funciones.R'), encoding = 'UTF-8')
source(file.path('src', '690_paquetesSimulacion.R'), encoding = 'UTF-8')

# Vector de parámetros
p_DF <- read_csv(file.path('data', 'populationParameters.txt'))
p <- setNames(p_DF$value, p_DF$parameter)
p['SCRMGDL'] = 0.54
p['f'] = 0.2

#-------------------------------------------------------------------------------#
# 1. Determinación de N.Dosis para Estado Estacionario -------------------
#-------------------------------------------------------------------------------#
# Administración

nDosis <- 1:20
nDosis_ls1 <- vector(mode = 'list', length = length(nDosis))
nDosis_ls2 <- vector(mode = 'list', length = length(nDosis))
N = 5e2

MIC_vector <- c(1 * (2 ^ (seq(-9, 9, 1))))
#  

pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: (:eta)",
                       total = length(nDosis))

for (i in seq_along(nDosis)) {
  
  pb$tick()
  
  if (i <= 2) {
    warning("Se pasa el elemento")
  } else {
    adm <- list(time = seq(0, 8 * i, by = 8), amount = 2000, tinf = 0.5)
    
    out <- list(name = 'y_1', 
                time = seq(tail(adm$time, 4)[1], tail(adm$time, 4)[4], 
                           length.out = 50))
    
    # print(out)
    group <- list(size = N, level = 'individual',
                  treatment = adm, parameter = p)
  
    res <- simulx(model  = mdir, output = out,
                  group  = group)
  
    expo1 <- res$y_1 %>%
      as_tibble() %>%
      conv_tibble(values = y_1) %>%
      UDF_exposure()
    
    expo2 <- res$y_1 %>%
      conv_matriz(id_col = time,
                  names = id,
                  values = y_1) %>% 
      pta_verificador(crit = 0.5, treshold = 1, MIC = MIC_vector)
  
    
    nDosis_ls1[[i]] <- expo1
    nDosis_ls2[[i]] <- expo2
  }
}

nDosis_ls1[sapply(nDosis_ls1, is.null)] <- NULL
nDosis_ls2[sapply(nDosis_ls2, is.null)] <- NULL

nDosis_DF1 <- nDosis_ls1 %>% 
  lapply(function(x) left_join(x$AUC_estm, x$Pmaximos, by = 'ID')) %>% 
  map_df( ~ .x, .id = 'grupo') %>%
  tibble() %>%
  mutate(nDosis = nDosis[as.double(grupo) + 2],
         nDosis = factor(nDosis))

# nDosis_DF2 <- nDosis_ls2 %>% 
#   lapply(function(x) x$PTA) %>% 
#   map_df( ~ .x, .id = 'grupo') %>% 
#   mutate(nDosis = nDosis[as.double(grupo) + 2],
#          nDosis = factor(nDosis)) %>% 
#   as_tibble()

nDosis_DF2 <- nDosis_ls2 %>% 
  lapply(function(x) x$fTmasMIC) %>% 
  map_df( ~ .x, .id = 'grupo') %>% 
  mutate(nDosis = nDosis[as.double(grupo) + 2],
         nDosis = factor(nDosis)) %>% 
  as_tibble() %>% 
  pivot_longer(cols = matches('^ID\\.')) %>% 
  pivot_wider(names_from = MIC, values_from = value)


g1 <- nDosis_DF1 %>% 
  ggplot() + 
  geom_boxplot(aes(x = nDosis, y = AUC)) + 
  xlab('N. Dosis') + ylab('AUC (mg*h/L)')

g2 <- nDosis_DF %>% 
  ggplot() + 
  geom_boxplot(aes(x = nDosis, y = Cmax)) + 
  xlab('N. Dosis') + ylab(expression(C[max]~(mg/L)))

g3 <- nDosis_DF2 %>% 
  ggplot() +
  geom_boxplot(aes(x = nDosis, y = `4`)) + 
  xlab('N. Dosis') + ylab('fT<MIC (4mg/L)')

g4 <- nDosis_DF2 %>% 
  ggplot() +
  geom_boxplot(aes(x = nDosis, y = `8`)) + 
  xlab('N. Dosis') + ylab('fT<MIC (8 mg/L)')


gT <- g1 + g2 + g3 + g4 + 
  plot_layout(ncol = 2) + 
  plot_annotation(
    title = 'N. de dosis para estado estacionario',
    subtitle = 'Simulaciones Cefepime 2g q8h tinf 30min')

ggsave('001_Convergencia_Ndosis.pdf', gT, 'pdf', 'figures', 1, 
       8, 6)

# Con 4 dosis se puede encontrar el estado estacionario para 1g q12h de VAN
# 
#-------------------------------------------------------------------------------#
# 2. N. de individuos para estimar estadístico ----------------
#-------------------------------------------------------------------------------#
N = 1e4

MIC_vector <- c(1 * (2 ^ (seq(1, 3, 1))))

adm <- list(time = seq(0, 8 * 6, by = 8), amount = 2000, tinf = 0.5)

out <- list(name = 'y_1', 
            time = seq(tail(adm$time, 4)[1], tail(adm$time, 4)[4], 
                       length.out = 50))

group <- list(size = N, level = 'individual',
              treatment = adm, parameter = p)

res <- simulx(model  = mdir, output = out,
              group  = group)

# expo1 <- res$y_1 %>%
#   as_tibble() %>%
#   conv_tibble(values = y_1) %>%
#   UDF_exposure()

expo2 <- res$y_1 %>%
  conv_matriz(id_col = time,
              names = id,
              values = y_1) %>% 
  pta_verificador(crit = 0.5, treshold = 1, MIC = MIC_vector)

dfExpo2 <- expo2$fTmasMIC %>% 
  as_tibble() %>% 
  pivot_longer(cols = matches('^ID\\.')) %>% 
  pivot_wider(names_from = MIC, values_from = value)

extraerMuestreos <- function(df, variable, N, repeticiones = 10) {
  
  repDF <- tibble(.rows = N) 
  
  for (i in 1:repeticiones) {
    # print(paste0('rep', i))
    # print(sample_n(df, N) %>% pull(variable))
    
    repDF[[paste0('rep', i)]] <- pull(slice_sample(df, n = N, replace = TRUE), variable)
  }
  
  return(repDF)
}


truePTA <- dfExpo2$`8` %>% mean()

exploracion1 <-
  tibble(N = c(25, 50, 100, 200, 300, 500, 1000, 1500, 2000, 3000, 4000, 5000, 7500, 
               10000)) %>% 
  mutate(M = map(N, ~extraerMuestreos(dfExpo2, '8', .x, repeticiones = 30)),
         S = map(M, ~summarise(.x, across(everything(), ~mean(.x))))) %>% 
  select(-M) %>% 
  unnest(S) %>% 
  pivot_longer(cols = matches('^rep'))

gexploracion11 <- exploracion1 %>% 
  ggplot(aes(x = N, y = value)) + 
  geom_point() + 
  ylab('PTA (fT>MIC; MIC = 8)') +
  geom_vline(xintercept = 1e3, lty = 'dashed', col = 'blue2') + 
  geom_hline(yintercept = truePTA, lty = 'dashed', col = 'black') + 
  annotate(geom = 'text', x = 1.4e3, y = 0.56, 
           label = 'N = 1000', angle = 90)

gexploracion12 <- exploracion1 %>% 
  ggplot(aes(x = factor(N), y = value)) + 
  geom_boxplot() + 
  scale_x_discrete(guide = guide_axis(n.dodge = 2)) + 
  ylab('PTA (fT>MIC; MIC = 8)') + 
  xlab('N')


exploracion2 <- exploracion1 %>%
  mutate(e = value - truePTA) %>%
  group_by(N) %>%
  summarise(MAE = mean(abs(e)),
            RMSE = sqrt(mean(e^2)))

gexploracion13 <- exploracion2 %>% 
  ggplot(aes(x = N, y = MAE)) +
  geom_point() + geom_line() + 
  geom_vline(xintercept = 1e3, lty = 'dashed', col = 'blue2') + 
  annotate(geom = 'text', x = 1.4e3, y = 0.02, 
           label = 'N = 1000', angle = 90)

gexploracion14 <- exploracion2 %>% 
  ggplot(aes(x = N, y = RMSE)) +
  geom_point() + geom_line() +
  geom_vline(xintercept = 1e3, lty = 'dashed', col = 'blue2') + 
  annotate(geom = 'text', x = 1.4e3, y = 0.02, 
           label = 'N = 1000', angle = 90)


gexploracion1T <- gexploracion11 + gexploracion12 + gexploracion13 + gexploracion14 +
  plot_layout(ncol = 2) + 
  plot_annotation(
    title = 'N. de individuos para encontrar PTA',
    subtitle = 'Simulaciones Cefepime 2g q8h tinf 30min')


ggsave('002_Convergencia_NpuntosPerfil.pdf', gexploracion1T, 'pdf', 'figures', 1, 
       8, 6)

# Se requieren por lo menos 1000
#'-------------------------------------------------------------------------------
# 3. Número de puntos en la malla ------------------
#'-------------------------------------------------------------------------------

N = 1e2

MIC_vector <- c(1 * (2 ^ (seq(1, 3, 1))))

vecMalla <- c(sapply(1:2, function(x) 1:9 * 10 ^ x))
adm <- list(time = seq(0, 8 * 6, by = 8), amount = 2000, tinf = 0.5)

group <- list()

for (i in seq_along(vecMalla)) {
  
  out <- list(name = 'y_1', 
              time = seq(tail(adm$time, 4)[1], tail(adm$time, 4)[4], 
                         length.out = vecMalla[i]))
  
  g1 <- list(size = N, level = 'individual',
             treatment = adm, parameter = p, output = out)
  
  group[[i]] <- g1
}

# Simulación
ptm <- proc.time()
res <- simulx(model = mdir, group  = group)
print(proc.time() - ptm)

nDosis_DF1 <- res$y_1 %>% 
  as_tibble() %>% 
  group_by(group) %>% 
  nest() %>% 
  mutate(mat = map(data, ~conv_matriz(.x, values = y_1))) %>% 
  mutate(PTA = map(mat, ~pta_verificador(.x, MIC_vector, crit = 0.5, treshold = 1))) %>% 
  mutate(
    fTmasMIC = map(PTA, ~ magrittr::use_series(.x, 'fTmasMIC')),
    PTA      = map(PTA, ~ magrittr::use_series(.x, 'PTA'))
  )

dfExpo2 <- nDosis_DF1 %>% 
  select(group, fTmasMIC) %>% 
  unnest(fTmasMIC) %>% 
  filter(MIC == 8) %>% 
  pivot_longer(cols = matches('^ID\\.')) %>% 
  pivot_wider(names_from = MIC, values_from = value)

truePTA <- dfExpo2 %>% filter(group == 12) %>% 
  pull(`8`) %>% mean()


exploracion1 <- dfExpo2 %>%
  group_by(group) %>%
  nest() %>%
  mutate(M = map(data, ~ extraerMuestreos(.x, '8', 30, 100)),
         S = map(M, ~ summarise(.x, across(
           everything(), ~ mean(.x)
         )))) %>% 
  select(-M) %>% 
  unnest(S) %>% 
  pivot_longer(cols = matches('^rep')) %>% 
  mutate(lMalla = vecMalla[as.integer(group)])

exploracion1

gexploracion11 <- exploracion1 %>% 
  ggplot(aes(x = lMalla, y = value)) + 
  geom_point() + 
  ylab('PTA (fT>MIC; MIC = 8)')  +
  xlab('N. línea malla')  +
  geom_vline(xintercept = 30, lty = 'dashed', col = 'blue2') + 
  geom_hline(yintercept = truePTA, lty = 'dashed', col = 'black') + 
  annotate(geom = 'text', x = 40, y = 0.42, 
           label = 'malla = 30', angle = 0, hjust = 0)


gexploracion11 +
  scale_x_continuous(trans = 'log10')

gexploracion12 <- exploracion1 %>% 
  ggplot(aes(x = group, y = value)) + 
  geom_boxplot() + 
  ylab('PTA (fT>MIC; MIC = 8)') +
  xlab('N. línea malla')  +
  
  exploracion2 <- exploracion1 %>%
  mutate(e = value - truePTA) %>%
  group_by(lMalla) %>%
  summarise(MAE = mean(abs(e)),
            RMSE = sqrt(mean(e^2)))


gexploracion13 <- exploracion2 %>%
  ggplot(aes(x = lMalla, y = MAE)) +
  geom_point() + geom_line()   +
  xlab('N. línea malla')  +
  geom_vline(xintercept = 30, lty = 'dashed', col = 'blue2') + 
  annotate(geom = 'text', x = 50, y = 0.1, 
           label = 'malla = 30', angle = 90, hjust = 0)

gexploracion14 <- exploracion2 %>%
  ggplot(aes(x = lMalla, y = RMSE)) +
  geom_point() + geom_line()   +
  xlab('N. línea malla')  +
  geom_vline(xintercept = 30, lty = 'dashed', col = 'blue2') + 
  annotate(geom = 'text', x = 50, y = 0.1, 
           label = 'malla = 30', angle = 90, hjust = 0)

gexploracion14


gexploracion1T <- gexploracion11 + gexploracion12 + gexploracion13 + gexploracion14 +
  plot_layout(ncol = 2) + 
  plot_annotation(
    title = 'N. de líneas en malla para  fT>MIC',
    subtitle = 'Simulaciones Cefepime 2g q8h tinf 30min')

gexploracion1T

ggsave('003_Convergencia_NMalla.pdf', gexploracion1T, 'pdf', 'figures', 1, 
       8, 6)