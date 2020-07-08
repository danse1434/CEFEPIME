##------------------------------------------------------------------------------#
## Nombre del Script: Análisis de PTA - regímenes de dosificación ---------------
##  
## Propósito del Script: analizar diferentes cumplimientos de objetivos farmaco-
## terapéuticos para diferentes regímenes de dosificación de cefepime en las 
## primeras 24h y en estado estacionario.
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  18-06-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
# Carga de paquetes
require(scales)
require(gt)
require(rlang)
require(Rcpp)
require(RcppArmadillo)
require(tidyverse)

# Lectura de archivo CPP
sourceCpp('src/91_percentiles_PTA.cpp')

#-------------------------------------------------------------------------------#
# PTA para 2000mg q8h tinf. 30min ----------------------------------------------
#-------------------------------------------------------------------------------#
# Lectura de resultados de simulación 
PTA_100fTmasMIC <- vector('list', 6L)
PTA_60fTmasMIC <- vector('list', 6L)

SCR_VEC <- c('07', '10', '13', '16', '19', '40')

for (i in 1:6) {
  PTA_100fTmasMIC[[i]] <-
    readRDS(paste0('results/L0', SCR_VEC[[i]], '/RES_PTA_100tmasMIC.rds'))
  
  PTA_60fTmasMIC[[i]] <-
    readRDS(paste0('results/L0', SCR_VEC[[i]], '/RES_PTA_60tmasMIC.rds'))
  
}

# Crear tabla unificada
PTA_100fTmasMIC <- map_dfr(.x = PTA_100fTmasMIC, .f = ~.x) %>% 
  add_column(Indicador = '100%fT>MIC')

PTA_60fTmasMIC <- map_dfr(.x = PTA_60fTmasMIC, .f = ~.x) %>% 
  add_column(Indicador = '60%fT>MIC')

RESPTA_1 <-
  bind_rows(PTA_100fTmasMIC, PTA_60fTmasMIC)

# Etiquetas para los grupos de simulación
labels <- read_csv('data/20_etiquetas_SIM.csv')

# Calcular valores de percentil
perc_vec = seq(0.1, 0.9, 0.1)

# Vector con puntos
MIC_vec = c(1 * (2 ^ (seq(-10, 10, 1))))

#-------------------------------------------------------------------------------#
# > fT > MIC -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Adición de etiquetas de identificación de grupos de simulación
#................................................................................
# 1 Adicionar columna para la identificación de lotes 'L'
# 2 Desagrupar data.frame
# 3 Cambiar la columna 'group' a número
# 4 Unir la tabla de etiquetas de identificación "labels"
# 5 Eliminar el indicativo de INI o SS en etiqueta
#................................................................................

vector_ID <-
  c(paste0('SCR', format(seq(0.2, 1, 0.2), digits = 1)), 'Desconocido')

RESPTA_1 <- RESPTA_1 %>% 
  add_column(ID = rep(rep(vector_ID, each=4),2), .before = 'group') %>% 
  ungroup()

#-------------------------------------------------------------------------------#
# > Percentiles Indice PK-PD por MIC  ------------------------------------------
#-------------------------------------------------------------------------------#
# 1 Seleccionar columnas de interés
# 2 Adicionar columna-lista con valores de percentil
# 3 Determinar valores de percentil por MIC
#................................................................................

RESPTA_2 <- RESPTA_1 %>%
  add_column(perc_vec = list(perc_vec)) %>% 
  mutate(percs = pmap(.l = list(data = fTmasMIC, prob_vec = perc_vec),
                      .f = perc_fT))

RESPTA_2_pta <- RESPTA_2 %>% 
  filter(group == 1) %>% 
  unnest(PTA) %>% 
  select(-data, -mat, -fTmasMIC, -perc_vec, -percs) %>% 
  mutate(ID = case_when(
    ID == 'Desconocido' ~ 'Unknown', 
    TRUE ~ ID
  )) %>% 
  mutate(ID = str_replace(ID, '^SCR', ''),
         ID = ifelse(ID != 'Unknown', paste0(ID, ' mg/dL'), ID))


#-------------------------------------------------------------------------------#
# > Gráfico -----------------------------------------------------
#-------------------------------------------------------------------------------#
G_compar_renal_PTA <- RESPTA_2_pta  %>% 
  ggplot(aes(MIC, PTA, group = ID)) +
  geom_line(col='gray60') +
  theme_bw() +
  geom_hline(yintercept = 0.9, lty = 'dashed') +
  geom_point(data = filter(RESPTA_2_pta, MIC %in% MIC_vec), aes(shape = ID)) +
  scale_x_continuous(
    trans = log2_trans(),
    breaks = MIC_vec,
    labels = function(x) 
    {format(x, drop0trailing = T, digits = 4, nsmall = 0, trim = T, 
            scientific = F)}) +
  coord_cartesian(xlim = c(2^-3, 2^6)) +
  facet_grid( ~ Indicador) +
  scale_shape_discrete(name = 'Serum Creatinine') +
  xlab('MIC (mg/L)') + ylab('PTA') +
  guides(shape = guide_legend(nrow=1)) +
  theme(legend.position = 'bottom', 
        legend.key = element_rect(fill = alpha('gray50', 0.2), colour = 'black')) 

G_compar_renal_PTA

ggsave('PTA_plot_2_article.pdf', G_compar_renal_PTA, 'pdf',
       '../2_Articulo_SIM', 1, 8, 4, 'in')
