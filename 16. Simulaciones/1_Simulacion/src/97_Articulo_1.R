##------------------------------------------------------------------------------#
## Nombre del Script: Gráfico para artículo de AAC ---------------
##  
## Propósito del Script: gráfico de artículo para AAC
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  04-07-2020
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
require(patchwork)

# Lectura de archivo CPP
sourceCpp('src/91_percentiles_PTA.cpp')

#-------------------------------------------------------------------------------#
# Preparación gráfico ---------------------------------------------------------
#-------------------------------------------------------------------------------#
# Lectura de resultados de simulación 
RESPTA_1 <- vector('list', 6L)

for (i in c(4:6)) {
  RESPTA_1[[i-3]] <- readRDS(paste0('results/L00', i, '/RES_PTA_60tmasMIC.rds'))
  RESPTA_1[[i]]   <- readRDS(paste0('results/L00', i, '/RES_PTA_100tmasMIC.rds'))
}

# Crear tabla unificada
RESPTA_1 <- map_dfr(.x = RESPTA_1, .f = ~.x, 
                    .id = 'Lote')

{labels <- tibble::tribble(
  ~Lote, ~group,               ~Dosis,   ~Indicador,
     1L,     1L,  "2g q8h inf. 30min",  "60%fT>MIC",
     1L,     2L,     "2g q8h inf. 2h",  "60%fT>MIC",
     1L,     3L,     "2g q8h inf. 4h",  "60%fT>MIC",
     1L,     4L,   "6g q24h inf. 24h",  "60%fT>MIC",
     2L,     1L, "2g q12h inf. 30min",  "60%fT>MIC",
     2L,     2L,  "1g q6h inf. 30min",  "60%fT>MIC",
     2L,     3L,   "4g q24h inf. 24h",  "60%fT>MIC",
     2L,     4L,    "2g q12h inf. 4h",  "60%fT>MIC",
     3L,     1L, "4g q12h inf. 30min",  "60%fT>MIC",
     3L,     2L,  "2g q6h inf. 30min",  "60%fT>MIC",
     3L,     3L,   "8g q24h inf. 24h",  "60%fT>MIC",
     3L,     4L,    "4g q12h inf. 4h",  "60%fT>MIC",
     4L,     1L,  "2g q8h inf. 30min", "100%fT>MIC",
     4L,     2L,     "2g q8h inf. 2h", "100%fT>MIC",
     4L,     3L,     "2g q8h inf. 4h", "100%fT>MIC",
     4L,     4L,   "6g q24h inf. 24h", "100%fT>MIC",
     5L,     1L, "2g q12h inf. 30min", "100%fT>MIC",
     5L,     2L,  "1g q6h inf. 30min", "100%fT>MIC",
     5L,     3L,   "4g q24h inf. 24h", "100%fT>MIC",
     5L,     4L,    "2g q12h inf. 4h", "100%fT>MIC",
     6L,     1L, "4g q12h inf. 30min", "100%fT>MIC",
     6L,     2L,  "2g q6h inf. 30min", "100%fT>MIC",
     6L,     3L,   "8g q24h inf. 24h", "100%fT>MIC",
     6L,     4L,    "4g q12h inf. 4h", "100%fT>MIC"
  )}

labels <- labels %>%
  mutate(Lote  = as.character(Lote),
         group = factor(group))



# Calcular valores de percentil
perc_vec = seq(0.1, 0.9, 0.1)

# Vector con puntos
MIC_vec = c(1 * (2 ^ (seq(-10, 10, 1))))

# > fT > MIC ----
# Adición de etiquetas de identificación de grupos de simulación
#................................................................................
# 1 Adicionar columna para la identificación de lotes 'L'
# 2 Desagrupar data.frame
# 3 Cambiar la columna 'group' a número
# 4 Unir la tabla de etiquetas de identificación "labels"
# 5 Eliminar el indicativo de INI o SS en etiqueta
#................................................................................

RESPTA_1 <- RESPTA_1 %>% 
  ungroup() %>% 
  left_join(labels, by = c('Lote', 'group')) %>% 
  select(-group, -data, -mat, -fTmasMIC)



RESPTA_2 <- RESPTA_1 %>% 
  unnest(PTA)


theme_set(theme_bw() +
            theme(panel.grid.minor = element_blank(), 
                  legend.position = 'left', 
                  legend.key = element_rect(colour = 'black', fill = 'white')))

RESPTA_2a <- RESPTA_2 %>%
  filter(Lote %in% c(1, 4))
RESPTA_2b <- RESPTA_2 %>%
  filter(Lote %in% c(2, 5))
RESPTA_2c <- RESPTA_2 %>%
  filter(Lote %in% c(3, 6))


PTA_P2a <- RESPTA_2a %>% 
  ggplot(aes(x = MIC, y = PTA, group = Dosis, shape = Dosis)) +
  geom_line(col = 'gray50') +
  geom_point(data = filter(RESPTA_2a, MIC %in% MIC_vec)) + 
  scale_shape_discrete(name = 'Dose')

PTA_P2b <- RESPTA_2b %>% 
  ggplot(aes(x = MIC, y = PTA, group = Dosis, shape = Dosis)) +
  geom_line(col = 'gray50') +
  geom_point(data = filter(RESPTA_2b, MIC %in% MIC_vec)) +
  scale_shape_discrete(name = 'Dose')

PTA_P2c <- RESPTA_2c %>% 
  ggplot(aes(x = MIC, y = PTA, group = Dosis, shape = Dosis)) +
  geom_line(col = 'gray50') +
  geom_point(data = filter(RESPTA_2c, MIC %in% MIC_vec)) +
  scale_shape_discrete(name = 'Dose')

PTA_P2_comp <- 
(PTA_P2b / PTA_P2a / PTA_P2c)  +
  plot_annotation(tag_levels = 'A') & 
  geom_hline(yintercept = 0.9, lty = 'dashed') & 
  scale_x_continuous(trans = log2_trans(), 
                     breaks = MIC_vec,
                     labels = function(x) format(x, drop0trailing = TRUE, digits = 4, 
                                                 nsmall = 0, trim = TRUE, 
                                                 scientific = FALSE),
                     guide = guide_axis(n.dodge = 2)) &
  coord_cartesian(xlim = c(2^-3, 2^6)) &
  facet_grid(~ Indicador) &
  xlab('MIC (mg/L)') 


ggsave('PTA_plot_1_article.pdf', PTA_P2_comp, 'pdf', '../2_Articulo_SIM', 
       1, 10*700/700, 10*600/700, 'in')
# 537