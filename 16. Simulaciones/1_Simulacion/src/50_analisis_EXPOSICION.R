##------------------------------------------------------------------------------#
## Nombre del Script: análisis de exposición de simulaciones --------------------
##  
## Propósito del Script: análisis de exposición de simulaciones
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  18-06-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
require(tidyverse)


EXPO1 <- readRDS('results/L001/EXPO1.rds')
EXPO2 <- readRDS('results/L002/EXPO1.rds')
EXPO3 <- readRDS('results/L003/EXPO1.rds')
EXPO4 <- readRDS('results/L004/EXPO1.rds')
EXPO5 <- readRDS('results/L005/EXPO1.rds')
EXPO6 <- readRDS('results/L006/EXPO1.rds')

labels <- read_csv('data/20_etiquetas_SIM.csv', quote = "\"")

EXPO <- bind_rows(EXPO1, EXPO2, EXPO3,
                  EXPO4, EXPO5, EXPO6)

EXPOt <- EXPO %>% 
  add_column(L = rep(x=paste0('L00', seq(1:6)), each=4), .before = 'group') %>% 
  ungroup() %>% 
  mutate(group = as.double(group)) %>% 
  left_join(labels, by = c('L', 'group')) %>% 
  mutate(etiqueta1 = str_replace(etiqueta, '\\sINI|\\sSS', ''))

#-------------------------------------------------------------------------------#
# Análisis de AUC ---------------------------------------------------------------
#-------------------------------------------------------------------------------#

G1 <- EXPOt %>% 
  select(-data, -mat, -tmax_cmax) %>%
  unnest(expo) %>%
  ggplot(aes(x = AUC, col = etiqueta1)) +
  geom_density() + 
  theme_bw() +
  facet_grid(estado ~ dd) +
  xlab('AUC (mg*h/L)') + ylab('') +
  theme(legend.title = element_blank(), 
        panel.grid = element_line(colour = 'gray98'))


#-------------------------------------------------------------------------------#
# Análisis de Cmax --------------------------------------------------------------
#-------------------------------------------------------------------------------#
G2 <- EXPOt %>% 
  select(-data, -mat, -expo) %>%
  unnest(tmax_cmax) %>%
  ggplot(aes(x = Cmax, col = etiqueta1)) +
  geom_density() + 
  theme_bw() +
  facet_grid(estado ~ dd) +
  xlab(expression(C[max]~"(mg/L)")) + ylab('') +
  theme(legend.title = element_blank(), 
        panel.grid = element_line(colour = 'gray98'))

#-------------------------------------------------------------------------------#
# Almacenamiento gráficos -----------------------------------------------------
#-------------------------------------------------------------------------------#
ggsave('50_análisis_AUC.pdf', G1, 'pdf', 'figures', 1, 8, 6)
ggsave('50_análisis_Cmax.pdf', G2, 'pdf', 'figures', 1, 8, 6)
