##------------------------------------------------------------------------------#
## Nombre del Script:  ----------------------------------------------------------
##  
## Propósito del Script:  
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
require(tidyverse)

PRCT_1 <- readRDS('results/L001/PRCT1.rds')$y
PRCT_2 <- readRDS('results/L002/PRCT1.rds')$y
PRCT_3 <- readRDS('results/L003/PRCT1.rds')$y
PRCT_4 <- readRDS('results/L004/PRCT1.rds')$y
PRCT_5 <- readRDS('results/L005/PRCT1.rds')$y
PRCT_6 <- readRDS('results/L006/PRCT1.rds')$y

labels <- read_csv('data/20_etiquetas_SIM.csv', quote = "\"")

PRCT <- bind_rows(PRCT_1, PRCT_2, PRCT_3,
                  PRCT_4, PRCT_5, PRCT_6)


PRCT_T <- PRCT %>%
  add_column(L = rep(x=paste0('L00', seq(1:6)), each=250*4), .before = 'time') %>% 
    mutate(group = as.double(group)) %>% 
    left_join(labels, by = c('L', 'group')) %>% 
  mutate(etiqueta1 = str_replace(etiqueta, '\\sINI|\\sSS', ''),
         etiqueta2 = glue::glue('q{frec}h, tinf {dur}h')) %>% 
  as_tibble(.name_repair = 'universal') 
  

G1 <- PRCT_T %>%
  ggplot(aes(x = time, y = ..50.)) + 
  geom_line(aes(col = etiqueta2, lty = etiqueta2)) +
  geom_ribbon(aes(ymin=..10., ymax=..90., fill=etiqueta2), alpha = 0.2) +
  theme_bw() + 
  facet_grid(dd ~ estado, scales = 'free_x', 
             labeller = labeller(dd = label_both, .cols = toupper)) + 
  xlab('Tiempo (h)') + ylab('Concentración (mg/L)') + 
  theme(legend.position = 'bottom', 
        legend.title = element_blank())

G2 <- PRCT_T %>%
  ggplot(aes(x = time, y = ..50.)) + 
  geom_line(aes(col = etiqueta2, lty = etiqueta2)) +
  geom_ribbon(aes(ymin=..10., ymax=..90., fill=etiqueta2), alpha = 0.2) +
  geom_ribbon(aes(ymin=..20., ymax=..80., fill=etiqueta2), alpha = 0.3) +
  geom_ribbon(aes(ymin=..30., ymax=..70., fill=etiqueta2), alpha = 0.4) +
  geom_ribbon(aes(ymin=..40., ymax=..60., fill=etiqueta2), alpha = 0.6) +
  theme_bw() + 
  facet_grid(dd ~ estado, scales = 'free_x', 
             labeller = labeller(dd = label_both, .cols = toupper)) + 
  xlab('Tiempo (h)') + ylab('Concentración (mg/L)') + 
  theme(legend.position = 'bottom', 
        legend.title = element_blank())

if (!file.exists('./figures/51_analisis_regimenes_1.pdf')) {
  ggsave('51_analisis_regimenes_1.pdf', G1, 'pdf', 'figures', 1, 6, 8)
  ggsave('51_análisis_regimenes_2.pdf', G2, 'pdf', 'figures', 1, 6, 8)
}


# require(plotly)
