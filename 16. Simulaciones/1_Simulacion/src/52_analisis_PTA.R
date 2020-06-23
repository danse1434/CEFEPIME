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
# Simulación 100%fT>MIC ---------------------------------------------------------
#-------------------------------------------------------------------------------#
# Lectura de resultados de simulación 
RESPTA_1 <- vector('list', 6L)

for (i in 1:6) {
  RESPTA_1[[i]] <- readRDS(paste0('results/L00', i, '/RES_PTA_100tmasMIC.rds'))
}

# Crear tabla unificada
RESPTA_1 <- map_dfr(.x = RESPTA_1, .f = ~.x)

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

RESPTA_1 <- RESPTA_1 %>% 
  add_column(L = rep(x=paste0('L00', seq(1:6)), each=4), .before = 'group') %>% 
  ungroup() %>% 
  mutate(group = as.double(group)) %>% 
  left_join(labels, by = c('L', 'group')) %>% 
  mutate(etiqueta1 = str_replace(etiqueta, '\\sINI|\\sSS', ''),
         etiqueta2 = glue::glue('q{frec}h, tinf {dur}h') ) 


#-------------------------------------------------------------------------------#
# >> Percentiles Indice PK-PD por MIC  ------------------------------------------
#-------------------------------------------------------------------------------#
# 1 Seleccionar columnas de interés
# 2 Adicionar columna-lista con valores de percentil
# 3 Determinar valores de percentil por MIC
#................................................................................

REPSTA_2 <- RESPTA_1 %>%
  select(L, group, fTmasMIC, dosis, frec, dur, estado, dd, etiqueta1, etiqueta2) %>%
  add_column(perc_vec = list(perc_vec)) %>% 
  mutate(percs = pmap(.l = list(data = fTmasMIC, prob_vec = perc_vec),
                      .f = perc_fT))

#-------------------------------------------------------------------------------#
# >> Gráfico de datos fT > MIC ------------------------------------------------
#-------------------------------------------------------------------------------#

G1 <- REPSTA_2 %>% 
  select(-fTmasMIC, -perc_vec) %>% 
  unnest(percs) %>% unnest(percs) %>% 
  ggplot(aes(x = MIC, y = p.5, group = etiqueta2)) +
  geom_line(aes(col = etiqueta2, lty = etiqueta2)) + 
  geom_ribbon(aes(ymin = p.1, ymax = p.9, col = etiqueta2, 
                  fill = etiqueta2), alpha = 0.2) +
  # geom_ribbon(aes(ymin = p.2, ymax = p.8,
  #                 fill = after_scale(alpha(group, 0.2)))) +
  # geom_ribbon(aes(ymin = p.3, ymax = p.7,
  #                 fill = after_scale(alpha(group, 0.2)))) +
  # geom_ribbon(aes(ymin = p.4, ymax = p.9,
  #                 fill = after_scale(alpha(group, 0.2)))) +
  theme_bw() + 
  scale_x_continuous(trans = 'log2', 
                     breaks = c(1,2,4,8,16,32,64)) + 
  coord_cartesian(xlim = c(0.25, 64)) +
  facet_grid(dd ~ estado, scales = 'free_x', 
             labeller = labeller(dd = label_both, .cols = toupper)) + 
  xlab('MIC (mg/L)') + ylab('fT > MIC') + 
  theme(legend.position = 'bottom', 
        legend.title = element_blank())  

ggsave('52_análisis_100fTmasMIC_1.pdf', G1, 'pdf', 'figures', 1, 6, 7)

#-------------------------------------------------------------------------------#
# > PTA      ------------------------------------------------------
#-------------------------------------------------------------------------------#

RESPTA_1_pta <- RESPTA_1 %>% 
  select(L, group, PTA, dosis, frec, dur, estado, dd, etiqueta1, etiqueta2) %>% 
  unnest(PTA)
#-------------------------------------------------------------------------------#
# >> Gráfico de datos PTA ----------------------------------------------------
#-------------------------------------------------------------------------------#
G2 <- RESPTA_1_pta %>% 
  ggplot(aes(x = MIC, y = PTA, col = etiqueta2)) + 
  geom_line() + 
  geom_point(data = filter(RESPTA_1_pta, MIC %in% MIC_vec)) +
  scale_x_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) + 
  theme_bw() +
  xlab('MIC (mg/L)') + ylab('PTA') + 
  geom_hline(yintercept = 0.9) +
  facet_grid(dd ~ estado, scales = 'free_x',
             labeller = labeller(dd = label_both, .cols = toupper))+ 
  theme(legend.position = 'bottom', 
        legend.title = element_blank()) 

# Almacenamiento
ggsave('53_análisis_100fTmasMIC_PTA.pdf', G2, 'pdf', 'figures', 1, 6, 7)

#-------------------------------------------------------------------------------#
# Simulación 50%fT> 4MIC ---------------------------------------------------------
#-------------------------------------------------------------------------------#
# Lectura de resultados de simulación 
RESPTA_5 <- vector('list', 6L)

for (i in 1:6) {
  RESPTA_5[[i]] <- readRDS(paste0('results/L00', i, '/RES_PTA_50tmas4MIC.rds'))
}

# Crear tabla unificada
RESPTA_5 <- map_dfr(.x = RESPTA_5, .f = ~.x)

RESPTA_5 <- RESPTA_5 %>% 
  add_column(L = rep(x=paste0('L00', seq(1:6)), each=4), .before = 'group') %>% 
  ungroup() %>% 
  mutate(group = as.double(group)) %>% 
  left_join(labels, by = c('L', 'group')) %>% 
  mutate(etiqueta1 = str_replace(etiqueta, '\\sINI|\\sSS', ''),
         etiqueta2 = glue::glue('q{frec}h, tinf {dur}h') ) 

REPSTA_5b <- RESPTA_5 %>%
  select(L, group, fTmasMIC, dosis, frec, dur, estado, dd, etiqueta1, etiqueta2) %>%
  add_column(perc_vec = list(perc_vec)) %>% 
  mutate(percs = pmap(.l = list(data = fTmasMIC, prob_vec = perc_vec),
                      .f = perc_fT))


G5b <- REPSTA_5b %>% 
  select(-fTmasMIC, -perc_vec) %>% 
  unnest(percs) %>% unnest(percs) %>% 
  ggplot(aes(x = MIC, y = p.5, group = etiqueta2)) +
  geom_line(aes(col = etiqueta2, lty = etiqueta2)) + 
  geom_ribbon(aes(ymin = p.1, ymax = p.9, col = etiqueta2, 
                  fill = etiqueta2), alpha = 0.2) +
  theme_bw() + 
  scale_x_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) + 
  coord_cartesian(xlim = c(2^-10, 64)) +
  facet_grid(dd ~ estado, scales = 'free_x', 
             labeller = labeller(dd = label_both, .cols = toupper)) + 
  xlab('MIC (mg/L)') + ylab('fT > 4 * MIC') + 
  theme(legend.position = 'bottom', 
        legend.title = element_blank())  

# Almacenamiento
ggsave('54_análisis_50fTmas4MIC_1.pdf', G5b, 'pdf', 'figures', 1, 6, 7)

#-------------------------------------------------------------------------------#
# > PTA      ------------------------------------------------------
#-------------------------------------------------------------------------------#
RESPTA_5_pta <- RESPTA_5 %>% 
  select(L, group, PTA, dosis, frec, dur, estado, dd, etiqueta1, etiqueta2) %>% 
  unnest(PTA)

#-------------------------------------------------------------------------------#
# >> Gráfico de datos PTA ----------------------------------------------------
#-------------------------------------------------------------------------------#
G5c <- RESPTA_5_pta %>% 
  ggplot(aes(x = MIC, y = PTA, col = etiqueta2)) + 
  geom_line() + 
  geom_point(data = filter(RESPTA_5_pta, MIC %in% MIC_vec)) +
  scale_x_continuous(trans = log2_trans(), 
                     breaks = trans_breaks("log2", function(x) 2^x),
                     labels = trans_format("log2", math_format(2^.x))) + 
  theme_bw() +
  xlab('MIC (mg/L)') + ylab('PTA') + 
  geom_hline(yintercept = 0.9) +
  facet_grid(dd ~ estado, scales = 'free_x',
             labeller = labeller(dd = label_both, .cols = toupper))+ 
  theme(legend.position = 'bottom', 
        legend.title = element_blank()) 

# Almacenamiento
ggsave('55_analisis_50fTmas4MIC_PTA.pdf', G5c, 'pdf', 'figures', 1, 6, 7)

#-------------------------------------------------------------------------------#
# Sección de tablas -----------------------------------------------------
#-------------------------------------------------------------------------------#

#' Función de formato condicional
#' Esta función sólo cuenta con dos reglas: (1) PTA > 0.90, (2) PTA > 0.80 y < 0.90
#'
#' @param data: datos
#' @param param: columna seleccionada para formato condicionado
#' @param crit1: criterio 1
#' @param crit2: criterio 2 
#'
#' @return
#' @export
#'
#' @examples
fun_param <- function(data, param, crit1 = 0.80, crit2 = 0.90) {
  param_quo <- rlang::ensym(param)
  
  z = data %>% 
    tab_style(style = list(cell_fill(color = 'lightgreen'),
                           cell_text(weight = 'bold')), 
              locations = cells_body(columns = vars(!!param_quo), 
                                     rows = (!!param_quo >= crit2))
    ) %>% 
    tab_style(style = list(cell_fill(color = 'thistle1'),
                           cell_text(weight = 'bold')), 
              locations = cells_body(columns = vars(!!param_quo), 
                                     rows = (!!param_quo >= crit1 & !!param_quo < crit2))
    )
  return(z)
}

#-------------------------------------------------------------------------------#
# > Simul. 100%fT>MIC-----------------------------------------------------
#-------------------------------------------------------------------------------#
FIM_RESPTA_1 <- filter(RESPTA_1_pta, 
       MIC %in% MIC_vec & MIC > 0.01) %>% 
  pivot_wider(id_cols = c('L', 'group','dd', 'dosis', 'etiqueta2', 'estado'), 
              names_from = MIC, 
              values_from = PTA) %>%
  arrange(dd, dosis, etiqueta2, estado) %>% 
  gt() %>% 
  tab_header(title = md('**Probabilidad de alcanzar el objetivo PK-PD (PTA)**'),
             subtitle = md('**Objetivo**: 100% <i>f</i>T <sub>>MIC</sub>; **Creatinina Sérica**: 0.54 mg/dL')) %>% 
  tab_spanner(label = md('**MIC (mg/L)**'), 7:23) %>%
  cols_merge(1:2) %>% 
  cols_label(
    dd        = 'DD (md)',
    dosis     = 'Dosis (mg)',
    etiqueta2 = 'Dosificación',
    estado    = 'Estado'
  ) %>% 
  fmt_number(7:23, decimals = 3) %>% 
  text_transform(locations = cells_body(vars(estado)),
                 fn = toupper) %>% 
  tab_options(
    column_labels.font.size = "smaller",
    table.font.size = "smaller",
    data_row.padding = px(3)
  ) %>% 
  fun_param(`0.015625`) %>% 
  fun_param(`0.03125`) %>% fun_param(`0.0625`) %>%
  fun_param(`0.125`) %>% fun_param(`0.25`) %>% 
  fun_param(`0.5`) %>%
  fun_param(`1`) %>% fun_param(`2`) %>% 
  fun_param(`4`) %>%
  fun_param(`8`) %>% fun_param(`16`) %>%
  cols_hide(vars('L'))

#-------------------------------------------------------------------------------#
# > Simul. 50%fT>4MIC-----------------------------------------------------
#-------------------------------------------------------------------------------#
FIM_RESPTA_2 <- filter(RESPTA_5_pta, 
       MIC %in% MIC_vec & MIC > 0.01) %>% 
  pivot_wider(id_cols = c('L', 'group','dd', 'dosis', 'etiqueta2', 'estado'), 
              names_from = MIC, 
              values_from = PTA) %>%
  arrange(dd, dosis, etiqueta2, estado) %>% 
  gt() %>% 
  tab_header(title = md('**Probabilidad de alcanzar el objetivo PK-PD (PTA)**'),
             subtitle = md('**Objetivo**: 50% <i>f</i>T <sub>>4xMIC</sub>; **Creatinina Sérica**: 0.54 mg/dL')) %>% 
  tab_spanner(label = md('**MIC (mg/L)**'), 7:23) %>%
  cols_merge(1:2) %>% 
  cols_label(
    dd        = 'DD (md)',
    dosis     = 'Dosis (mg)',
    etiqueta2 = 'Dosificación',
    estado    = 'Estado'
  ) %>% 
  fmt_number(7:23, decimals = 3) %>% 
  text_transform(locations = cells_body(vars(estado)),
                 fn = toupper) %>% 
  tab_options(
    column_labels.font.size = "smaller",
    table.font.size = "smaller",
    data_row.padding = px(3)
  ) %>% 
  fun_param(`0.015625`) %>% 
  fun_param(`0.03125`) %>% fun_param(`0.0625`) %>%
  fun_param(`0.125`) %>% fun_param(`0.25`) %>% 
  fun_param(`0.5`) %>%
  fun_param(`1`) %>% fun_param(`2`) %>% 
  fun_param(`4`) %>%
  fun_param(`8`) %>% fun_param(`16`) %>%
  cols_hide(vars('L'))


FIM_RESPTA_1 %>% 
  gtsave(filename = '56_PTA_1.html', path = file.path(getwd(), 'figures'))

FIM_RESPTA_2 %>% 
  gtsave(filename = '57_PTA_2.html', path = file.path(getwd(), 'figures'))
