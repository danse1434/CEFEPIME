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
PTA_50fTmas4MIC <- vector('list', 6L)

SCR_VEC <- c('07', '10', '13', '16', '19', '40')

for (i in 1:6) {
  PTA_100fTmasMIC[[i]] <-
    readRDS(paste0('results/L0', SCR_VEC[[i]], '/RES_PTA_100tmasMIC.rds'))
  
  PTA_50fTmas4MIC[[i]] <-
    readRDS(paste0('results/L0', SCR_VEC[[i]], '/RES_PTA.rds'))

}

# Crear tabla unificada
PTA_100fTmasMIC <- map_dfr(.x = PTA_100fTmasMIC, .f = ~.x) %>% 
  add_column(Indicador = '100%fT>MIC')

PTA_50fTmas4MIC <- map_dfr(.x = PTA_50fTmas4MIC, .f = ~.x) %>% 
  add_column(Indicador = '50%fT>4MIC')

RESPTA_1 <-
  bind_rows(PTA_100fTmasMIC, PTA_50fTmas4MIC)

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
  unnest(PTA)

#-------------------------------------------------------------------------------#
# > Gráfico -----------------------------------------------------
#-------------------------------------------------------------------------------#
G_compar_renal_PTA <- RESPTA_2_pta %>%
  ggplot(aes(MIC, PTA, col = ID)) +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept = 0.9, lty = 'dashed') +
  geom_point(data = filter(RESPTA_2_pta, MIC %in% MIC_vec)) +
  scale_x_continuous(
    trans = log2_trans(),
    breaks = trans_breaks("log2", function(x)
      2 ^ x),
    labels = trans_format("log2", math_format(2 ^ .x))
  ) +
  coord_cartesian(xlim = c(2 ^ -10, 64)) +
  facet_grid( ~ Indicador) +
  xlab('MIC (mg/L)') + ylab('PTA') +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) 

ggsave('58_comparativo_func_renal_PTA_CEP2gq8h.pdf', G_compar_renal_PTA, 'pdf', 
       'figures', 1, 7, 5)

#-------------------------------------------------------------------------------#
# > Indicador PK-PD ------------------------------------------------
#-------------------------------------------------------------------------------#
RESPTA_2_indic <- RESPTA_2 %>% 
  filter(group == 1) %>% 
  unnest(percs) %>% unnest(percs) 

G_compar_renal_Indic <- RESPTA_2_indic %>% 
  ggplot(aes(MIC, p.5, col = ID)) +
  geom_line() +
  theme_bw() +
  geom_hline(yintercept = 0.9, lty = 'dashed') +
  # geom_ribbon(aes(ymin = p.1, ymax = p.9, col = ID, fill = ID), alpha = 0.2) +
  geom_point(data = filter(RESPTA_2_indic, MIC %in% MIC_vec)) +
  scale_x_continuous(
    trans = log2_trans(),
    breaks = trans_breaks("log2", function(x)
      2 ^ x),
    labels = trans_format("log2", math_format(2 ^ .x))
  ) +
  coord_cartesian(xlim = c(2 ^ -10, 64)) +
  facet_grid( ~ Indicador) +
  xlab('MIC (mg/L)') + ylab('Indicador') +
  theme(legend.position = 'bottom',
        legend.title = element_blank()) 

ggsave('59_comparativo_func_renal_Indic_CEP2gq8h.pdf', G_compar_renal_Indic, 'pdf', 
       'figures', 1, 7, 5)

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
FIM_RESPTA_2 <- RESPTA_2_pta %>% 
  select(ID, group, MIC, PTA, Indicador) %>% 
  filter((MIC %in% MIC_vec) & (MIC > 0.01)) %>% 
  pivot_wider(id_cols = c('Indicador', 'ID','group'), 
              names_from = MIC, 
              values_from = PTA) %>% 
  gt() %>% 
  tab_options(
    column_labels.font.size = "smaller",
    table.font.size = "smaller",
    data_row.padding = px(3)
  ) %>% 
  tab_header(title = md('**Probabilidad de alcanzar el objetivo PK-PD (PTA)**'),
             subtitle = md('**Régimen**: Cefepime 2000mg q8h en bolo itermitente (tinf. 30 min)')) %>%
  tab_spanner(label = md('**MIC (mg/L)**'), 4:20) %>%
  cols_label(ID=md('S<sub>CR</sub> (mg/dL)')) %>% 
  fmt_number(4:20, decimals = 3) %>%
  fun_param(`0.015625`) %>% 
  fun_param(`0.03125`) %>% fun_param(`0.0625`) %>%
  fun_param(`0.125`) %>% fun_param(`0.25`) %>% 
  fun_param(`0.5`) %>%
  fun_param(`1`) %>% fun_param(`2`) %>% 
  fun_param(`4`) %>%
  fun_param(`8`) %>% fun_param(`16`) %>%
  cols_hide(vars('group')) %>% 
  text_transform(locations = cells_body(vars(ID)),
                 fn = function(x) str_replace(x, "SCR", "")) 

FIM_RESPTA_2 %>% 
  gtsave(filename = '60_PTA_SCR.html', path = file.path(getwd(), 'figures'))


 