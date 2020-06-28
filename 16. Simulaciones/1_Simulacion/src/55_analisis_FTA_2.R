##------------------------------------------------------------------------------#
## Nombre del Script: Comparación de FTA terapia dirigida vs empírica --------
##  
## Propósito del Script: análisis de FTA (Alcance de Objetivo Fraccional)
## En este script se utilizaron los datos de PTA obtenidos con SCR desconocido
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 24-06-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
# Carga de paquetes
require(gt)
require(magrittr)
require(scales)
require(readxl)
require(tidyverse)

# Vector con puntos
MIC_vec = c(1 * (2 ^ (seq(-10, 10, 1))))

#-------------------------------------------------------------------------------#
# Lectura EUCAST ------------------------------------------------------
#-------------------------------------------------------------------------------#
# Lectura de archivo de EUCAST para cefepime
eucast_df <- read_excel('data/EUCAST/EUCAST_Wild_Type.xlsx', 
                        sheet = 'EUCAST_CEP', skip = 4, na = 'ND', 
                        .name_repair = "minimal")


# Cálculo de Fracción en tratamiento empírico -----
#' 1 Colapsar tabla de datos con columna de MIC, y columna de N
#' 2 Agrupar por microorganismos
#' 3 Cambiar el MIC por la secuencia de potencias de 2 correcta _MIC_
#' 4 Calcular frecuencia de aislamientos por microorganismo vs MIC _Frc_
#' 
eucast_df1 <- eucast_df %>% 
  pivot_longer(cols = matches("[0-9]"), 
               names_to = 'MIC', values_to = 'Number') %>% 
  group_by(Microorganism) %>% 
  mutate(MIC = 2^(-9:9),
         Frc = Number/sum(Number))

#-------------------------------------------------------------------------------#
# Tabla de Puntos de Corte
break_points_EUCAST <- tibble::tribble(
           ~Microorganismo, ~P_Corte_EUCAST_S, ~P_Corte_EUCAST_R,
        "Escherichia coli",                 1,                4L,
   "Klebsiella pneumoniae",                 1,                4L,
  "Pseudomonas aeruginosa",             0.001,                8L
)

# Cálculo de Fracción en tratamiento directo -----
#' 1 Agrupar por variables de interés
#' 2 Anidar
#' 3 Unir los valores de punto de corte de MIC
#' 4 Filtrar registros en que no haya _puntos de corte_ (se eliminar microorganismos 
#' no deseados)
#' 5 Filtrar tablas para eliminar MIC que sean iguales o mayores al punto de corte
#' 6 Desanidar la tabla filtrada _data1_  
#' 7 Agrupar por microorganismo
#' 8 Calcular la fracción 
eucast_df2 <- eucast_df1 %>% 
  group_by(Microorganism, ECOFF, Distributions, Observations) %>% 
  nest() %>% 
  left_join(break_points_EUCAST, by = c("Microorganism" = "Microorganismo")) %>% 
  filter(!is.na(P_Corte_EUCAST_S)) %>%
  mutate(data1 = map(data, ~filter(., MIC < P_Corte_EUCAST_R))) %>% 
  unnest(data1) %>% 
  group_by(Microorganism, .drop = TRUE) %>% 
  mutate(Frc = Number/sum(Number))
  
#-------------------------------------------------------------------------------#
# Lectura de archivos de PTA ------------------------------------------------
#-------------------------------------------------------------------------------#

RESPTA_1a <- readRDS(paste0('results/L040/RES_PTA_100tmasMIC.rds'))
RESPTA_2a <- readRDS(paste0('results/L040/RES_PTA_60tmasMIC.rds'))

RESPTA <- bind_rows(RESPTA_1a, RESPTA_2a) %>% 
  add_column(Indicador = rep(c('100%fTmasMIC', '60%fTmasMIC'), each = 4))

RESPTA_1b <- RESPTA %>% 
  ungroup() %>% 
  mutate(group = as.double(group)) %>% 
  unnest(PTA) %>% 
  filter(MIC %in% MIC_vec)

#-------------------------------------------------------------------------------#
# Manipulación de tablas -----------------------------------------------------
#-------------------------------------------------------------------------------#


eucast_ls1 <- eucast_df1 %>% 
  filter(Microorganism %in% break_points_EUCAST$Microorganismo) %>% 
  pivot_wider(id_cols = MIC, names_from = Microorganism, values_from = Frc) 

eucast_ls2 <- eucast_df2 %>% 
  pivot_wider(id_cols = MIC, names_from = Microorganism, values_from = Frc) 

eucast_ls3a <- RESPTA_1b[, c('Indicador', 'group', 'MIC', 'PTA')] %>% 
  left_join(eucast_ls1, by = 'MIC') %>% 
  mutate(across(.cols = matches("^\\w+\\s\\w+$"), .fns = ~.x * PTA)) %>% 
  group_by(Indicador, group, .drop = TRUE) %>% 
  summarise(across(matches("^\\w+\\s\\w+$"), ~sum(.x, na.rm = TRUE))) %>% 
  pivot_longer(cols = matches("^\\w+\\s\\w+$"), 
               names_to = 'Microorganism', 
               values_to = 'CFR')

eucast_ls3b <- RESPTA_1b[, c('Indicador', 'group', 'MIC', 'PTA')] %>% 
  left_join(eucast_ls2, by = 'MIC') %>% 
  mutate(across(.cols = matches("^\\w+\\s\\w+$"), .fns = ~.x * PTA)) %>% 
  group_by(Indicador, group, .drop = TRUE) %>% 
  summarise(across(matches("^\\w+\\s\\w+$"), ~sum(.x, na.rm = TRUE))) %>% 
  pivot_longer(cols = matches("^\\w+\\s\\w+$"), 
               names_to = 'Microorganism', 
               values_to = 'CFR')

df1 <- eucast_ls3a %>% 
  left_join(eucast_ls3b, by = c('Indicador', 'group', 'Microorganism')) %>% 
  mutate(etiqueta = case_when(
    group == 1 ~ '2g q8h tinf 30min',
    group == 2 ~ '2g q8h tinf 2h',
    group == 3 ~ '2g q8h tinf 4h',
    group == 4 ~ '6g q24h tinf 24h',
    TRUE ~ NA_character_
  )) %>% 
  rename('Empírico' = CFR.x,
         'Objetivo' = CFR.y) %>% 
  pivot_longer(cols      = c(Empírico, Objetivo), 
               names_to  = 'Terapia', 
               values_to = 'CFR')


#-------------------------------------------------------------------------------#
# Sección de tablas -----------------------------------------------------
#-------------------------------------------------------------------------------#
#' Función de formato condicional
#' Esta función sólo cuenta con dos reglas: (1) CFR > 0.85, (2) CFR > 0.70 y CFR < 0.90
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
fun_param <- function(data, param, crit1 = 0.70, crit2 = 0.85) {
  param_quo <- rlang::ensym(param)
  
  z = data %>% 
    tab_style(style = list(cell_fill(color = alpha('#00ff00',0.5))), 
              locations = cells_body(columns = vars(!!param_quo), 
                                     rows = (!!param_quo >= crit2))
    ) %>% 
    tab_style(style = list(cell_fill(color = alpha('#ffff66',0.5))), 
              locations = cells_body(columns = vars(!!param_quo), 
                                     rows = (!!param_quo >= crit1 & !!param_quo < crit2))
    ) %>% 
    tab_style(style = list(cell_fill(color = alpha('#e6b8b7',0.5))), 
              locations = cells_body(columns = vars(!!param_quo), 
                                     rows = (!!param_quo < crit1))
    )
  return(z)
}



FIM_RESPTA_1 <- df1 %>% 
  pivot_wider(names_from  = c('Microorganism', 'Terapia'), 
              values_from = 'CFR') %>% 
  select(-group) %>% 
  gt() %>% 
  fmt_number(columns = 3:8, decimals = 3) %>% 
  cols_label(
    "etiqueta"                        = 'Régimen',
    "Escherichia coli_Empírico"       = 'Empírico',
    "Escherichia coli_Objetivo"       = 'Objetivo',
    "Klebsiella pneumoniae_Empírico"  = 'Empírico',
    "Klebsiella pneumoniae_Objetivo"  = 'Objetivo',
    "Pseudomonas aeruginosa_Empírico" = 'Empírico',
    "Pseudomonas aeruginosa_Objetivo" = 'Objetivo'
  ) %>% 
  tab_spanner( label = md('*E. coli*'), 
               columns = starts_with('Escherichia') ) %>% 
  tab_spanner( label = md('*K. pneumoniae*'), 
               columns = starts_with('Klebsiella') ) %>% 
  tab_spanner( label = md('*P. aeruginosa*'), 
               columns = starts_with('Pseudomonas') ) %>% 
  tab_header(title = md('**Fracción Acumulada de Respuesta (CFR)**'),
             subtitle = md('**Creatinina Sérica**: Desconocida')) %>% 
  tab_options(
    column_labels.font.size = "smaller",
    table.font.size = "smaller",
    data_row.padding = px(3)
  ) %>% 
  fun_param("Escherichia coli_Empírico") %>% 
  fun_param("Escherichia coli_Objetivo") %>% 
  fun_param("Klebsiella pneumoniae_Empírico") %>% 
  fun_param("Klebsiella pneumoniae_Objetivo") %>% 
  fun_param("Pseudomonas aeruginosa_Empírico") %>% 
  fun_param("Pseudomonas aeruginosa_Objetivo")

# Almacenamiento de tablas
FIM_RESPTA_1 %>% 
  gtsave(filename = '69_CFR_Objetivo_Empírico.html', path = file.path(getwd(), 'figures'))


