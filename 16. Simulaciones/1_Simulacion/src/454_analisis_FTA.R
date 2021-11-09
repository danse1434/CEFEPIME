##------------------------------------------------------------------------------#
## Nombre del Script: Análisis de FTA -----------------------------------------
##  
## Propósito del Script: análisis de FTA (Alcance de Objetivo Fraccional)
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


# > Gráficos distribuciones -----------------------------------------------

G1 <- eucast_df %>% 
  ggplot(aes(x = Observations, y = fct_reorder(Microorganism, Observations))) +
  geom_bar(stat = 'identity') +
  theme_bw() + 
  xlab('Observaciones') + ylab('Microrganismo') +
  labs(title = 'Número de observaciones EUCAST', 
       subtitle = 'Fármaco: Cefepime')

G2 <- eucast_df1 %>% 
  filter(Observations > 1000) %>% 
  filter(MIC > 0.001 & MIC < 513) %>% 
  mutate(MIC = round(MIC, digits = 3),
         ECOFF = round(ECOFF, digits = 3)) %>% 
  ggplot(aes(x = factor(MIC), y = Frc)) +
  geom_bar(stat = 'identity', fill = 'gray', col = 'black') +
  theme_bw() +
  geom_vline(aes(xintercept = factor(ECOFF)), lty= 'dashed', col='blue4') + 
  facet_wrap(. ~ Microorganism) +
  # scale_x_discrete(guide = guide_axis(n.dodge = 2)) +
  # scale_x_continuous(trans = log2_trans(), 
  #                    breaks = trans_breaks("log2", function(x) 2^x),
  #                    labels = trans_format("log2", math_format(2^.x))) +
  ylab('Proporción de aislamientos') +
  labs(title = 'Distribuciones EUCAST frente a cefepime',
       subtitle = '') + 
  theme(axis.text.x = element_text(angle = 90))

# # Almacenamiento de gráficos
# ggsave('61_numero_aislamientos_EUCAST_CEP.pdf', G1, 'pdf', 'figures', 1, 
#        10, 9)
# ggsave('62_distribuciones_EUCAST_CEP.pdf', G2, 'pdf', 'figures', 1, 6*2, 5*2)

#-------------------------------------------------------------------------------#
# Indicador 1: 100%fT>MIC ------------------------------
#-------------------------------------------------------------------------------#
# Lectura de resultados de simulación 
RESPTA_1 <- vector('list', 6L)

id_SRC <- c('L001', 'L002', 'L401', 'L004', 'L005', 'L402')

for (i in seq_along(id_SRC)) {
  RESPTA_1[[i]] <- readRDS(paste0('results/', id_SRC[i], '/RES_PTA_100tmasMIC.rds'))
}

# Crear tabla unificada
RESPTA_1 <- map_dfr(.x = RESPTA_1, .f = ~.x, .id = 'ID')

# Etiquetas para los grupos de simulación
labels <- read_csv('data/20_etiquetas_SIM.csv')

# Calcular valores de percentil
perc_vec = seq(0.1, 0.9, 0.1)

#-------------------------------------------------------------------------------#
#' Crear a la tabla **RESPTA_1b** con identificadores
#' 1 Adicionar columna de identificación de lotes
#' 2 Desagrupar data.frame
#' 3 Convertir a la columna _group_ en número
#' 4 Unir a los identificadores de régimen de dosificación en *labels*
#' 5 Crear la _etiqueta1_ sin identificación de estado
#' 6 Crear la _etiqueta2_ con sólo el fraccionamiento de la dosis diaria
#' 7 Desanidar a _PTA_
#' 8 Filtrar a MIC en valores múltiplos de 2^
#' 
RESPTA_1b <- RESPTA_1 %>%
  mutate(L = id_SRC[as.integer(ID)]) %>% 
  ungroup() %>% 
  mutate(group = as.double(group)) %>% 
  left_join(labels, by = c('L', 'group')) %>% 
  mutate(etiqueta1 = str_replace(etiqueta, '\\sINI|\\sSS', ''),
         etiqueta2 = glue::glue('q{frec}h, tinf {dur}h') ) %>% 
  unnest(PTA) %>% 
  filter(MIC %in% MIC_vec)

#-------------------------------------------------------------------------------#
#' Crear tabla **CFR_MO** con filas con régimen de dosificación, columnas por 
#' microorganismo, y valores de CFR en las celdas.
#' 
#' 1 Expandir a la tabla *eucast_df1* con los _microorganismos_ en columnas, 
#' _MIC_ en las filas, y _Frc_ en las celdas (fracción poblacional de MO)
#' 2 Unir por derecha a la tabla *RESPTA_1b* por _MIC_
#' 3 Eliminar a columnas-listas _-data_, _-mat_, y _-fTmasMIC_
#' 4 Calcular _CFRi_, al multiplicar frecuencias (en la columna de cada MO) 
#' por _PTA_
#' 5 Agrupar por identificadores de régimen de dosificación
#' 6 Sumar todos _CFRi_ para obtener _CFR_
#' 
CFR_MO <- eucast_df1 %>% 
  pivot_wider(id_cols = c(MIC),
              names_from = Microorganism, values_from = Frc) %>% 
  right_join(RESPTA_1b, by = 'MIC') %>% 
  select(-data, -mat, -fTmasMIC) %>% 
  # Cálculo de _CFR_
  mutate(across(.cols = matches("^\\w+\\s\\w+$"), .fns = ~.x * PTA)) %>% 
  group_by(etiqueta, etiqueta1, etiqueta2, dd, dosis, frec, dur, estado) %>% 
  summarise(across(matches("^\\w+\\s\\w+$"), ~sum(.x, na.rm = TRUE)))

# Vector con órden de microorganismo en la tabla
orden_MO <- c("Escherichia coli", "Klebsiella pneumoniae", "Klebsiella aerogenes", 
              "Enterobacter cloacae", "Moraxella catarrhalis", 
              "Pseudomonas aeruginosa", "Acinetobacter spp", 
              "Acinetobacter baumannii", "Stenotrophomonas maltophilia", 
              "Staphylococcus aureus", "Staphylococcus epidermidis", 
              "Streptococcus agalactiae", "Enterococcus faecalis", 
              "Enterococcus faecium")

#-------------------------------------------------------------------------------#
#' Crear tabla **RESPTA_1c** con filas con microorganismos, columnas por 
#' régimen de dosificación, y valores de CFR en las celdas.
#' 
#' 1 Desagrupar a *CFR_MO*
#' 2 Filtrar sólo regímenes de dosificación con estado estacionario "ss"
#' 3 Eliminar variables redundantes
#' 4 Colapsar la tabla con una columna para los microorganismos, y otras para 
#' los CFR
#' 5 Expandir la tabla con los regímenes como columnas y valores de celda igual 
#' a CFR
#' 6 Seleccionar a los microorganismos en *orden_MO* como valores de _MO_
#' 7 Convertir a _MO_ en un factor ordenado
#' 8 Ordenar al data.frame por _MO_
#'
RESPTA_1c <- ungroup(CFR_MO) %>% 
  filter(estado == 'ss') %>% 
  select(-etiqueta, -etiqueta2, -estado, -dd) %>% 
  pivot_longer(cols = matches("^\\w+\\s\\w+$"), 
               names_to = 'MO', values_to = 'CFR') %>% 
  pivot_wider(id_cols = 'MO',
              names_from = etiqueta1, values_from = CFR) %>% 
  filter(MO %in% orden_MO) %>% 
  mutate(MO = factor(MO, levels = orden_MO)) %>% 
  arrange(MO) 

#-------------------------------------------------------------------------------#
# Indicador 2: 60%fT>MIC  ------
#-------------------------------------------------------------------------------#
# Lectura de resultados de simulación 
RESPTA_60a <- vector('list', 6L)

for (i in seq_along(id_SRC)) {
  RESPTA_60a[[i]] <- readRDS(paste0('results/', id_SRC[i], '/RES_PTA_60tmasMIC.rds'))
}

# Crear tabla unificada
RESPTA_60a <- map_dfr(.x = RESPTA_60a, .f = ~.x, .id = 'ID')

#-------------------------------------------------------------------------------#
#' Crear a la tabla **RESPTA_1b** con identificadores
#' 1 Adicionar columna de identificación de lotes
#' 2 Desagrupar data.frame
#' 3 Convertir a la columna _group_ en número
#' 4 Unir a los identificadores de régimen de dosificación en *labels*
#' 5 Crear la _etiqueta1_ sin identificación de estado
#' 6 Crear la _etiqueta2_ con sólo el fraccionamiento de la dosis diaria
#' 7 Desanidar a _PTA_
#' 8 Filtrar a MIC en valores múltiplos de 2^
#' 
RESPTA_60b <- RESPTA_60a %>% 
  mutate(L = id_SRC[as.integer(ID)]) %>% 
  ungroup() %>% 
  mutate(group = as.double(group)) %>% 
  left_join(labels, by = c('L', 'group')) %>% 
  mutate(etiqueta1 = str_replace(etiqueta, '\\sINI|\\sSS', ''),
         etiqueta2 = glue::glue('q{frec}h, tinf {dur}h') ) %>% 
  unnest(PTA) %>% 
  filter(MIC %in% MIC_vec)

#-------------------------------------------------------------------------------#
#' Crear tabla **CFR_60_MO** con filas con régimen de dosificación, columnas por 
#' microorganismo, y valores de CFR en las celdas.
#' 
#' 1 Expandir a la tabla *eucast_df1* con los _microorganismos_ en columnas, 
#' _MIC_ en las filas, y _Frc_ en las celdas (fracción poblacional de MO)
#' 2 Unir por derecha a la tabla *RESPTA_1b* por _MIC_
#' 3 Eliminar a columnas-listas _-data_, _-mat_, y _-fTmasMIC_
#' 4 Calcular _CFRi_, al multiplicar frecuencias (en la columna de cada MO) 
#' por _PTA_
#' 5 Agrupar por identificadores de régimen de dosificación
#' 6 Sumar todos _CFRi_ para obtener _CFR_
#' 
CFR_60_MO <- eucast_df1 %>% 
  pivot_wider(id_cols = c(MIC),
              names_from = Microorganism, values_from = Frc) %>% 
  right_join(RESPTA_60b, by = 'MIC') %>% 
  select(-data, -mat, -fTmasMIC) %>% 
  # Cálculo de _CFR_
  mutate(across(.cols = matches("^\\w+\\s\\w+$"), .fns = ~.x * PTA)) %>% 
  group_by(etiqueta, etiqueta1, etiqueta2, dd, dosis, frec, dur, estado) %>% 
  summarise(across(matches("^\\w+\\s\\w+$"), ~sum(.x, na.rm = TRUE)))

#-------------------------------------------------------------------------------#
#' Crear tabla **RESPTA_60c** con filas con microorganismos, columnas por 
#' régimen de dosificación, y valores de CFR en las celdas.
#' 
#' 1 Desagrupar a *CFR_MO*
#' 2 Filtrar sólo regímenes de dosificación con estado estacionario "ss"
#' 3 Eliminar variables redundantes
#' 4 Colapsar la tabla con una columna para los microorganismos, y otras para 
#' los CFR
#' 5 Expandir la tabla con los regímenes como columnas y valores de celda igual 
#' a CFR
#' 6 Seleccionar a los microorganismos en *orden_MO* como valores de _MO_
#' 7 Convertir a _MO_ en un factor ordenado
#' 8 Ordenar al data.frame por _MO_
#'
RESPTA_60c <- ungroup(CFR_60_MO) %>% 
  filter(estado == 'ss') %>% 
  select(-etiqueta, -etiqueta2, -estado, -dd) %>% 
  pivot_longer(cols = matches("^\\w+\\s\\w+$"), 
               names_to = 'MO', values_to = 'CFR') %>% 
  pivot_wider(id_cols = 'MO',
              names_from = etiqueta1, values_from = CFR) %>% 
  filter(MO %in% orden_MO) %>% 
  mutate(MO = factor(MO, levels = orden_MO)) %>% 
  arrange(MO) 

#-------------------------------------------------------------------------------#
# Indicador 3: 50%fT>4MIC ------------------------------
#-------------------------------------------------------------------------------#
# Lectura de resultados de simulación 
RESPTA_2 <- vector('list', 6L)

for (i in seq_along(id_SRC)) {
  RESPTA_2[[i]] <- readRDS(paste0('results/', id_SRC[i], '/RES_PTA_50tmas4MIC.rds'))
}

# Crear tabla unificada
RESPTA_2 <- map_dfr(.x = RESPTA_2, .f = ~.x, .id = 'ID')

RESPTA_2b <- RESPTA_2 %>% 
  mutate(L = id_SRC[as.integer(ID)]) %>% 
  ungroup() %>% 
  mutate(group = as.double(group)) %>% 
  left_join(labels, by = c('L', 'group')) %>% 
  mutate(etiqueta1 = str_replace(etiqueta, '\\sINI|\\sSS', ''),
         etiqueta2 = glue::glue('q{frec}h, tinf {dur}h') ) %>% 
  unnest(PTA) %>% 
  filter(MIC %in% MIC_vec)

CFR_MO2 <- eucast_df1 %>% 
  pivot_wider(id_cols = c(MIC),
              names_from = Microorganism, values_from = Frc) %>% 
  right_join(RESPTA_2b, by = 'MIC') %>% 
  select(-data, -mat, -fTmasMIC) %>% 
  mutate(across(.cols = matches("^\\w+\\s\\w+$"), ~.x * PTA)) %>% 
  group_by(etiqueta, etiqueta1, etiqueta2, dd, dosis, frec, dur, estado) %>% 
  summarise(across(matches("^\\w+\\s\\w+$"), ~sum(.x, na.rm = TRUE)))


RESPTA_2c <- ungroup(CFR_MO2) %>% 
  filter(estado == 'ss') %>% 
  select(-etiqueta, -etiqueta2, -estado, -dd) %>% 
  pivot_longer(cols = matches("^\\w+\\s\\w+$"), 
               names_to = 'MO', values_to = 'CFR') %>% 
  pivot_wider(id_cols = 'MO',
              names_from = etiqueta1, values_from = CFR) %>% 
  filter(MO %in% orden_MO) %>% 
  mutate(MO = factor(MO, levels = orden_MO)) %>% 
  arrange(MO) 

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

# Tabla Ind.1 ------

dosisDefinidas <- colnames(RESPTA_1c)
dosisDefinidas <- dosisDefinidas[!str_detect(dosisDefinidas, 'MO')]


FIM_RESPTA_1 <- RESPTA_1c %>% 
  select(c("MO", dosisDefinidas)) %>% 
  gt::gt() %>% 
  tab_header(title = md('**Fracción Acumulada de Respuesta (CFR)**'),
             subtitle = md('**Objetivo**: 100% <i>f</i>T <sub>>MIC</sub>; **Creatinina Sérica**: 0.54 mg/dL')) %>% 
  tab_spanner(label = md('**Regímenes (mg/L)**'), 2:13) %>% 
  tab_spanner(label = md('**DD = 2g**'), 2:(2+3)) %>% 
  tab_spanner(label = md('**DD = 4g**'), 6:(6+3)) %>% 
  tab_spanner(label = md('**DD = 6g**'), 10:(10+3)) %>% 
  fmt_number(2:13, decimals = 3) %>% 
  cols_label(MO = md('**Microorganismo**')) %>% 
  tab_options(
    column_labels.font.size = "smaller",
    table.font.size = "smaller",
    data_row.padding = px(3)
  ) %>% 
  text_transform(locations = cells_body(columns = vars(MO)), 
                 fn = function(x) {gt::md(paste0("<em>",x,"</em>"))} ) %>% 
  cols_align(columns = vars(MO), align = "left") %>%
  fun_param("500mg q6h inf. 30min") %>% 
  fun_param("1g q12h inf. 30min") %>% 
  fun_param("2g q24h inf. 24h") %>% 
  fun_param("1g q12h inf. 4h") %>%
  fun_param("2g q12h inf. 30min") %>% 
  fun_param("1g q6h inf. 30min") %>% 
  fun_param("4g q24h inf. 24h") %>% 
  fun_param("2g q12h inf. 4h") %>% 
  fun_param("2g q8h inf. 30min") %>% 
  fun_param("2g q8h inf. 2h") %>% 
  fun_param("2g q8h inf. 4h") %>% 
  fun_param("6g q24h inf. 24h")
  
FIM_RESPTA_1

# Tabla Ind.3 ------

FIM_RESPTA_2 <- RESPTA_2c %>% 
  select(c("MO", dosisDefinidas)) %>%  
  gt::gt() %>% 
  tab_header(title = md('**Fracción Acumulada de Respuesta (CFR)**'),
             subtitle = md('**Objetivo**: 50% <i>f</i>T <sub>>4 x MIC</sub>; **Creatinina Sérica**: 0.54 mg/dL')) %>% 
  tab_spanner(label = md('**Regímenes (mg/L)**'), 2:13) %>% 
  tab_spanner(label = md('**DD = 4g**'), 2:(2+3)) %>% 
  tab_spanner(label = md('**DD = 6g**'), 6:(6+3)) %>% 
  tab_spanner(label = md('**DD = 8g**'), 10:(10+3)) %>% 
  fmt_number(2:13, decimals = 3) %>% 
  cols_label(MO = md('**Microorganismo**')) %>% 
  tab_options(
    column_labels.font.size = "smaller",
    table.font.size = "smaller",
    data_row.padding = px(3)
  ) %>%
  fun_param("500mg q6h inf. 30min") %>% 
  fun_param("1g q12h inf. 30min") %>% 
  fun_param("2g q24h inf. 24h") %>% 
  fun_param("1g q12h inf. 4h") %>%
  fun_param("2g q12h inf. 30min") %>% 
  fun_param("1g q6h inf. 30min") %>% 
  fun_param("4g q24h inf. 24h") %>% 
  fun_param("2g q12h inf. 4h") %>% 
  fun_param("2g q8h inf. 30min") %>% 
  fun_param("2g q8h inf. 2h") %>% 
  fun_param("2g q8h inf. 4h") %>% 
  fun_param("6g q24h inf. 24h") %>% 
  text_transform(locations = cells_body(columns = vars(MO)), 
                 fn = function(x) {gt::md(paste0("<em>",x,"</em>"))} ) %>% 
  cols_align(columns = vars(MO), align = "left")


# Tabla Ind.2 ---------------------------------------------------------------------------------

FIM_RESPTA_60 <- RESPTA_60c %>% 
  select(c("MO", dosisDefinidas)) %>%   
  gt::gt() %>% 
  tab_header(title = md('**Fracción Acumulada de Respuesta (CFR)**'),
             subtitle = md('**Objetivo**: 60% <i>f</i>T <sub>>MIC</sub>; **Creatinina Sérica**: 0.54 mg/dL')) %>% 
  tab_spanner(label = md('**Regímenes (mg/L)**'), 2:13) %>% 
  tab_spanner(label = md('**DD = 4g**'), 2:(2+3)) %>% 
  tab_spanner(label = md('**DD = 6g**'), 6:(6+3)) %>% 
  tab_spanner(label = md('**DD = 8g**'), 10:(10+3)) %>% 
  fmt_number(2:13, decimals = 3) %>% 
  cols_label(MO = md('**Microorganismo**')) %>% 
  tab_options(
    column_labels.font.size = "smaller",
    table.font.size = "smaller",
    data_row.padding = px(3)
  ) %>% 
  fun_param("500mg q6h inf. 30min") %>% 
  fun_param("1g q12h inf. 30min") %>% 
  fun_param("2g q24h inf. 24h") %>% 
  fun_param("1g q12h inf. 4h") %>%
  fun_param("2g q12h inf. 30min") %>% 
  fun_param("1g q6h inf. 30min") %>% 
  fun_param("4g q24h inf. 24h") %>% 
  fun_param("2g q12h inf. 4h") %>% 
  fun_param("2g q8h inf. 30min") %>% 
  fun_param("2g q8h inf. 2h") %>% 
  fun_param("2g q8h inf. 4h") %>% 
  fun_param("6g q24h inf. 24h") %>% 
  text_transform(locations = cells_body(columns = vars(MO)), 
                 fn = function(x) {gt::md(paste0("<em>",x,"</em>"))} ) %>% 
  cols_align(columns = vars(MO), align = "left")

# Almacenamiento de tablas

FIM_RESPTA_1 %>% 
  gtsave(filename = '463_CFR_Indicador_100fTmasMIC.html', path = file.path(getwd(), 'figures'))

FIM_RESPTA_2 %>% 
  gtsave(filename = '464_CFR_Indicador_50fTmas4MIC.html', path = file.path(getwd(), 'figures'))

FIM_RESPTA_60 %>% 
  gtsave(filename = '468_CFR_Indicador_60fTmasMIC.html', path = file.path(getwd(), 'figures'))


  
  
  
  
  
  
  