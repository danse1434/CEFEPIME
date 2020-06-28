##------------------------------------------------------------------------------#
## Nombre del Script: Validación externa de modelos de referencia --------------
##  
## Propósito del Script: en este script se calculan varios parámetros que muestran 
## la capacidad de modelos externos para explicar los datos de este estudio.
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  21-06-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
# Carga de paquetes
require(gt)
require(tidyverse)

#-------------------------------------------------------------------------------#
# Introducción -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Lectura de archivos con PRED vs OBS
pred_list <- list()
pred_list[[1]] <-
  read_csv(
    '1_Rhodes/1_Rhodes_TSFD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )
pred_list[[2]] <-
  read_csv(
    '2_Whited/1_Whited_TSFD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )
pred_list[[3]] <-
  read_csv(
    '1_Rhodes/2_Rhodes_TAD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )
pred_list[[4]] <-
  read_csv(
    '2_Whited/2_Whited_TAD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )

#-------------------------------------------------------------------------------#
# Lista de modelos probados
model_list <- tribble(
  ~Modelo, ~Nombre,
    '1', "Rhodes et al (TSFD)",
    '2', "Whited et al (TSFD)",
    '3', "Rhodes et al (TAD)",
    '4', "Whited et al (TAD)")

#-------------------------------------------------------------------------------#
# Análisis de parámetros GOF  ---------------------------------------------------
#-------------------------------------------------------------------------------#
#' DataFrame con análisis de predicción
#................................................................................
#' 1 Convertir lista a data.frame
#' 2 Seleccionar columnas específicas
#' 3 Calcular error de predicción _PE_ y promedio de predicción _PR_
#' 4 Unir lista de modelos *model_list*
#' 5 Ordenar por _Nombre_
#................................................................................

pred_df <- pred_list %>% 
  map_dfr(~.x, .id = 'Modelo') %>% 
  select(Modelo, ID, time, y_1, popPred) %>% 
  mutate(PE = (popPred - y_1)/y_1,
         PR = (popPred + y_1)/2) %>% 
  left_join(model_list, by = 'Modelo') %>% 
  arrange(Nombre)

#-------------------------------------------------------------------------------#
# Resumen de parámetros
#' Calcular resumen de parámetros de predicción
#................................................................................
#' Agrupar por _Nombre_
#' Resumir por _n_, error medio de predicción _MPE_, y error estándar medio 
#' cuadrático _RMSE_
#................................................................................
pred_df1 <- pred_df %>%
  group_by(Nombre) %>%
  summarise(
    n = n(),
    MPE = mean(PE),
    RMSE = sqrt(sum(PE ^ 2) / n()),
    MPE_li = MPE - (1.96 * sd(PE)),
    MPE_ls = MPE + (1.96 * sd(PE))
  )

#-------------------------------------------------------------------------------#
# Creación de tabla resultados GOF
nota_pie <- glue::glue("
[1] Whited L, Grove M, Rose D, et al. Pharmacokinetics of Cefepime in Patients with Cancer and Febrile Neutropenia in the Setting of Hematologic Malignancies or Hematopoeitic Cell Transplantation. Pharmacother J Hum Pharmacol Drug Ther 2016; 36: 1003–1010.
[2] Rhodes NJ, Grove ME, Kiel PJ, et al. Population pharmacokinetics of cefepime in febrile neutropenia: implications for dose-dependent susceptibility and contemporary dosing regimens. Int J Antimicrob Agents 2017; 50: 482–486.")

pred_tbl <- pred_df1 %>% 
  gt() %>% 
  fmt_number(columns = 3:6, decimals = 3) %>% 
  tab_header(
    title = md('**Resultados de validación externa de modelos de referencia**')
  ) %>% 
  cols_move_to_end(columns = vars(RMSE)) %>% 
  cols_merge(columns = vars(MPE_li, MPE_ls), pattern = "[{1}, {2}]") %>% 
  cols_merge(columns = vars(MPE, MPE_li), pattern = "{1} {2}") %>% 
  cols_label(
    Nombre = md("**Estudio**"),
    MPE    = md("**MPE [IC95%]**"),
    n      = md("**N**"),
    RMSE   = md("**RMSE**")
  ) %>% 
  tab_footnote(
    locations = cells_column_labels(columns = vars(Nombre)),
    footnote = nota_pie
  ) %>% 
  opt_row_striping(row_striping = TRUE) %>% 
  tab_options(
    heading.title.font.size = px(18),
    footnotes.font.size = px(10), 
    container.width = pct(100), 
    container.overflow.x = TRUE)

gtsave(pred_tbl,
       "1_tabla_resultados.html",
       file.path(getwd(), "resultados"), inline_css = TRUE)
  
#-------------------------------------------------------------------------------#
# Gráfico Bland-Altman -----------------------------------------------------
#-------------------------------------------------------------------------------#
G_BA <- pred_df %>% 
  ggplot(aes(x = PR, y = PE * 100, col = Nombre)) +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept = 0) +
  geom_hline(data = pred_df1, aes(yintercept = MPE * 100), 
             col = "red", size = 0.5, lty = 'dashed') +
  geom_hline(data = pred_df1, aes(yintercept = MPE_li * 100), 
             col = "red", size = 0.5, lty = 'dashed') +
  geom_hline(data = pred_df1, aes(yintercept = MPE_ls * 100), 
             col = "red", size = 0.5, lty = 'dashed') +
  facet_wrap(~Nombre) +
  labs(title = 'Gráfico de Bland-Altman para modelos de Referencia', 
       subtitle = 'Predicciones con TSFD vs TAD') +
  xlab('Promedio de Pred y Obs (mg/L)') + 
  ylab('Error de predicción (%)') +
  coord_cartesian(ylim = c(-400, 600)) +
  theme(panel.border = element_rect(fill = NA, colour = 'black'),
        legend.position = "none") 
  
ggsave('2_graf_Bland-Altman.pdf', G_BA, 'pdf', 'resultados', 1, 7, 5)

#-------------------------------------------------------------------------------#
# Gráfico de predicciones -----------------------------------------------------
#-------------------------------------------------------------------------------#
G_PRED <- pred_df %>% 
  ggplot(aes(x = PE * 100, y = Nombre)) +
  geom_boxplot(fill = alpha('gray', 0.4)) +
  theme_bw() + 
  geom_vline(xintercept = c(-20, +20), lty = 'dashed', col = 'red3') +
  coord_cartesian(xlim = c(-400, 600)) +
  xlab('Error de predicción (%)') +
  theme(axis.title.y = element_blank())

ggsave('3_graf_predicciones.pdf', G_PRED, 'pdf', 'resultados', 1, 6.5, 4)
