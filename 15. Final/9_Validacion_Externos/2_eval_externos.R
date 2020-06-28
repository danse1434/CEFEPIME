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
require(boot)
require(gt)
require(tidyverse)

#-------------------------------------------------------------------------------#
# Introducción -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Lectura de archivos con PRED vs OBS
pred_list <- list()

pred_list[[1]] <-
  read_csv(
    '1_Rhodes/2_Rhodes_TAD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )
pred_list[[2]] <-
  read_csv(
    '2_Whited/2_Whited_TAD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )
pred_list[[3]] <-
  read_csv(
    '3_Lee/3_Lee_TAD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )
pred_list[[4]] <-
  read_csv(
    '4_Sime/4_Sime_TAD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )
pred_list[[5]] <-
  read_csv(
    '5_Roos_UCI/5_Roos_UCI_TAD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )
pred_list[[6]] <-
  read_csv(
    '6_Delattre_UCI/6_Delattre_UCI_TAD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )
pred_list[[7]] <-
  read_csv(
    '7_Georges_UCI/7_Georges_UCI_TAD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )
pred_list[[8]] <-
  read_csv(
    '8_Nicasio_UCI/8_Nicasio_UCI_TAD/ChartsData/ObservationsVsPredictions/y_1_obsVsPred.txt'
  )


#-------------------------------------------------------------------------------#
# Lista de modelos probados
model_list <- tribble(
  ~Modelo, ~Nombre,
  '1', "Rhodes et al (NF)",
  '2', "Whited et al (NF)",
  '3', "Lee et al (NF)",
  '4', "Sime et al (NF)",
  '5', "Roos et al (UCI)",
  '6', "Delattre et al (UCI)",
  '7', "Georges et al (UCI)",
  '8', "Nicasio et al (UCI)"
  )

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
  mutate(Modelo = factor(Modelo, paste0(1:8))) %>% 
  arrange(Modelo)
#-------------------------------------------------------------------------------#
# Resumen de parámetros
#' Calcular resumen de parámetros de predicción
#'  funciones para calcular bootstrap
MPE_b <- function(data, indices) {
  d <- data[indices, ]
  m <- mean(d$PE)
  return(m)
}

RMSE_b <- function(data, indices) {
  d <- data[indices, ]
  s <- sqrt(sum(d$PE^2))
  return(s)
}

#-------------------------------------------------------------------------------#
#' Manipulación de tabla *pred_df*
#................................................................................
#' 1 Agrupar por _Nombre_
#' 2 Anidar por _Nombre_
#' 3 Realizar un bootstrap no paramétrico con una distribución de *MPE_b* y 
#' *RMSE_b*
#' 4 Obtención de intervalos de confianza de percentil del 95%
#' 5 Calculo de tamaño muestral
#' 6 Obtención de _MPE_ y _RMSE_
#' 7 Desanidar de MPE_CI95, y RMSE_CI95
#' 8 Adicionar columna con límites inferior y superior
#' 9 Convertir filas de CI95% a columna
#................................................................................

pred_df1 <- pred_df %>%
  group_by(Modelo, Nombre) %>%
  nest() %>% 
  mutate(
    MPE_bt    = map(data, ~boot(data = .x, statistic = MPE_b, R = 1E3)),
    RMSE_bt   = map(data, ~boot(data = .x, statistic = RMSE_b, R = 1E3)),
    MPE_CI95  = map(MPE_bt,  ~boot.ci(.x, type = 'perc')$perc[4:5]),
    RMSE_CI95 = map(RMSE_bt, ~boot.ci(.x, type = 'perc')$perc[4:5]),
    n         = map_dbl(data, ~ length(.x$PE)),
    MPE       = map_dbl(data, ~ mean(.x$PE)),
    RMSE      = map_dbl(data, ~ sqrt(sum(.x$PE ^ 2)))) %>%
  unnest(c(MPE_CI95, RMSE_CI95)) %>% 
  add_column(L = rep(c('li', 'ls'), 8)) %>%
  pivot_wider(names_from = L, values_from = c(MPE_CI95, RMSE_CI95))

#-------------------------------------------------------------------------------#
# Creación de tabla resultados GOF
nota_pie <- glue::glue("
[1] Whited L, Grove M, Rose D, et al. Pharmacokinetics of Cefepime in Patients with Cancer and Febrile Neutropenia in the Setting of Hematologic Malignancies or Hematopoeitic Cell Transplantation. Pharmacother J Hum Pharmacol Drug Ther 2016; 36: 1003–1010.
[2] Rhodes NJ, Grove ME, Kiel PJ, et al. Population pharmacokinetics of cefepime in febrile neutropenia: implications for dose-dependent susceptibility and contemporary dosing regimens. Int J Antimicrob Agents 2017; 50: 482–486.
[3] Lee DG, Choi SM, Yoo JH, et al. Population pharmacokinetics of cefepime in febrile neutropenic patients. J Korean Soc Clin Pharmacol Ther 2003; 11: 23–29.
[4] Sime FB, Roberts MS, Tiong IS, et al. Adequacy of High-Dose Cefepime Regimen in Febrile Neutropenic Patients with Hematological Malignancies. Antimicrob Agents Chemother 2015; 59: 5463–5469.
[5] Roos JF, Bulitta J, Lipman J, et al. Pharmacokinetic-pharmacodynamic rationale for cefepime dosing regimens in intensive care units. J Antimicrob Chemother 2006; 58: 987–993.
[6] Delattre IK, Musuamba FT, Jacqmin P, et al. Population pharmacokinetics of four β-lactams in critically ill septic patients comedicated with amikacin. Clin Biochem 2012; 45: 780–786.
[7] Georges B, Conil J-M, Seguin T, et al. Cefepime in intensive care unit patients: Validation of a population pharmacokinetic approach and influence of covariables. Int J Clin Pharmacol Ther 2008; 46: 157–164.
[8] Nicasio AM, Ariano RE, Zelenitsky SA, et al. Population pharmacokinetics of high-dose, prolonged-infusion cefepime in adult critically 111 patients with ventilator-associated pneumonia. Antimicrob Agents Chemother 2009; 53: 1476–1481.")

pred_tbl <- pred_df1 %>% 
  ungroup() %>% 
  select(-Modelo, -data, -MPE_bt, -RMSE_bt) %>% 
  gt() %>% 
  fmt_number(columns = 3:8, decimals = 3) %>% 
  tab_header(
    title = md('**Resultados de validación externa de modelos de referencia**')
  ) %>% 
  # cols_move_to_end(columns = vars(RMSE)) %>% 
  cols_merge(columns = vars(MPE_CI95_li, MPE_CI95_ls), pattern = "({1}, {2})") %>% 
  cols_merge(columns = vars(MPE, MPE_CI95_li), pattern = "{1} {2}") %>%
  
  cols_merge(columns = vars(RMSE_CI95_li, RMSE_CI95_ls), pattern = "({1}, {2})") %>% 
  cols_merge(columns = vars(RMSE, RMSE_CI95_li), pattern = "{1} {2}") %>% 
  cols_label(
    Nombre = md("**Estudio**"),
    MPE    = md("**MPE [IC95%]**"),
    RMSE   = md("**RMSE [IC95%]**"),
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
    container.width = pct(100))

gtsave(pred_tbl,
       "4_tabla_resultados.html",
       file.path(getwd(), "resultados"), inline_css = TRUE)
  
#-------------------------------------------------------------------------------#
# Gráfico Bland-Altman -----------------------------------------------------
#-------------------------------------------------------------------------------#
QUA_b <- function(data, indices) {
  d <- data[indices, ]
  rang <- quantile(d$PE, probs = c(0.05, 0.50, 0.95))
  return(rang)
}

pred_df2 <- pred_df %>%
  group_by(Modelo, Nombre) %>%
  nest() %>% 
  mutate(
    QUA_bt    = map(data, ~boot(data = .x, statistic = QUA_b, R = 1E3)),
    QUA_P1    = map_dbl(QUA_bt,  ~boot.ci(.x, type='perc', index=1)$t0[[1]]),
    QUA_P2    = map_dbl(QUA_bt,  ~boot.ci(.x, type='perc', index=2)$t0[[1]]),
    QUA_P3    = map_dbl(QUA_bt,  ~boot.ci(.x, type='perc', index=3)$t0[[1]]),
    QUA_P1_CI = map(QUA_bt,  ~boot.ci(.x, type='perc', index=1)$perc[4:5]),
    QUA_P2_CI = map(QUA_bt,  ~boot.ci(.x, type='perc', index=2)$perc[4:5]),
    QUA_P3_CI = map(QUA_bt,  ~boot.ci(.x, type='perc', index=3)$perc[4:5])) %>%
  unnest(c(QUA_P1_CI, QUA_P2_CI, QUA_P3_CI)) %>% 
  add_column(L = rep(c('li', 'ls'), 8)) %>%
  pivot_wider(names_from = L, values_from = c(QUA_P1_CI, QUA_P2_CI, QUA_P3_CI))

G_BA <- pred_df %>% 
  mutate(Nombre = fct_reorder(Nombre, as.integer(Modelo), .desc = T)) %>% 
  ggplot(aes(x = PR, y = PE * 100, col = Nombre)) +
  geom_point() +
  theme_classic() +
  geom_hline(yintercept = 0) +
  geom_hline(data = pred_df2, aes(yintercept = QUA_P1 * 100), 
             col = "red4", size = 0.5, lty = 'dashed') +
  geom_hline(data = pred_df2, aes(yintercept = QUA_P2 * 100), 
             col = "red4", size = 0.5, lty = 'dashed') +
  geom_hline(data = pred_df2, aes(yintercept = QUA_P3 * 100), 
             col = "red4", size = 0.5, lty = 'dashed') +
  
  geom_rect(data = pred_df2, 
              aes(ymin = QUA_P1_CI_li * 100, ymax = QUA_P1_CI_ls * 100, 
                  xmin = 0, xmax = 120), 
              fill = alpha("red", 0.25), inherit.aes = FALSE) +
  geom_rect(data = pred_df2, 
            aes(ymin = QUA_P2_CI_li * 100, ymax = QUA_P2_CI_ls * 100, 
                xmin = 0, xmax = 120), 
            fill = alpha("red", 0.25), inherit.aes = FALSE) +
  geom_rect(data = pred_df2, 
            aes(ymin = QUA_P3_CI_li * 100, ymax = QUA_P3_CI_ls * 100, 
                xmin = 0, xmax = 120), 
            fill = alpha("red", 0.25), inherit.aes = FALSE) +
  
  facet_wrap(~Nombre) +
  labs(title = 'Gráfico de Bland-Altman para modelos de Referencia', 
       subtitle = 'Predicciones con TSFD vs TAD') +
  xlab('Promedio de Pred y Obs (mg/L)') + 
  ylab('Error de predicción (%)') +
  coord_cartesian(xlim = c(0,120), ylim = c(-100, 600), expand = FALSE) +
  theme(panel.border = element_rect(fill = NA, colour = 'black'),
        legend.position = "none") 

G_BA  

ggsave('5_graf_Bland-Altman.pdf', G_BA, 'pdf', 'resultados', 1, 7, 5)

#-------------------------------------------------------------------------------#
# Gráfico de predicciones -----------------------------------------------------
#-------------------------------------------------------------------------------#
G_PRED <- pred_df %>% 
  mutate(Nombre = fct_reorder(Nombre, as.integer(Modelo), .desc = T)) %>% 
  ggplot(aes(x = PE * 100, y = Nombre)) +
  geom_boxplot(fill = alpha('gray', 0.4)) +
  theme_bw() + 
  geom_vline(xintercept = c(-20, +20), lty = 'dashed', col = 'red3') +
  coord_cartesian(xlim = c(-100, 600)) +
  xlab('Error de predicción (%)') +
  theme(axis.title.y = element_blank())
G_PRED

ggsave('6_graf_predicciones.pdf', G_PRED, 'pdf', 'resultados', 1, 6.5, 4)

