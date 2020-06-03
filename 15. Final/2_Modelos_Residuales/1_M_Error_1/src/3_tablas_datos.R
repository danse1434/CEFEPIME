##------------------------------------------------------------------------------#
## Nombre del Script: Lectura de tablas de resultados de Monolix ---
##  
## Propósito del Script: Leer tablas de datos de resultados de Monolix
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 04-feb-2020  
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#

# Carga de paquetes
require(gt)
require(tidyverse)

# Especificación de carpeta del modelo
auxdir <- file.path('BASE_MODEL')

#-------------------------------------------------------------------------------#
# Matriz FIM -----------------------------------------------------
#-------------------------------------------------------------------------------#
# Matriz de correlación

FIM1 <-
  read_csv(file.path(auxdir, 'FisherInformation', 'correlationEstimatesSA.txt'), 
           col_names = FALSE)


colnames(FIM1) = c("Parametro", FIM1$X1)

#-------------------------------------------------------------------------------#
fun_param <- function(data, param) {
  param_quo <- rlang::ensym(param)

  z = data %>% 
    tab_style(style = list(cell_fill(color = 'yellow4'),
                           cell_text(weight = 'bold')), 
    locations = cells_body(columns = vars(!!param_quo), 
                           rows = (!!param_quo >= 0.5 & !!param_quo < 0.8))
    ) %>% 
    tab_style(style = list(cell_fill(color = 'orange'),
                           cell_text(weight = 'bold')), 
    locations = cells_body(columns = vars(!!param_quo), 
                           rows = (!!param_quo >= 0.8 & !!param_quo < 0.9))
    ) %>% 
    tab_style(style = list(cell_fill(color = 'red4'),
                           cell_text(weight = 'bold')), 
    locations = cells_body(columns = vars(!!param_quo), 
                             rows = (!!param_quo >= 0.9))
    )
  
  return(z)
}

FIMGT1 <- gt(FIM1) %>%
  tab_header(title = 'Matriz de correlación',
             subtitle = 'Modelo Final de Covariables') %>%
  fmt_number(columns = 2:10, decimals = 3) %>%
  fun_param('Cl_pop') %>%
  fun_param('beta_Cl_tSCRMGDL') %>%
  fun_param('V1_pop') %>%
  fun_param('V2_pop') %>%
  fun_param('omega_Cl') %>%
  fun_param('omega_V1') %>%
  fun_param('omega_Q') %>%
  fun_param('omega_V2') %>%
  fun_param('a') #%>%
  # fun_param('b')

#-------------------------------------------------------------------------------#
# Matriz de COVARIANZA

FIM2 <-
  read_csv(file.path(auxdir, 'FisherInformation', 'covarianceEstimatesSA.txt'), 
           col_names = FALSE)

colnames(FIM2) = c("Parametro", FIM2$X1)

FIMGT2 <- gt(FIM2) %>% 
  tab_header(title = 'Matriz de covarianzas', 
             subtitle = 'Modelo Final de Covariables') %>% 
  fmt_number(columns = 2:10, decimals = 3) 


#-------------------------------------------------------------------------------#
# Pruebas estadísticas  ---------------------------------------------------------
#-------------------------------------------------------------------------------#

# Prueba de correlación a parámetros del modelo de covariables
testdf1 <- read_csv(file.path(auxdir, 'Tests', 
                   'correlationIndividualParametersCovariates.txt'))

testGT1 <- testdf1 %>% 
  gt() %>% 
  tab_header(title = 'Prueba de correlación de Pearson', 
             subtitle = 'Modelo Final de Covariables') %>% 
  cols_label(parameter = md('*Parámetro*'),
             covariate = 'Covariable', 
             statistics = 'Estadístico', 
             `p-value` = 'Valor p', 
             test = 'Prueba')

# Se considera que el valor del parámetro tSCRMGDL no es significativamente 
# diferente de cero a un nivel de significancia de 0.05.

#-------------------------------------------------------------------------------#
# Prueba de correlación efectos aleatorios en distribución conjunta de covariables

testdf2 <- read_csv(file.path(auxdir, 'Tests', 
                              'correlationRandomEffects.txt'))
testGT2 <- testdf2 %>% 
  gt() %>% 
  cols_label(
    eta1       = md("**Eta**"),
    eta2       = md(""),
    coeff      = 'Coeficiente',
    statistics = "Estadístico",
    `p-value`  = 'Valor p') %>%
  tab_header(title = 'Prueba de correlación entre distribución conjunta de efectos aleatorios',
             subtitle = 'Modelo Final de Covariables') %>% 
  fmt_number(vars(statistics,`p-value`), decimals = 4) %>%
  tab_style(
    style = list(cell_fill(color = '#F0EDE6'),
                 cell_text(weight = 'bold')), 
    locations = cells_body(columns = vars(`p-value`), 
                           rows = (`p-value` < 0.05))
  ) %>% 
  cols_hide(vars(coeff))


#-------------------------------------------------------------------------------#
# Efectos aleatorios vs covariables 
# Esto sirve para evaluar si se deben adicionar estas covariables 

testdf3 <- read_csv(file.path(auxdir, 'Tests', 
                              'correlationRandomEffectsCovariates.txt'))

testGT3 <- testdf3 %>% 
  gt(groupname_col = 'covariate') %>% 
  cols_label(eta = "Eta", 
             covariate = 'Covariable', statistics = 'Estadístico', 
             `p-value` = 'Valor p', test = 'Prueba') %>% 
  tab_header(title = 'Efectos aleatorios vs covariables',
             subtitle = 'Modelo Final de Covariables') %>% 
  fmt_number(columns = vars(statistics, `p-value`), decimals = 3) 

#-------------------------------------------------------------------------------#
# Test de Wald

testdf4 <- read_csv(file.path(auxdir, 'Tests', 'fixedEffects.txt'))

testGT4 <- testdf4 %>% 
  gt() %>% 
  tab_header(title    = 'Test de Wald', 
             subtitle = 'Modelo Final de Covariables') %>% 
  cols_label(parameter        = md('*Parámetro*'),
             `statistics(sa)` = 'Estadístico', 
             `p-value(sa)`    = 'Valor p') %>% 
  tab_style(
    style = list(cell_fill(color = '#F0EDE6'),
                 cell_text(weight = 'bold')), 
    locations = cells_body(columns = vars(`p-value(sa)`), 
                           rows = (`p-value(sa)` < 0.05))
  )

#-------------------------------------------------------------------------------#
# Prueba de normalidad en parámetros farmacocinéticos

testdf5 <- read_csv(file.path(auxdir, 'Tests', 
                              'normalityIndividualParameters.txt'))
testGT5 <- testdf5 %>% 
  gt() %>% 
  cols_label(parameter = 'Parámetro', 
             statistics = 'Estadístico', 
             `p-value` = 'Valor p', 
             test = 'Prueba') %>% 
  tab_header(title = 'Prueba de normalidad de parámetros PK',
             subtitle = 'Modelo Final de Covariables') %>% 
  tab_footnote('Esta prueba se realiza en el dominio transformado',
               cells_column_labels(columns = vars(test))) %>% 
  fmt_number(vars(statistics,`p-value`), decimals = 3)

#-------------------------------------------------------------------------------#
# Prueba de normalidad en efectos aleatorios
# Supuesto importante en modelos PK

testdf6 <- read_csv(file.path(auxdir, 'Tests', 
                              'normalityRandomEffects.txt'))

testGT6 <- gt(testdf6) %>% 
  cols_label(eta        = md('**Eta**'), 
             statistics = md('**Estadístico**'), 
             `p-value`  = md('**Valor p**')) %>% 
  tab_header(title = 'Prueba de normalidad en efectos aleatorios',
             subtitle = 'Modelo Final de Covariables') %>% 
  fmt_number(vars(statistics,`p-value`), decimals = 4) 


#-------------------------------------------------------------------------------#
# Prueba de normalidad en residuales

testdf7 <- read_csv(file.path(auxdir, 'Tests', 'normalityResiduals.txt'))

testGT7 <- testdf7 %>% 
  mutate(residuals = str_replace(residuals, '\\_y\\_1', '')) %>% 
  gt() %>% 
  cols_label(residuals  = 'Residuales', 
             statistics = 'Estadístico', 
             `p-value`  = 'Valor p', 
             test       = 'Prueba') %>% 
  tab_header(title    = 'Prueba de normalidad en residuales',
             subtitle = 'Modelo Final de Covariables') %>% 
  fmt_number(vars(statistics,`p-value`), decimals = 4) %>% 
  tab_style(
    style = list(cell_fill(color = '#F0EDE6'),
                 cell_text(weight = 'bold')), 
    locations = cells_body(columns = vars(`p-value`), 
                           rows = (`p-value` < 0.05))
  )

#-------------------------------------------------------------------------------#
# Pruebas de simetría en residuales

testdf8 <- read_csv(file.path(auxdir, 'Tests', 'symmetryResiduals.txt'))

testGT8 <- testdf8 %>% 
  mutate(test = 'Simetría', 
         residuals = str_replace(residuals, '\\_y\\_1', '')) %>% 
  gt() %>% 
  cols_label(residuals  = 'Residuales', 
             statistics = 'Estadístico', 
             `p-value`  = 'Valor p', 
             test       = 'Prueba') %>% 
  tab_header(title    = 'Pruebas de simetría en residuales',
             subtitle = 'Modelo Final de Covariables') %>% 
  fmt_number(vars(statistics,`p-value`), decimals = 4) %>% 
  tab_style(
    style = list(cell_fill(color = '#F0EDE6'),
                 cell_text(weight = 'bold')), 
    locations = cells_body(columns = vars(`p-value`), 
                           rows = (`p-value` < 0.05))
    )

#-------------------------------------------------------------------------------#
# Almacenamiento de tablas en formato HTML

testGT1 %>% 
  gtsave(filename = 'testGT1.html', path = file.path(getwd(), 'tables'))

testGT2 %>% 
  gtsave(filename = 'testGT2.html', path = file.path(getwd(), 'tables'))

testGT3 %>% 
  gtsave(filename = 'testGT3.html', path = file.path(getwd(), 'tables'))

testGT4 %>% 
  gtsave(filename = 'testGT4.html', path = file.path(getwd(), 'tables'))

testGT5 %>% 
  gtsave(filename = 'testGT5.html', path = file.path(getwd(), 'tables'))

testGT6 %>% 
  gtsave(filename = 'testGT6.html', path = file.path(getwd(), 'tables'))

testGT7 %>% 
  gtsave(filename = 'testGT7.html', path = file.path(getwd(), 'tables'))

testGT8 %>% 
  gtsave(filename = 'testGT8.html', path = file.path(getwd(), 'tables'))

