##------------------------------------------------------------------------#
## Nombre del Script: Anólisis de distribución de errores residuales ---
##  
## Proposito del Script: Anólisis de sesgo con mótodos de censura de datos
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion:  06-feb-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
# Carga de paquetes
require(tidyverse)
require(rlang)
require(grid)
require(gridExtra)
require(readxl)
##########################################################################-
# Selección de directorio de trabajo
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '04. Residual'))

# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

##########################################################################-
# Grófico de residuales vs TAD --------------------------------------------
##########################################################################-
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Apertura de archivo de datos con los datos originales, y que tiene 
##  delimitador ";", este fue usado para realizar el modelamiento en só mismo
##  2 Crear un archivo *last_dose* que contiene la óltima hora de adminis-
##  tración del antibiótico en el archivo *data* para cada individuo. 
##    2a Agrupar por ID
##    2b Filtrar los datos con EVID igual a 1 (sólo los eventos de 
##    administración)
##    2c Seleccionar los óltimos en cada grupo (los datos ya estón ordenados 
##    de menor a mayor)
##    2d Seleccionar sólo las columnas ID y TIME
##  3 Modificar el archivo de residuales *y_1_residuals* al unir los datos 
##  de óltimnas observaciones calculadas.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data <- read_delim("Modelo_1/interv_censored.csv",
                   delim = ";", 
                   escape_double = FALSE,
                   trim_ws = TRUE)

last_dose <- data %>% 
  group_by(ID) %>% 
  filter(EVID == 1) %>% 
  slice(n()) %>% 
  select(ID, TIME)

##########################################################################-
# Apertura de archivos de residuales

res_1 <- read_csv(file.path('Modelo_1', 'RES_M1', 'ChartsData', 
                            'ScatterPlotOfTheResiduals','y_1_residuals.txt'))

res_2 <- read_csv(file.path('Modelo_2', 'RES_M1', 'ChartsData', 
                            'ScatterPlotOfTheResiduals','y_1_residuals.txt'))

res_3 <- read_csv(file.path('Modelo_3', 'RES_M1', 'ChartsData', 
                            'ScatterPlotOfTheResiduals','y_1_residuals.txt'))

res_4 <- read_csv(file.path('Modelo_4', 'RES_M1', 'ChartsData', 
                            'ScatterPlotOfTheResiduals','y_1_residuals.txt'))

res_5 <- read_csv(file.path('Modelo_5', 'RES_M1', 'ChartsData', 
                            'ScatterPlotOfTheResiduals','y_1_residuals.txt'))

res_6 <- read_csv(file.path('Modelo_6', 'RES_M1', 'ChartsData', 
                            'ScatterPlotOfTheResiduals','y_1_residuals.txt'))

res_7 <- read_csv(file.path('Modelo_7', 'RES_M1', 'ChartsData', 
                            'ScatterPlotOfTheResiduals','y_1_residuals.txt'))


# Adicionar columna TAD
TAD_add <- function(data) {
  data = data %>%
    left_join(., last_dose, by = 'ID') %>%
    mutate(TAD = time - TIME)
  return(data)
}

res_1 <- TAD_add(res_1)
res_2 <- TAD_add(res_2)
res_3 <- TAD_add(res_3)
res_4 <- TAD_add(res_4)
res_5 <- TAD_add(res_5)
res_6 <- TAD_add(res_6)
res_7 <- TAD_add(res_7)

##########################################################################-


RES_C <- function(data, x, y, xlab, ylab) {
  ##########################################################################-
  # Volver las variables x y y en expresiones para evaluación tardóa
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ##  Calculo de percentiles empiricos
  ##  1 Tomar la variable de datos a explorar.
  ##  2 Crear una variable que tome 6 grupos a partir de los datos ordenados 
  ##  con la función dplyr::ntile().
  ##  3 Agrupar por la variable reción creado.
  ##  4 Resumir por la media de TAD, y los percentiles P5%, P50%, y P95% de 
  ##  los datos de residuales.
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  A <- data %>%
    mutate(gr = ntile(TAD, 6)) %>% 
    group_by(gr) %>% 
    summarise(TIME = mean(TAD), 
              ME = quantile(x = !!y, probs = 0.50), 
              LI = quantile(x = !!y, probs = 0.05), 
              LS = quantile(x = !!y, probs = 0.95))
  ##########################################################################-
  # Crear el grófico de residuales con las especificaciones mostradas
  data %>% 
    ggplot(mapping = aes(x = !!x, y = !!y)) +
    geom_hline(yintercept = 0) +
    geom_point(col = '#4682B4') +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    xlab(xlab) + ylab(ylab) +
    stat_smooth(method = 'loess', formula = y ~ x, se = FALSE, 
                col = '#EBA213', lty = 'solid', size = 1) +
    geom_line(data = A, aes(TIME, ME), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(data = A, aes(TIME, LI), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(data = A, aes(TIME, LS), 
              col = '#EBA213', lty = 'dashed', size = 1.2) %>% 
    return(.)
}

G_RES_C_PWRES <- list()

G_RES_C_PWRES[[1]] <- RES_C(res_1, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')
G_RES_C_PWRES[[2]] <- RES_C(res_2, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')
G_RES_C_PWRES[[3]] <- RES_C(res_3, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')
G_RES_C_PWRES[[4]] <- RES_C(res_4, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')
G_RES_C_PWRES[[5]] <- RES_C(res_5, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')
G_RES_C_PWRES[[6]] <- RES_C(res_6, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')
G_RES_C_PWRES[[7]] <- RES_C(res_7, x = TAD, y = pwRes, xlab = 'TAD', ylab = 'PWRES')

G_RES_C_IWRES <- list()
G_RES_C_IWRES[[1]] <- RES_C(res_1, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')
G_RES_C_IWRES[[2]] <- RES_C(res_2, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')
G_RES_C_IWRES[[3]] <- RES_C(res_3, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')
G_RES_C_IWRES[[4]] <- RES_C(res_4, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')
G_RES_C_IWRES[[5]] <- RES_C(res_5, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')
G_RES_C_IWRES[[6]] <- RES_C(res_6, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')
G_RES_C_IWRES[[7]] <- RES_C(res_7, x = TAD, y = iwRes_mean, xlab = 'TAD', ylab = 'IWRES')

G_RES_C_NPDE <- list()
G_RES_C_NPDE[[1]] <- RES_C(res_1, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')
G_RES_C_NPDE[[2]] <- RES_C(res_2, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')
G_RES_C_NPDE[[3]] <- RES_C(res_3, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')
G_RES_C_NPDE[[4]] <- RES_C(res_4, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')
G_RES_C_NPDE[[5]] <- RES_C(res_5, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')
G_RES_C_NPDE[[6]] <- RES_C(res_6, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')
G_RES_C_NPDE[[7]] <- RES_C(res_7, x = TAD, y = npde, xlab = 'TAD', ylab = 'NPDE')

pdf(file = './SCRIPT/PWRES_C.pdf', width = 6, height = 7);{
gridExtra::grid.arrange(G_RES_C_PWRES[[1]], G_RES_C_PWRES[[2]], G_RES_C_PWRES[[3]], 
                        G_RES_C_PWRES[[4]], G_RES_C_PWRES[[5]], G_RES_C_PWRES[[6]],
                        G_RES_C_PWRES[[7]])
};dev.off()

pdf(file = './SCRIPT/IWRES_C.pdf', width = 6, height = 7);{
gridExtra::grid.arrange(G_RES_C_IWRES[[1]], G_RES_C_IWRES[[2]], G_RES_C_IWRES[[3]], 
                        G_RES_C_IWRES[[4]], G_RES_C_IWRES[[5]], G_RES_C_IWRES[[6]],
                        G_RES_C_IWRES[[7]])
};dev.off()

pdf(file = './SCRIPT/NPDE_C.pdf', width = 6, height = 7);{
gridExtra::grid.arrange(G_RES_C_NPDE[[1]], G_RES_C_NPDE[[2]], G_RES_C_NPDE[[3]], 
                        G_RES_C_NPDE[[4]], G_RES_C_NPDE[[5]], G_RES_C_NPDE[[6]],
                        G_RES_C_NPDE[[7]])
};dev.off()








RES_PRE <- function(data, x, y, xlab, ylab) {
  ##########################################################################-
  # Volver las variables x y y en expresiones para evaluación tardóa
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  ##  Calculo de percentiles empiricos
  ##  1 Tomar la variable de datos a explorar.
  ##  2 Crear una variable que tome 6 grupos a partir de los datos ordenados 
  ##  con la función dplyr::ntile().
  ##  3 Agrupar por la variable reción creado.
  ##  4 Resumir por la media de TAD, y los percentiles P5%, P50%, y P95% de 
  ##  los datos de residuales.
  ##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  A <- data %>%
    mutate(gr = ntile(!!x, 6)) %>% 
    group_by(gr) %>% 
    summarise(PRED = mean(!!x), 
              ME = quantile(x = !!y, probs = 0.50), 
              LI = quantile(x = !!y, probs = 0.05), 
              LS = quantile(x = !!y, probs = 0.95))
  ##########################################################################-
  # Crear el grófico de residuales con las especificaciones mostradas
  data %>% 
    ggplot(mapping = aes(x = !!x, y = !!y)) +
    geom_hline(yintercept = 0) +
    geom_point(col = '#4682B4') +
    coord_cartesian(ylim = c(-2.5, 2.5)) +
    xlab(xlab) + ylab(ylab) +
    stat_smooth(method = 'loess', formula = y ~ x, se = FALSE, 
                col = '#EBA213', lty = 'solid', size = 1) +
    geom_line(data = A, aes(PRED, ME), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(data = A, aes(PRED, LI), 
              col = '#EBA213', lty = 'dashed', size = 1.2) +
    geom_line(data = A, aes(PRED, LS), 
              col = '#EBA213', lty = 'dashed', size = 1.2) %>% 
    return(.)
}

G_RES_C_PWRES <- list()
G_RES_C_PWRES[[1]] <- RES_PRE(res_1, x = prediction_pwRes, y = pwRes, xlab = 'PRED', ylab = 'PWRES')
G_RES_C_PWRES[[2]] <- RES_PRE(res_2, x = prediction_pwRes, y = pwRes, xlab = 'PRED', ylab = 'PWRES')
G_RES_C_PWRES[[3]] <- RES_PRE(res_3, x = prediction_pwRes, y = pwRes, xlab = 'PRED', ylab = 'PWRES')
G_RES_C_PWRES[[4]] <- RES_PRE(res_4, x = prediction_pwRes, y = pwRes, xlab = 'PRED', ylab = 'PWRES')
G_RES_C_PWRES[[5]] <- RES_PRE(res_5, x = prediction_pwRes, y = pwRes, xlab = 'PRED', ylab = 'PWRES')
G_RES_C_PWRES[[6]] <- RES_PRE(res_6, x = prediction_pwRes, y = pwRes, xlab = 'PRED', ylab = 'PWRES')
G_RES_C_PWRES[[7]] <- RES_PRE(res_7, x = prediction_pwRes, y = pwRes, xlab = 'PRED', ylab = 'PWRES')


G_RES_C_IWRES <- list()
G_RES_C_IWRES[[1]] <- RES_PRE(res_1, x = prediction_iwRes_mean, y = iwRes_mean, xlab = 'PRED', ylab = 'IWRES')
G_RES_C_IWRES[[2]] <- RES_PRE(res_2, x = prediction_iwRes_mean, y = iwRes_mean, xlab = 'PRED', ylab = 'IWRES')
G_RES_C_IWRES[[3]] <- RES_PRE(res_3, x = prediction_iwRes_mean, y = iwRes_mean, xlab = 'PRED', ylab = 'IWRES')
G_RES_C_IWRES[[4]] <- RES_PRE(res_4, x = prediction_iwRes_mean, y = iwRes_mean, xlab = 'PRED', ylab = 'IWRES')
G_RES_C_IWRES[[5]] <- RES_PRE(res_5, x = prediction_iwRes_mean, y = iwRes_mean, xlab = 'PRED', ylab = 'IWRES')
G_RES_C_IWRES[[6]] <- RES_PRE(res_6, x = prediction_iwRes_mean, y = iwRes_mean, xlab = 'PRED', ylab = 'IWRES')
G_RES_C_IWRES[[7]] <- RES_PRE(res_7, x = prediction_iwRes_mean, y = iwRes_mean, xlab = 'PRED', ylab = 'IWRES')



G_RES_C_NPDE <- list()
G_RES_C_NPDE[[1]] <- RES_PRE(res_1, x = prediction_npde, y = npde, xlab = 'PRED', ylab = 'NPDE')
G_RES_C_NPDE[[2]] <- RES_PRE(res_2, x = prediction_npde, y = npde, xlab = 'PRED', ylab = 'NPDE')
G_RES_C_NPDE[[3]] <- RES_PRE(res_3, x = prediction_npde, y = npde, xlab = 'PRED', ylab = 'NPDE')
G_RES_C_NPDE[[4]] <- RES_PRE(res_4, x = prediction_npde, y = npde, xlab = 'PRED', ylab = 'NPDE')
G_RES_C_NPDE[[5]] <- RES_PRE(res_5, x = prediction_npde, y = npde, xlab = 'PRED', ylab = 'NPDE')
G_RES_C_NPDE[[6]] <- RES_PRE(res_6, x = prediction_npde, y = npde, xlab = 'PRED', ylab = 'NPDE')
G_RES_C_NPDE[[7]] <- RES_PRE(res_7, x = prediction_npde, y = npde, xlab = 'PRED', ylab = 'NPDE')



pdf(file = './SCRIPT/PWRES_C.pdf', width = 6, height = 7);{
  gridExtra::grid.arrange(G_RES_C_PWRES[[1]], G_RES_C_PWRES[[2]], G_RES_C_PWRES[[3]], 
                          G_RES_C_PWRES[[4]], G_RES_C_PWRES[[5]], G_RES_C_PWRES[[6]],
                          G_RES_C_PWRES[[7]])
};dev.off()

pdf(file = './SCRIPT/IWRES_C.pdf', width = 6, height = 7);{
  gridExtra::grid.arrange(G_RES_C_IWRES[[1]], G_RES_C_IWRES[[2]], G_RES_C_IWRES[[3]], 
                          G_RES_C_IWRES[[4]], G_RES_C_IWRES[[5]], G_RES_C_IWRES[[6]],
                          G_RES_C_IWRES[[7]])
};dev.off()

pdf(file = './SCRIPT/NPDE_C.pdf', width = 6, height = 7);{
  gridExtra::grid.arrange(G_RES_C_NPDE[[1]], G_RES_C_NPDE[[2]], G_RES_C_NPDE[[3]], 
                          G_RES_C_NPDE[[4]], G_RES_C_NPDE[[5]], G_RES_C_NPDE[[6]],
                          G_RES_C_NPDE[[7]])
};dev.off()



##########################################################################-
# Pruebas de normalidad en los residuales ---------------------------------
##########################################################################-
# A continuación, se realizan pruebas de normalidad univariadas para los 
# residuales en cada uno de los modelos.

##########################################################################-
# Función que aplica una bateróa de pruebas de normalidad a los datos
## *data* ingrese el objeto de datos a analizar.
## *vector* ingrese la posición de las columnas donde se encuentran los 
## residuales de interós.
## *alpha* ingrese el valor de probabilidad para los tests. 
normtest_batery = function(data, vector, alpha){
  
  df = matrix(nrow = length(unique(vector)), ncol = 7)
  for (j in vector) {
    X = dplyr::pull(data,j) # Selecciona como un vector atómico a una columna 
    i = match(j,vector) # Encuentra la posición en el vectror
    df[i,1] = colnames(data[,j])
    df[i,2] = ifelse(shapiro.test(X)$p.value < alpha, '+', '-') # Shapiro-Wilks
    df[i,3] = ifelse(nortest::ad.test(X)$p.value < alpha, '+', '-') #Anderson-Darling
    df[i,4] = ifelse(nortest::cvm.test(X)$p.value < alpha, '+', '-') #Cramer von Mises
    df[i,5] = ifelse(nortest::lillie.test(X)$p.value < alpha, '+', '-') #Liliefors
    df[i,6] = ifelse(nortest::pearson.test(X)$p.value < alpha, '+', '-') #Pearson
    df[i,7] = ifelse(nortest::sf.test(X)$p.value < alpha, '+', '-') #Shapiro Francia
  }
  colnames(df) = c('Variable','Shapiro', 'Anderson_Darling', 'Cramer_von_Mises',
                   'Liliefors','Pearson','Shapiro_Francia')
  return(df)
}

##########################################################################-
# Aplicación de bateróa de pruebas de normalidad a la tabla
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Crear un objeto de tipo lista para almacenar los resultados de pruebas 
##  de normalidad..
##  2 Para todos los elementos de la lista se asigna uno de los data frames 
##  y se aplica la bateróa de pruebas de normalidad a los residuales que se 
##  encuentran en las posiciones 4 (pwRes), 7 (iwRes), y 13 (Npde) del 
##  dataframe.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

norm_resid_list = list()

for (i in 1:7) {
  norm_resid_list[[i]] = normtest_batery(data = get(paste0("res_", i)),
                                         vector = c(4, 7, 13),
                                         alpha = 0.05) %>% 
    as.data.frame()
}

norm_resid_list %>% 
  map_dfr( ~ as.data.frame(.x), .id = 'ID') %>% 
  write_csv(., path = 'SCRIPT/resid_list.txt')

norm_resid_list[[1]]
norm_resid_list[[2]]
norm_resid_list[[3]]
norm_resid_list[[4]]
norm_resid_list[[5]]
norm_resid_list[[6]]
norm_resid_list[[7]]

##########################################################################-
# Función para obtener grófico de tipo QQ con lónea y bandas, dado el modelo 
# el parómetro deseado y color opcional.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Convertir param en expresión
##  2 Crear grófico
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

normal_dist_plot <- function(model, param, colour) {
  param = rlang::ensym(param)
  if(is_missing(colour)){colour = 'blue4'}
  
  g1 <- paste0('res_', i) %>%
    get(.) %>%
    ggplot(mapping = aes(sample = !!param)) +
    stat_qq() +
    qqplotr::stat_qq_line(distribution = 'norm',
                          col = colour) +
    qqplotr::stat_qq_band(distribution = 'norm',
                          fill = colour,
                          alpha = 0.1) +
    geom_abline(slope = 1, intercept = 0, lty = 'dashed') +
    theme_bw() + 
    coord_cartesian(xlim = c(-3, 3), ylim = c(-3, 3)) + 
    annotate(geom = 'label', x = -1.5, y = 2, 
             label = paste0('Modelo: (',letters[model], ')')) +
    xlab('Cuantiles teóricos') + ylab('Cuantiles de muestreo')
  
  return(g1)
}

require(patchwork)

plot_resid_comp <- 
(normal_dist_plot(1, 'iwRes_mean', colour = 'red3') + 
normal_dist_plot(2, 'iwRes_mean', colour = 'red3') + 
normal_dist_plot(3, 'iwRes_mean', colour = 'red3') ) / 
(normal_dist_plot(4, 'iwRes_mean', colour = 'red3') + 
normal_dist_plot(5, 'iwRes_mean', colour = 'red3') + 
normal_dist_plot(6, 'iwRes_mean', colour = 'red3') ) / 
(normal_dist_plot(7, 'iwRes_mean', colour = 'red3')  +
plot_spacer() + plot_spacer()) +
plot_annotation(tag_levels = 'A')

ggsave('IWRES_comp.pdf', plot_resid_comp, 'pdf', 
       './SCRIPT/qqplot', 1, 8, 8, 'in')


for (i in 1:7) {
  normal_dist_plot(i, 'iwRes_mean', colour = 'red3')
  ggsave(filename = paste0('./SCRIPT/qqplot/iwres_', i, '.pdf'), 
        device = 'pdf', width = 2.5, height = 2.5, units = 'in')
}

for (i in 1:7) {
  normal_dist_plot(i, 'pwRes', colour = 'blue3')
  ggsave(filename = paste0('./SCRIPT/qqplot/pwRes_', i, '.pdf'), 
         device = 'pdf', width = 4, height = 4, units = 'in')
}

for (i in 1:7) {
  normal_dist_plot(i, 'npde', colour = 'green3')
  ggsave(filename = paste0('./SCRIPT/qqplot/npde_', i, '.pdf'), 
         device = 'pdf', width = 4, height = 4, units = 'in')
}


gridExtra::grid.arrange(normal_dist_plot(1, 'iwRes_mean'),
                        normal_dist_plot(2, 'iwRes_mean'),
                        normal_dist_plot(3, 'iwRes_mean'),
                        normal_dist_plot(4, 'iwRes_mean'),
                        normal_dist_plot(5, 'iwRes_mean'),
                        normal_dist_plot(6, 'iwRes_mean'),
                        normal_dist_plot(7, 'iwRes_mean'))



# res_1 %>% 
#   filter(iwRes_mean == max(iwRes_mean)) %>% View()

nortest::pearson.test(res_7$iwRes_mean)$p.value
# nortest::lillie.test(res_3$iwRes_mean)$p.value

##########################################################################-
# Apertura de archivos con forma de la distribución -----------------------
##########################################################################-
resid_guides = list()

for (i in 1:7) {
  resid_guides[[i]] <-
    read_csv(file.path(paste0('Modelo_', i), 'RES_M1', 
                       'ChartsData', 'DistributionOfTheResiduals', 
                       'theoreticalGuides.txt'))
}

resid_guid1 <- 
  do.call(rbind.data.frame, resid_guides) %>% 
  add_column(Modelo = rep(1:7, each = 128))


ggplot(resid_guid1, aes(x = abscissa, y = cdf, group = Modelo,
                                 col = Modelo)) +
  geom_line()

##########################################################################-
# 

resid_pdf = list()

for (i in 1:6) {
  resid_pdf[[i]] <-
    read_csv(file.path(paste0('Modelo_', i), 'RES_M1', 
                       'ChartsData', 'DistributionOfTheResiduals', 
                       'y_1_pdf.txt'))
}

resid_pdf[[7]] <- read_csv(file.path('Modelo_7', 'RES_M7', 'ChartsData', 
                      'DistributionOfTheResiduals', 'y_1_pdf.txt'))

# resid_pdf1 <-
#   do.call(rbind.data.frame, resid_pdf) %>%
#   add_column(Modelo = rep(1:7, each = 128))


  





resid_pdf %>% 
  map(~ magrittr::extract(.x)) %>% 
  map(~ dim(.x))




##########################################################################-
# 
resid_cdf = list()

for (i in 1:6) {
  resid_cdf[[i]] <-
    read_csv(file.path(paste0('Modelo_', i), 'RES_M1', 
                       'ChartsData', 'DistributionOfTheResiduals', 
                       'y_1_cdf.txt'))
}

resid_cdf[[7]] <- read_csv(file.path('Modelo_7', 'RES_M7', 'ChartsData', 
                                     'DistributionOfTheResiduals', 'y_1_cdf.txt'))

resid_cdf %>%
  map_dfr( ~ as.data.frame(.x), .id = 'id') %>%
  ggplot(mapping = aes(x = pwRes_abscissa, y = pwRes_cdf, group = id)) +
  geom_line(aes(col = as.factor(id))) +
  scale_color_viridis_d()


##########################################################################-
# Adición de datos de residuales
resid_ls = list()

for (i in 1:7) {
  resid_ls[[i]] <-
    read_csv(file.path(paste0('Modelo_', i), 'RES_M1', 
                       'ChartsData', 'ScatterPlotOfTheResiduals',
                       'y_1_residuals.txt'))
}

xtolet <- function(x) letters[x]
xtolet <- Vectorize(xtolet, vectorize.args = c('x'))

resid_df <- resid_ls %>% 
  map_dfr(~ as.data.frame(.x), .id = 'Model') %>%
  mutate(Model = xtolet(as.numeric(Model))) %>%
  mutate(Model = fct_relabel(Model, ~ paste0('(',  .x, ')'))) %>% 
  as_tibble()

resid_df_SW <- resid_df %>% 
  select(Model, pwRes, iwRes_mean, npde) %>% 
  pivot_longer(cols = pwRes:npde, names_to = 'Residual', 
               values_to = 'Val') %>% 
  group_by(Model, Residual) %>% 
  nest() %>% 
  mutate(SW = map(data, ~shapiro.test(.x$Val)),
         pval = map_dbl(SW, 'p.value'))

##########################################################################-
# Función para generar grófico de distribución de residuales --------------
##########################################################################-
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Conversión de param en una expresión no condicionada
##  2 Creación de tabla con resid_df, que permita obtener media y desvest 
##  para cada modelo en el data.frame.
##  3 Creación de lista Y con data.frames que tienen una abcisa con valores 
##  de -5 a 5, identificador de modelo, medias y desvest diferentes para 
##  cada elemento de la lista.
##  4 Creación de data frame Z que tiene todos los elementos de la lista Y, 
##  con valores de densidad de acuerdo a media y desvest.M.
##  5 Creación de grófico con facet, y colores para cada subgrupo, integra 
##  la información creada por Z.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

resid_dist_plot <- function(data, param) {
  param <- rlang::ensym(param)
  
  X = data %>%
    group_by(Model) %>%
    summarise(mn = mean(!!param),
              sd = sd(!!param))
  
  Y = list()
  for (i in 1:(dim(X)[1])) {
    Y[[i]] = data.frame(
      val = seq(-5, 5, 0.01),
      Model = rep(X[[i, 'Model']], 1001),
      mn = rep(X[[i, 'mn']], 1001),
      sd = rep(X[[i, 'sd']], 1001)
    )
  }
  
  Z = Y %>%
    map_dfr( ~ as.data.frame(.x), .id = 'ID') %>%
    mutate(density = dnorm(x = val, mean = mn, sd = sd)) %>%
    as_tibble()
  
  g1 <-
    ggplot(data, mapping = aes(x = !!param, group = Model)) +
    geom_vline(xintercept = 0) +
    geom_histogram(aes(y = ..density.., fill = Model),
                   bins = 20,
                   col = 'black') +
    geom_line(Z, mapping = aes(x = val, y = density, group = Model)) +
    facet_wrap(. ~ Model, ncol = 3) +
    coord_cartesian(xlim = c(-4, 4), ylim = c(0, 0.8)) +
    scale_fill_viridis_d(option = 'D') +
    theme_bw() +
    theme(legend.position = 'none')
  
  return(g1)
}






##########################################################################-
# Almacenamiento de los gróficos en formato pdf
filter(resid_df, Model != '(g)') %>%
  resid_dist_plot(param = pwRes) + 
  geom_label(
    data = filter(resid_df_SW, Residual == 'pwRes' && Model != '(g)'),
    mapping = aes(
      label = paste('Shapiro-Wilk:', '\n', 
                    'p = ', formatC(pval, 2, format = 'e'))),
    x = 2.0, y = 0.7, size = 4) +
  xlab('WRES') + ylab('Densidad')

ggsave('SCRIPT/Res_PWRES_Dist.pdf', device = 'pdf', width = 5*1.5, height = 4*1.5)


# IWRES_mean
filter(resid_df, Model != '(g)') %>%
  resid_dist_plot(param = iwRes_mean) + 
  geom_label(
    filter(resid_df_SW, Residual == 'iwRes_mean' && Model != '(g)'),
    mapping = aes(
      label = paste('Shapiro-Wilk:', '\n', 
                    'p = ', formatC(pval, 2, format = 'e'))),
    x = 2.0, y = 0.7, size = 4) +
  xlab('IWRES') + ylab('Densidad')

ggsave('SCRIPT/Res_IWRES_Dist.pdf', device = 'pdf', width = 5*1.5, height = 4*1.5)

# NPDE
filter(resid_df, Model != '(g)') %>%
  resid_dist_plot(param = npde) + 
  geom_label(
    filter(resid_df_SW, Residual == 'npde' && Model != '(g)'),
    mapping = aes(
      label = paste('Shapiro-Wilk:', '\n', 
                    'p = ', formatC(pval, 2, format = 'e'))),
    x = 2.0, y = 0.7, size = 4) +
  xlab('NPDE') + ylab('Densidad')

ggsave('SCRIPT/Res_NPDE_Dist.pdf', device = 'pdf', width = 5*1.5, height = 4*1.5)





##########################################################################-
# Grófico de residuales ---------------------------------------------------
##########################################################################-

acp <- function(n, resid) {
  if (n >= 6) {
    stop ('Demasiado Lag')
  }
  
  resid = rlang::ensym(resid) 
  resid1 = paste0(resid, 1)
  
  x = res_1 %>%
    mutate(ID = factor(ID)) %>%
    group_by(ID) %>%
    select(ID, !!resid) %>%
    mutate(resid1 = lead(!!resid, n = n))
  
  y = cor(x[, 2], x[, 3], use = 'complete.obs')
  
  return(y)
}

acp = Vectorize(acp, vectorize.args = c("n"))

ACF <- data.frame('Lag' = 1:5) %>%
  mutate(pwRes = acp(n = Lag, resid = 'pwRes')) %>%
  mutate(iwRes = acp(n = Lag, resid = 'iwRes_mean')) %>%
  mutate(NPDE = acp(n = Lag, resid = 'npde'))


ACF %>%
  gather(`pwRes`, `iwRes`, `NPDE`, key = Residual, value = ACF) %>% 
  ggplot(mapping = aes(x = Lag,
                       y = ACF,
                       group = Residual,
                       fill = Residual)) + 
  geom_bar(stat = 'identity') + 
  facet_grid(.~Residual)

ggsave(filename = './SCRIPT/ACF_PLOT.pdf', device = 'pdf', 
       width = 4, height = 4)


































