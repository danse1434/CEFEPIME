##=========================================================================#
## Nombre del Script: Obtencion de gráficos a partir de datos de figuras---
## generados por Monolix GUI. 
##  
## Proposito del Script: crear y almacenar gráficos generados a partir de 
##  estudios de la distribución de efectos aleatorios de parámetros. Estos 
##  se estudian mediante simulación con MCMC desde la distribución condicional. 
##  Se muestran correlogramas y gráficos de distribución mediante diagrama de 
##  cajas.
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion: 03-mar-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##=========================================================================#
##########################################################################-
# Introducción -----------------------------------------------------
##########################################################################-
# Carga de paquetes
require(rlang)
require(tidyverse)
require(grid)
##########################################################################-
# Selección de directorio principal
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '06. Correlacion_ETA', 'BASE_CONTROL'))


##########################################################################-
# Gráficos de correlación de Parámetros -----------------------------------
##########################################################################-
corr_dir <- file.path('ChartsData','CorrelationBetweenRandomEffects')

# Tabla con valores eta simulados
data_corr_1 <- 
  read_csv(file.path(corr_dir, 'simulatedEta.txt'))

# Tabla con valores de eta porcada individuo
data_corr_2 <- 
  read_csv(file.path(corr_dir, 'eta.txt'))

# Tabla con ayudas visuales a gráficos de correlación
data_corr_3 <- 
  read_csv(file.path(corr_dir, 'visualGuides.txt'))

# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

# Función de creación de diagramas de dispersión
# Se especifican las variables a modelar
corr_plot <- function(x, y) {
  x = rlang::ensym(x)
  y = rlang::ensym(y)
  G = data_corr_1 %>%
    mutate(ID = factor(ID)) %>%
    ggplot(aes(x = !!x, y = !!y)) +
    geom_point(aes(col = ID)) +
    stat_smooth(method = 'lm') +
    scale_color_viridis_d() +
    xlab('') + ylab('') +
    theme(legend.position = "none")
  return(G)
}
  
##########################################################################-
# Función de creación de histogramas
# Se especifican las variables a modelar
# Se crea una serie distribuciones para cada simulación de parámetro eta 
# por MCMC.
# Se calculan las medias y desviaciones estándar por grupo y se adicionan 
# al gráfico en forma de un loop.

hist_plot <- function(x, parameter) {
  x = rlang::ensym(x)
  G = data_corr_1 %>% 
    mutate(ID = factor(ID)) %>% 
    ggplot(aes(x = !!x)) + 
    geom_histogram(aes(y = ..density.., fill = ID)) +  
    xlab('') + ylab('') +
    theme(legend.position = "none") + 
    scale_fill_viridis_d()
  
  Y = list()
  
  for (i in 1:15) {
    X = data_corr_1 %>%
      filter(ID == i) %>%
      summarise(mean = mean(!!x),
                sd = sd(!!x))
    
    
    Y[[i]] = stat_function(
      fun = dnorm,
      args = list(mean = X[['mean']],
                  sd = X[['sd']]),
      size = 0.7
    )
  }
  
  H = G + Y + annotate(geom = "text",
    label = parameter,
    x = 0.8 * mean(pull(data_corr_1, !!x)),
    y = Inf, hjust = 0, vjust = 1)
  
  return(H)
}

# Ejemplo de generación de gráfico con título en letras griegas
# hist_plot(eta_V1_simulated, 
#           expression( eta~(V[1]))) 
# Ejemplo de funcionamiento de gráfico de dispersión.
# corr_plot(eta_Cl_simulated, eta_Q_simulated)
# dev.off()


##########################################################################-
# Creación de gráfico
# Se utilizan las funciones creadas con anterioridad para generar histogramas 
# y diagramas de dispersión con regresión lineal incorporada.

pdf(file = 'FIGURAS/Correlation_Plot.pdf', width = 6, height = 5);{
grid::pushViewport(grid::viewport(layout = grid.layout(4, 4)))
vplayout <- function(x, y) grid::viewport(layout.pos.row = x, layout.pos.col = y)

# Función para especificar gráfico y posición en la malla cartesiana. 
pr_f <- function(plot, x, y) {
  print(plot, vp = vplayout(x, y))
}

# Especificación de gráficos y posiciones. 
pr_f(hist_plot(eta_Cl_simulated,
                  expression(eta ~ (Cl))) , 1, 1)
pr_f(hist_plot(eta_V1_simulated,
                  expression(eta ~ (V[1]))) , 2, 2)
pr_f(hist_plot(eta_Q_simulated,
                  expression(eta ~ (Q))) , 3, 3)
pr_f(hist_plot(eta_V2_simulated,
                  expression(eta ~ (V[2]))) , 4, 4)
pr_f(corr_plot(eta_Cl_simulated, eta_V1_simulated), 2, 1)
pr_f(corr_plot(eta_Cl_simulated, eta_Q_simulated), 3, 1)
pr_f(corr_plot(eta_Cl_simulated, eta_V2_simulated), 4, 1)
pr_f(corr_plot(eta_V1_simulated, eta_Q_simulated), 3, 2)
pr_f(corr_plot(eta_V1_simulated, eta_V2_simulated), 4, 2)
pr_f(corr_plot(eta_Q_simulated, eta_V2_simulated), 4, 3)
pr_f(corr_plot(eta_V1_simulated, eta_Cl_simulated), 1, 2)
pr_f(corr_plot(eta_Q_simulated, eta_Cl_simulated), 1, 3)
pr_f(corr_plot(eta_V2_simulated, eta_Cl_simulated), 1, 4)
pr_f(corr_plot(eta_Q_simulated, eta_V1_simulated), 2, 3)
pr_f(corr_plot(eta_V2_simulated, eta_V1_simulated), 2, 4)
pr_f(corr_plot(eta_V2_simulated, eta_Q_simulated), 3, 4)
}; dev.off()

# ggsave('FIGURAS/Correlation_Plot.pdf', device = 'pdf', width = 6,height =4)


##########################################################################-
# Distribución de Eta -----------------------------------------------------
##########################################################################-
# Especificación del directorio de extracción de datos
esteta_dir <-
  file.path('ChartsData', 'DistributionOfTheStandardizedRandomEffects')

# Lectura de archivos de datos contenidos en el directorio
data_esteta_1 <-
  read_csv(file.path(esteta_dir, 'simulatedStandardizedEta.txt'))

data_esteta_2 <-
  read_csv(file.path(esteta_dir, 'StandardizedEta.txt'))

##########################################################################-
# Especificación de parámetros de gráfico
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Eliminar columna filtro
##  2 Convertir en formato colapsado
##  3 Eliminar las instancias en donde aparece "stand" o "_simulated"
##  4 Convertir la variable parameter en un factor ordenado
##  5 Generar graficos
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

G_boxplot <-
data_esteta_1 %>%
  select(-filter) %>%
  gather(matches("_simulated"), key = 'parameter', value = 'eta_value') %>%
  mutate(parameter = str_replace(parameter, "stand", "")) %>%
  mutate(parameter = str_replace(parameter, "_simulated", "")) %>%
  mutate(parameter = factor(parameter, levels = c('Eta_Cl', 'Eta_V1', 'Eta_Q', 'Eta_V2'))) %>%
  ggplot(aes(x = parameter, y = eta_value, col = parameter)) +
  geom_hline(yintercept = 0) +
  # IC95% para distribución normal estándar
  geom_hline(yintercept = c(-1.960, 1.960), lty = 'dashed', size = 0.5) +
  # IC99% para distribución normal estándar
  geom_hline(yintercept = c(-2.576, 2.576), lty = 'dashed', size = 0.25) +
  # Violin, Puntos, y 
  geom_violin() +
  geom_jitter(shape = 16, width = 0.1, height = 0, size = 1) + 
  geom_boxplot(
    width = .1,
    position = "dodge",
    outlier.shape = 4,
    outlier.colour = "red",
    colour = "black"
  ) +
  theme(legend.position = "none") + 
  coord_cartesian(ylim = c(-5, 5)) +
  ylab(expression(eta["std"])) + xlab('') + 
  scale_color_viridis_d(end = 0.9)

ggsave(G_boxplot,
  filename = 'FIGURAS/boxplot_eta_est.pdf',
  device = 'pdf', width = 5/1.5, height = 4/1.5)










