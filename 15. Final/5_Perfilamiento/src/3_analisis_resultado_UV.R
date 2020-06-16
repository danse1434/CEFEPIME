##========================================================================#
## Nombre del Script: Análisis de mapeo de la función objetivo
##  
## Propósito del Script: El mapeo de OFV permite conocer si se ha alcanzado 
## un mínimo global del modelo, además permite conocer IC no asintóticos 
## alrededor de los parámetros del modelo. En este script se analizan los 
## resultados obtenidos mediante la evaluación de verosimilitud con la suite 
## de Monolix. En este script se leen datos y se produce un gráfico con 
## funciones de optimización.

## Autor: Daniel S. Parra González 
## Fecha de creación: 07-03-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##========================================================================#

#-------------------------------------------------------------------------------#
# Introducción ------------------------------------------------------------
#-------------------------------------------------------------------------------#
# Carga de paquetes
require(patchwork)
require(rlang)
require(tidyverse)

# Definición de directorio principal
aux_dir =  file.path(getwd(), '1_univariado')

# Lectura de archivo funciones
source('src/6_funciones.R', encoding = 'UTF-8')

# Apertura de archivo de datos de parámetros de modelo base final
populationParameters <-
  read_csv("data/populationParameters.txt")

popParam <- populationParameters %>% 
  select(parameter, value) %>% 
  pivot_wider(names_from = parameter, values_from = value)

#-------------------------------------------------------------------------------#
# Lectura de archivos de verosimilitud ------------------------------------
#-------------------------------------------------------------------------------#

df_Clpop      <- extractor('Cl_pop')
df_V1pop      <- extractor('V1_pop')
df_V2pop      <- extractor('V2_pop')
df_beta_Clpop <- extractor('beta_Cl_tSCRMGDL')
df_a          <- extractor('a')

# Cálculo de valores mínimos y adición de 3.84

df_Clpop1 <- df_Clpop %>%
  xval.func(Cl_pop, LL1)
# Mínimo: (19.3, 0.00); Corte: (18.9, 3.84) y (19.8, 3.84); t0: (20.6)

df_V1pop1 <- df_V1pop %>%
  xval.func(V1_pop, LL1)
# Mínimo: (18.9, 0.00); Corte: (18.0, 3.84) y (19.8, 3.84); t0: (23.8)

df_V2pop1 <- df_V2pop %>%
  xval.func(V2_pop, LL1)
# Mínimo: (6.92, 0.00); Corte: (6.65, 3.84) y (8.15, 3.84); t0: (13.2)

df_beta_Clpop1 <- df_beta_Clpop %>%
  xval.func(beta_Cl_tSCRMGDL, LL1)
# Mínimo: (-0.468, 0.00); Corte: (-0.489, 3.84) y (-0.446, 3.84); t0: (-0.415)

# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

# Gráficos
aux_plot <- list(geom_line(),
                 ylab('LL - min(LL)'),
                 geom_hline(yintercept = 3.84, lty = 'dotted'))

G_ll_1 <- ggplot(df_Clpop, aes(x = Cl_pop, y = LL1)) +
  geom_point(data = df_Clpop1$Minimo, col = 'red3', shape = 8) +
  geom_point(data = df_Clpop1$Punto_Corte, col = 'green1') +
  geom_vline(data = popParam, aes(xintercept = Cl_pop), lty='dashed', col='blue4') + 
  xlab(expression('Cl:'~theta[0]~'(L/h)')) + aux_plot

G_ll_2 <- ggplot(df_beta_Clpop, aes(x = beta_Cl_tSCRMGDL, y = LL1)) +
  geom_point(data = df_beta_Clpop1$Minimo, col = 'red3', shape = 8) +
  geom_point(data = df_beta_Clpop1$Punto_Corte, col = 'green1') +
  geom_vline(data = popParam, aes(xintercept = beta_Cl_tSCRMGDL), lty='dashed', col='blue4') + 
  xlab(expression('Cl:'~theta[1])) + aux_plot

G_ll_3 <- ggplot(df_V1pop, aes(x = V1_pop, y = LL1)) +
  geom_point(data = df_V1pop1$Minimo, col = 'red3', shape = 8) +
  geom_point(data = df_V1pop1$Punto_Corte, col = 'green1') +
  geom_vline(data = popParam, aes(xintercept = V1_pop), lty='dashed', col='blue4') + 
  xlab(expression(V[1]~'(L)')) + aux_plot

G_ll_4 <- ggplot(df_V2pop, aes(x = V2_pop, y = LL1)) +
  geom_point(data = df_V2pop1$Minimo, col = 'red3', shape = 8) +
  geom_point(data = df_V2pop1$Punto_Corte, col = 'green1') +
  geom_vline(data = popParam, aes(xintercept = V2_pop), lty='dashed', col='blue4') + 
  xlab(expression(V[2]~'(L)')) + aux_plot

G_ll <- ((G_ll_1 + G_ll_2) / (G_ll_3 + G_ll_4)) + 
  plot_annotation(tag_levels = 'A')

# Almacenamiento de gráfico
ggsave(G_ll, filename = 'figures/1_perfiles_LL_UV.pdf', device = 'pdf', 
       width = 7, height = 5)




df_Qpop      <- extractor('Q_pop')

ggplot(df_Qpop, aes(x = Q_pop, y = LL1)) + 
  geom_line() +
  geom_vline(data = popParam, aes(xintercept = Q_pop), lty='dashed', col='blue4')
