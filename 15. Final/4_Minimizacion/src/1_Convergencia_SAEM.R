##------------------------------------------------------------------------------#
## Nombre del Script: Análisis de datos de convergencia SAEM --------------------
##  
## Propósito del Script: se realiza un análisis de convergencia del algoritmo 
## SAEM con diversos valores iniciales de los parámetros del modelo. Esto para el 
## modelo final con covariables Cl-SCR.
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 13-06-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
require(rlang)
require(scales)
require(tidyverse)

# Vector con órden de parámetros
level_par <- c('Cl_pop', 'beta_Cl_tSCRMGDL', 'V1_pop', 'Q_pop', 'V2_pop', 
               'omega_Cl', 'omega_V1', 'omega_V2', 'omega_Q', 'a')

# Selección directorio auxiliar
aux_dir <- file.path(getwd(), 'FACTORIAL')

#-------------------------------------------------------------------------------#
# 1 Resumen Parámetros de Diseño Factorial --------------------------------------
#-------------------------------------------------------------------------------#
# Lista vacía *d_fct_ls*
d_fct_ls <- vector('list', 27)
# Lectura de datos de parámetros estimados
for (i in 1:27) {
  d_fct_ls[[i]] <-
    read_csv(file.path(aux_dir, paste0('A_', i), 'populationParameters.txt'),
             col_types = cols())
}
# Conversión en data.frame unificado
d_fct_ls1 <- d_fct_ls %>%
  map_dfr(~ as_tibble(.x), .id = 'A') %>%
  mutate(
    A = as.double(A),
    parameter = factor(parameter, levels = level_par))

# Ajuste configuración de gráfico
theme_set(theme_bw())

# Gráfico de resumen de parámetros por corrida
G1 <- d_fct_ls1 %>% 
  ggplot(aes(x = A, y = value, group = A)) +
  geom_point(colour = "#1c86ee") + 
  geom_errorbar(aes(ymin = value-se_sa, ymax = value+se_sa),colour = "#1c86ee") +
  facet_wrap(~parameter, scales = 'free_y') +
  # scale_color_gradientn(colours = rainbow(7) ) +
  theme(legend.position = 'none', 
        axis.title = element_blank())

G1

# Almacenamiento de gráficos
ggsave('1_parm_indic.pdf', G1, 'pdf', 'figures', 1, 7, 5)
ggsave('1_parm_indic.png', G1, 'png', 'figures', 1.5, dpi = 'print')

#-------------------------------------------------------------------------------#
# 2 Convergencia Diseño Factorial -----------------------------------------------
#-------------------------------------------------------------------------------#
# Creación de vector vacío
s_fct_ls <- vector('list', 27)
# Lectura de datos de convergencia SAEM
for (i in 1:27) {
  s_fct_ls[[i]] <-
    read_csv(file.path(aux_dir, paste0('A_', i), 'ChartsData', 'Saem', 
                       'CvParam.txt'), col_types = cols())
}
# Conversión en data.frame unificado
s_fct_ls1 <- s_fct_ls %>%
  map_dfr(~ as_tibble(.x), .id = 'I') %>%
  rename(LL = convergenceIndicator) %>% 
  pivot_longer(
    cols = matches('\\_pop|omega\\_|^a$|beta|^LL$'),
    names_to = 'parameter',
    values_to = 'value'
  ) %>% 
  mutate(I = as.integer(I), 
         parameter = factor(parameter, levels = c(level_par, 'LL')))

# Gráfico de convergencia SAEM por corrida


G2 <- s_fct_ls1 %>%
  filter((iteration + 9) %% 10 == 0) %>%
  ggplot(aes(x = iteration, y = value, group = I)) +
  geom_line(colour = 'gray50', alpha = 0.4) +
  facet_wrap( ~ parameter, scales = 'free_y') +
  # scale_color_gradientn(colours = rainbow(7) ) +
  scale_x_continuous(
    labels = function(x)
      format(x, scientific = TRUE)
  ) + theme(
    legend.position = 'none',
    axis.title = element_blank(),
    axis.text.x = element_text(angle = -45),
    panel.grid = element_blank()
  )

G2

# Almacenamiento de gráficos
ggsave('2_conv_indic.pdf', G2, 'pdf', 'figures', 1, 7, 5)
ggsave('2_conv_indic.png', G2, 'png', 'figures', 
       dpi = 'print', width = 7, height = 5)

#-------------------------------------------------------------------------------#
# 3 Resumen parámetros herramienta Monolix --------------------------------------
#-------------------------------------------------------------------------------#
aux_dir <- file.path(getwd(), 'BASE_MODEL', 'Assessment')
# Lista vacía *d_mon_ls*
d_mon_ls <- vector('list', 20)
# Lectura de datos de parámetros estimados
for (i in 1:20) {
  d_mon_ls[[i]] <-
    read_csv(file.path(aux_dir, paste0('Run', i), 'populationParameters.txt'),
             col_types = cols())
}
# Creación tabla unificada
d_mon_ls1 <- d_mon_ls %>%
  map_dfr(~ as_tibble(.x), .id = 'I') %>% 
  mutate(I = as.integer(I),
         parameter = factor(parameter, levels = level_par))

# Gráfico de resumen por cada corrida
G3 <- d_mon_ls1 %>% 
  ggplot(aes(x = I, y = value, col = I)) +
  geom_point(colour = "#1c86ee") + 
  geom_errorbar(aes(ymin = value-se_lin, ymax = value+se_lin), colour = "#1c86ee") +
  facet_wrap(~parameter, scales = 'free_y') +
  scale_color_gradientn(colours = rainbow(7) ) +
  theme(legend.position = 'none', 
        axis.title = element_blank())

# Almacenamiento de gráficos
ggsave('3_parm_indic.pdf', G3, 'pdf', 'figures', 1, 7, 5)
ggsave('3_parm_indic.png', G3, 'png', 'figures', 1.5, dpi = 'print')
