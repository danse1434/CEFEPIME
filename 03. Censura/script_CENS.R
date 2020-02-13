##------------------------------------------------------------------------#
## Nombre del Script: Análisis de sesgo con métodos de censura de datos ---
##  
## Proposito del Script: Análisis de sesgo con métodos de censura de datos
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
require(grid)
require(gridExtra)
require(readxl)
##########################################################################-
# Selección de directorio de trabajo
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '03. Censura'))
##########################################################################-
# Apertura de archivos
Modelos_Censura <- read_excel("Modelos_Censura.xlsx", 
                              sheet = "2")

##########################################################################-
# Preparar la tabla de sesgos de modelos
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Seleccionar modelo.
##  2 Crear una variable dummy denominada *Ref*.
##  3 Para las variables con el patrón 'M\\d\\_CENS' se recalculan errores 
##  relativo a la variable Ref.
##  4 Eliminar las variable Ref, y M1_CENS.
##  5 Renombrar `Project Name` como Parametro.
##  6 Ordenar los valores de parámetros de acuerdo a un vector de ordenacion.
##  7 Colapsar la tabla en forma de Método con valor de sesgo.
##  8 Eliminar la expresión \\_CENS de las variables en la columna Metodo.
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

Modelo <-
  Modelos_Censura %>% 
  mutate(Ref = M1_CENS) %>% 
  mutate_at(vars(matches("M\\d\\_CENS")), list( ~ ((. - Ref) / Ref))) %>%
  select(-Ref, -M1_CENS) %>%
  rename(Parametro = `Project Name`) %>% 
  mutate(Parametro = factor(Parametro,
    levels = c('Cl_pop','V1_pop','Q_pop','V2_pop',
               'omega_Cl','omega_V1','omega_Q','omega_V2',
               'a') ) ) %>% 
  gather('M2_CENS', 'M3_CENS', 'M4_CENS', key = 'Metodo', value = 'Sesgo') %>% 
  mutate(Metodo = str_replace(Metodo, "\\_CENS", "")) 

##########################################################################-
# Crear un gráfico con los sesgos en la estimación de modelos
G1 <-
  Modelo %>% 
  ggplot(mapping = aes(x = Metodo, y = Sesgo, fill = Metodo)) + 
  theme_bw() +
  coord_cartesian(ylim = c(-0.2, 0.2))+
  geom_bar(stat = 'identity', colour = 'black') +
  facet_wrap(.~Parametro, ncol = 4) +
  geom_hline(yintercept = 0) +
  theme(legend.position = 'none', panel.grid = element_blank())

##########################################################################-
# Almacenar gráfico en pdf
ggsave('./FIGURAS/efecto_censura_sesgo.pdf', G1, device = 'pdf', 
       width = 6, height = 5)











