##------------------------------------------------------------------------#
## Nombre del Script: An?lisis de sesgo con m?todos de censura de datos ---
##  
## Proposito del Script: An?lisis de sesgo con m?todos de censura de datos
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
require(ggrepel)
##########################################################################-
# Selecci?n de directorio de trabajo
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
##  3 Para las variables con el patr?n 'M\\d\\_CENS' se recalculan errores 
##  relativo a la variable Ref.
##  4 Eliminar las variable Ref, y M1_CENS.
##  5 Renombrar `Project Name` como Parametro.
##  6 Ordenar los valores de par?metros de acuerdo a un vector de ordenacion.
##  7 Colapsar la tabla en forma de M?todo con valor de sesgo.
##  8 Eliminar la expresi?n \\_CENS de las variables en la columna Metodo.
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
# Crear un gr?fico con los sesgos en la estimaci?n de modelos
G1 <- Modelo %>% 
  ggplot(mapping = aes(x = Metodo, y = Sesgo)) + 
  theme_classic() +
  geom_bar(stat = 'identity', colour = 'black', 
           aes(fill = ifelse(abs(Sesgo) > 0.15, 'out', 'in'))) +
  facet_wrap(.~Parametro, ncol = 4) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-0.15, +0.15), lty = 'dashed') +
  geom_hline(yintercept = c(-0.20, +0.20), lty = 'dashed') +
  geom_label_repel(
    data = filter(Modelo, Sesgo <= 0),
    mapping = aes(y = Sesgo * 1.5, label = round(Sesgo, 2),
      col = ifelse(abs(Sesgo) > 0.15, 'O', 'I')),
    fill = 'white', direction = 'y', ylim = c(NA, 0)) +
  
  geom_label_repel(
    data = filter(Modelo, Sesgo > 0),
    mapping = aes(y = Sesgo * 1.5, label = round(Sesgo, 2),
                  col = ifelse(abs(Sesgo) > 0.15, 'O', 'I')),
    fill = 'white', direction = 'y', ylim = c(0, NA)) +
  
  scale_fill_manual(values = c('blue', 'red')) +
  scale_color_manual(values = c('black', 'red')) +
  xlab('Método') +
  coord_cartesian(ylim = c(-0.30, 0.30))+
  theme(legend.position = 'none',
        panel.border = element_rect(fill = NA, colour = 'black'))

G1

##########################################################################-
# Almacenar gr?fico en pdf
ggsave('./FIGURAS/efecto_censura_sesgo.pdf', G1, device = 'pdf', 
       width = 6*1.2, height = 5*1.2)











