##------------------------------------------------------------------------#
## Nombre del Script: Análisis de evaluación de convergencia de modelo base 
##  de Monolix  -----------------------------------------------------------
##  
## Propósito del Script: generar gráficas para el análisis de la evaluación 
##  de convergencia de modelo base de Monolix.  
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion: 04/03/2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
# Carga de paquetes 
require(rlang)
require(tidyverse)

##########################################################################-
# Introducción -----------------------------------------------------
##########################################################################-
# Selección de directorio principal
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '07. Minimizacion', '1_MINIMIZACION', 'BASE_CONTROL', 
                'Assessment'))

##########################################################################-
# Apertura de archivos de parámetros poblacional
# Se abren los archivos para las 30 evaluaciones de convergencia realizadas 
# con diversos parámetros del modelo, y semillas.
popparamL <- vector(mode = "list", length = 30)
# Apertura con for loop se accede a las carpetas mediante la iteración
for (i in 1:30) {
  popparamL[[i]] <-
    read_csv(file.path(paste0('Run', i), 'populationParameters.txt'))
}

# Transformación de lista con tablas a una tabla única
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Convertir lista en una tabla única, se adiciona una columna que indica 
##  la iteración correspondiente.
##  2 Transformar en una tibble
##  3 Convertir la iteración en tipo número
##  4 Tranformar la columna parameter en un factor ordenado
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

popparam <- popparamL %>% 
  map_dfr( ~ as.data.frame(.x), .id = 'Iteration') %>% 
  as_tibble(.) %>% 
  mutate(Iteration = as.double(Iteration)) %>% 
  mutate(parameter = 
           factor(
             parameter,
             levels = c('Cl_pop', 'V1_pop', 'Q_pop', 'V2_pop', 
                        'omega_Cl', 'omega_V1', 'omega_Q', 'omega_V2',
                        'a', 'b')))

##########################################################################-
# Gráfico con iteraciones -------------------------------------------------
##########################################################################-
# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

##########################################################################-
# Resumen de medianas para cada parámetro en específico.
# Agrupar por parámetro; resumir por mediana de valor en todas las iteraciones
popparamDF <- popparam %>%
  group_by(parameter) %>%
  summarise(mn = mean(value))

# Gráfico de distribución de convergencia con barras de error SE
G1 <- popparam %>% 
  ggplot(aes(x = Iteration, col = factor(Iteration))) + 
  geom_errorbar(aes(ymin = value - se_sa, ymax = value + se_sa)) +
  geom_point(aes(y = value)) + 
  geom_hline(data = popparamDF, mapping = aes(yintercept = mn)) +
  facet_wrap(~parameter, ncol = 4, scales = "free_y") + 
  theme(legend.position = "none") +
  xlab('Iteración') + ylab('')

# Almacenamiento de Gráfico
ggsave('eval_parametros.pdf', plot = G1, device = 'pdf', 
       width = 8, height = 6)

##########################################################################-
# Criterios de información ------------------------------------------------
##########################################################################-
# Apertura iterativa de archivos con criterios de información
convassesL <- vector(mode = "list", length = 30)

for (i in 1:30) {
  convassesL[[i]] <-
    read_csv(file.path(paste0('Run', i), 'LogLikelihood', 'logLikelihood.txt'))
}

##########################################################################-
# Conversión en archivo de gráfico
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Convertir en una tabla, con columna que indica iteración de origen
##  2 Convertir en tibble
##  3 Cambiar iteración por variable numérica
##  4 Renombrar la columna IS
##  5 Eliminar el criterio 'stdError'
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

convassesDF <- convassesL %>% 
  map_dfr( ~ as.data.frame(.x), .id = 'Iteration') %>% 
  as_tibble(.) %>% 
  mutate(Iteration = as.double(Iteration)) %>% 
  rename(IS = importanceSampling) %>% 
  filter(criteria != 'stdError')

# Gráfico con -2LL alcanzada por cambio de valores iniciales
# Sólo se muestra este criterio porque los demás son idénticos y no aportan 
# información sobre mejorías por el cambio en el conjunto de valores iniciales. 
# En adición, se muestra una línea que indica la disminución de 3.84 desde el 
# valor máximo de -2LL (esta es el nivel de disminución referido como 
# significativo en un test LRT con 1DF a alpha de 0.05).

G2 <- convassesDF %>%
  filter(!(criteria %in% c('AIC', 'BIC', 'BICc'))) %>% 
  ggplot(aes(x = Iteration, y = IS)) +
  geom_point(aes(col = factor(Iteration))) + 
  # geom_line() + 
  theme(legend.position = "none") +
  geom_hline(yintercept = 513.0 - 3.84, col = 'gray80', lty = 'dashed') +
  xlab('Iteración') + ylab('-2LL')

# Almacenamiento de Gráfico
ggsave('eval_OFV.pdf', plot = G2, device = 'pdf', 
       width = 6, height = 4)




