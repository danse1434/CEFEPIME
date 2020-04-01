##------------------------------------------------------------------------#
## Nombre del Script: Revisión de trayectoria de algoritmos de evaluación 
##  covariables  ----------------------------------------------------
##  
## Propósito del Script: revisión de comportamiento de algoritmo COSSAC con 
##  seguimiento de BICc con restricciones en las covariables.
##  
## Autor: Daniel S. Parra González 
## Fecha de creación: 31-03-2020 
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
# Selección de directorio de trabajo
setwd(file.path('F:','Documentos','(Proyecto)_Estudio_PKPD','CEFEPIME', 
                '12. Covariables', 'BASE_MODEL', 'M_COSSAC_LRT'))

# Carga de paquetes
require(tidyverse)
require(ggplot2)
require(gridExtra)

##########################################################################-
# Introducción ------------------------------------------------------------
##########################################################################-
# Lectura de archivo de datos
# El archivo de datos original, se llama "arguments.dat" y tiene una 
# codificación desconocida, por lo cual se pasó a una hoja de cálculo, y se
# realizó un preprocesamiento con separación de las líneas por "|" y 
# conversión a un archivo "*.csv".

data <- read_csv(".Internals/Argumentos_1.csv",
                 skip = 12)

# Se eliminan las primeras 12 filas del archivo

##########################################################################-
# Modificación de archivo de datos ----------------------------------------
##########################################################################-
# Eliminación de espacios blancos en nombres de la tabla
dcol <- colnames(data) %>% 
  str_replace_all(" ", "") # Remoción espacios

# Desplazamiento una celda a la derecha de los nombres
for (i in 2:29) {
  colnames(data)[1] <- 'Intro' 
  colnames(data)[i] <- dcol[i-1]
}

##########################################################################-
# Creación de objeto *data1*
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Detección de filas con letras en la columna X1
##  2 Renombrar a X1 como la variable Parametro
##  3 Remoción de columna Intro
##  4 Eliminación de espacios blancos en la variable parámetro
##  5 Eliminación de encabezados con la palabra clave !AGEA
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data1 <- data %>% 
  filter(str_detect(X1, "\\w")) %>% 
  rename(Parametro = X1) %>% 
  select(-Intro) %>% 
  mutate(Parametro = str_replace_all(Parametro, "\\s", "")) %>% 
  filter(!str_detect(Parametro, "AGEA")) 
  
##########################################################################-
# Modificación a objeto *data2*
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Añadir columna con valores de iteración repetiendo 4 veces (por cada 
##  parámetro), la secuencia que indica la iteración.
##  2 Seleccionar covariables evaluadas en el assesment
##  3 Cambiar por 0 aquellos lugares donde exista NA en la tabla (covariable 
##  ausente para el parámetro), y 1 donde no (covariable presente para el 
##  parámetro).
##  4 Colapsar la tabla por covariables y cumplimiento de criterio, se tienen 
##  ahora dos columnas que indican el para cov-param involucrado por 
##  iteración.
##  5 Concatenar las columnas Parámetro y Covariable en el *par* involucrado
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data2 <- data1 %>% 
  add_column(Iteracion = rep(seq(1, dim(data1)[1]/4, 1), each = 4), 
             .before = "Parametro") %>% 
  select(matches("logt|ANTU|LLP|LMP|SEXF|Iteracion|Parametro")) %>% 
  mutate_at(vars(matches("logt|ANTU|LLP|LMP|SEXF")), 
            ~ ifelse(is.na(.), 0, 1)) %>% 
  gather("Covariable", "Criterio", matches("logt|ANTU|LLP|LMP|SEXF")) %>% 
  unite(Parametro, Covariable, col = "Par", sep = '_')

##########################################################################-
# Gráfico con pares estudiados ----------------------------------------------
##########################################################################-
# Modificación de objeto a *data3*
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Modificación de *Par* a factor.
##  2 Creación de *Par1* que anonimiza los factores en enteros.
##  3 Modificación de *Criterio* a factor.
##  4 Separación de Par a dos columnas teniendo en cuenta el separador "_"
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data3 <- data2 %>%
  mutate(Par = factor(Par),
         Par1 = as.integer(fct_anon(Par)),
         Criterio = factor(Criterio)) %>% 
  separate(Par, into = c('Parametro', 'Covariable'), sep = "\\_", 
           remove = FALSE)

##########################################################################-
# Color asignado de acuerdo a criterio incluido (1) o no (0)
Criterio <- c("1" = 'red', "0" = NA) # Cuando no se incluye el par no se muestra

# Gráfico 1 - Muestra las covariables incluidas en el modelo

G1 <- data3 %>%
  ggplot(aes(x = Iteracion, y = Par1)) + 
  geom_tile(aes(fill = Criterio)) +
  theme_classic() +
  xlab("Iteración") + ylab('Parámetro - Covariable') +
  scale_y_continuous(breaks = seq(0, 60, 10), 
                     minor_breaks = seq(1, 60, 1), 
                     sec.axis = dup_axis(breaks = NULL, name = NULL)) +
  scale_x_continuous(position = 'top',
                     breaks = round(seq(0, max(data3$Iteracion), length.out = 20)),
                     sec.axis = dup_axis(breaks = NULL, name = NULL)) +
  scale_fill_manual(values = Criterio) +
  theme(legend.position = "none", 
        axis.text.x = element_text(angle = 90),
        panel.grid.major.x = element_line(colour = 'gray80', size = 0.01),
        panel.grid.minor.y = element_line(colour = 'gray95', size = 0.01))

# ggsave(filename = "Heatmap.pdf", plot = G1, device = 'pdf', width = 8, height = 8)  
  

##########################################################################-
# Gráfico con trayectorias de convergencia --------------------------------
##########################################################################-
# Modificación de objeto a *data4*
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Selección de la primera columna
##  2 Filtrar filas que contengan las expresiones "L L" y "B I", estas 
##  tienen los criterios de convergencia alcanzados en cada paso.
##  3 Eliminar espacios blancos en columna
##  4 Separar la variable en *Parametro* y *Valor*
##  5 Convertir la variable *valor* en número
##  6 Adicionar una columna que indique la iteración, como se eliminó la 
##  primera iteración con el comando skip en la lectura se arranca desde 2.
##  7 Convertir los valores `-2*LL` en LRT en la variable *Parametro*
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data4 <- data %>% 
  select(Intro) %>% 
  filter(str_detect(Intro, "L\\sL|B\\sI")) %>% 
  mutate(Intro = str_replace_all(Intro, "\\s", "")) %>% 
  separate(Intro, into = c('Parametro', 'Valor'), "\\:") %>% 
  mutate(Valor = as.double(Valor)) %>% 
  add_column(Iteracion = rep(seq(2, (dim(.)[1] / 2)+1, 1), each = 2), 
             .before = "Parametro") %>% 
  mutate(Parametro = if_else(Parametro == "BICc", 'BICc', 'LRT'))

##########################################################################-
# Gráfico 2 - Muestra las trayectorias de los criterios de convergencia
# Selección mejor modelo
h <- 10

G2 <- data4 %>% 
  ggplot(aes(x = Iteracion, y = Valor, colour = Parametro)) + 
  geom_line() +
  theme_bw() +
  xlab('Iteración') + ylab('Valor') +
  geom_vline(xintercept = h, col = 'black') +
  scale_colour_manual(values = c('blue', 'red'), name = 'Parametro') +
  scale_x_continuous(breaks = round(seq(0, max(data4$Iteracion), length.out = 10)), 
                     sec.axis = dup_axis(breaks = NULL, name = NULL)) + 
  theme(legend.position = c(0.2, 0.8))


# Creación de archivo en forma de Arreglo de Grobs

pdf(file = 'Trayectoria_Convergencia.pdf', width = 6, height = 8); {
  gridExtra::grid.arrange(G1,G2)
}; dev.off()

##########################################################################-
# Verificación de distancia entre criterios -------------------------------
##########################################################################-
# Creación de objeto *data5* 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Expansión de variable *Parametro* en columnas
##  2 Calcular la diferencia entre los dos criterios
##  3 Colapsar las columnas en la variable *Parámetro*
##  4 COnvertir la variable *Parametro* en un factor ordenado
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

data5 <- data4 %>% 
  spread(Parametro, Valor) %>% 
  mutate(d = BICc - LRT) %>% 
  gather(Parametro, Valor, -matches("Iteracion")) %>% 
  mutate(Parametro = factor(Parametro, levels = c('BICc', 'LRT', 'd')))

# Gráfico 3 - Muestra diferencias en las trayectorias de convergencia

G3 <- data5 %>% 
  ggplot(aes(x = Iteracion, y = Valor, colour = Parametro)) + 
  geom_line() +
  theme_bw() +
  xlab('Iteración') + ylab('Valor') +
  geom_vline(xintercept = h, col = 'black') + 
  scale_colour_manual(values = c('blue', 'red', 'green4'), name = 'Parametro') +
  scale_x_continuous(breaks = round(seq(0, max(data4$Iteracion), length.out = 10)), 
                     sec.axis = dup_axis(breaks = NULL, name = NULL)) + 
  facet_wrap(Parametro ~ ., ncol = 1, scales = 'free') +
  theme(legend.position = "none")

# Almacenamiento de gráfico G3
ggsave(filename = "Tray_Difer.pdf", plot = G3, device = 'pdf', 
       width = 5, height = 7)  






