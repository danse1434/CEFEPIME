##------------------------------------------------------------------------------#
## Nombre del Script: Análisis de comportamiento de algoritmo Simplex 
## en estimación de tiempos de muestreo 
##  
## Propósito del Script: En este análisis se muestra la trayectoria de un 
## algoritmo simplex en la estimación de diseño de estudio óptimo mediante 
## PFIM basado en el diseño de cefepime.
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  03-04-2021
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#
require(tidyverse)
require(gganimate)
require(plotly)

source(file.path('src', '400_fun_acumular.R'), encoding = 'UTF-8')

plasmaDF <-
  read_csv2(file.path('data', '002_perfilPlasma1.csv'),
            col_names = c('t', 'C'))

#-------------------------------------------------------------------------------#
# 1. Manipulación de datos -----------------------------------------------------
#-------------------------------------------------------------------------------#

df_ls <- read_lines(file.path('01', '01_outFIM.txt'))

# Se buscan sólo lineas con resultados numéricos a los que precede "[1]"
df_ls1 <- df_ls %>% 
  {`[`(., str_detect(., "^\\[1\\]"))} 

# Se filtran las líneas que contienen números
df_ls1a <- df_ls1 %>% 
  `[`(., str_detect(., "\\d\\.\\d+\\s"))

df_ls1b <- df_ls1a %>% 
  map_df(~tibble(.x), .id = 'A') %>% 
  mutate(`.x` = str_replace(.x, '^\\[1\\]\\s', ''),
         `.x` = str_split(.x, '\\s')) %>% 
  unnest(.x) %>% 
  mutate(x = as.double(.x),
         P = rep(1:6 %>% paste0('P', .), 292)) %>% 
  pivot_wider(id_cols = A, names_from = P, values_from = x) %>% 
  mutate(A = as.numeric(A))

# Colocar sólo iteraciones tras inicio de algoritmo Simplex
lec <- 0
df_ls1d <- list()

for (i in 1:length(df_ls)) {
  d <- df_ls[i]
  if (str_detect(d, 'entering fun.amoeba')) {
    lec = 1
  }
  
  if (lec == 1) {
    if (str_detect(d, "^crit")) {
      df_ls1d <- append(df_ls1d, d)
    }
  }
}

df_ls1d <- df_ls1d %>%
  str_replace(., '^crit\\=\\s', '') %>%
  str_trim() %>%
  tibble(crit = .) %>%
  rownames_to_column('count') %>%
  mutate(across(c(count, crit), as.numeric))
  
# `[`(., str_detect(., "^\\[1\\]\\s\\d{1,2}\\.\\d+$"))

df_ls1b <- df_ls1b %>% 
  left_join(df_ls1d, by = c('A' = 'count'))

# Escribir csv
write_csv(df_ls1b, file.path('data', '001_iteraciones_Simplex_FEP.csv'))

df_ls1b1 <- df_ls1b %>% 
  mutate(iteracion = A) %>% 
  accumulate_by(~iteracion)

#-------------------------------------------------------------------------------#
# 2. Gráfico con diseño experimental ----------------------------
#-------------------------------------------------------------------------------#
fig <- df_ls1b %>% 
  mutate(frame = as.double(A)) %>% 
  plot_ly(data = ., type = 'scatter', mode = 'lines') %>% 
  add_segments(x = ~ P1, xend = ~ P1, y = 0, yend = 32, 
               frame = ~ frame, name = 'P1', line=list(dash='dot'))  %>% 
  add_segments(x = ~ P2, xend = ~ P2, y = 0, yend = 32, 
               frame = ~ frame, name = 'P2', line=list(dash='dot'))  %>% 
  add_segments(x = ~ P3, xend = ~ P3, y = 0, yend = 32, 
               frame = ~ frame, name = 'P3', line=list(dash='dot'))  %>% 
  add_segments(x = ~ P4, xend = ~ P4, y = 0, yend = 32, 
               frame = ~ frame, name = 'P4', line=list(dash='dot'))  %>% 
  add_segments(x = ~ P5, xend = ~ P5, y = 0, yend = 32, 
               frame = ~ frame, name = 'P5', line=list(dash='dot'))  %>% 
  add_segments(x = ~ P6, xend = ~ P6, y = 0, yend = 32, 
               frame = ~ frame, name = 'P6', line=list(dash='dot')) %>% 
  add_lines(x = ~t, y = ~C, data=plasmaDF, mode='lines', type='scatter',
            name='Cp', inherit = FALSE) %>% 
  layout(
    xaxis = list(title = 'Tiempo (hrs)', range = c(0, 8)),
    yaxis = list(title = 'Conc. plasmática (mg/L)')
  ) %>% 
  animation_opts(
    frame = 10, 
    transition = 0, 
    redraw = FALSE
  )

htmlwidgets::saveWidget(
  fig,
  file.path('figures', '001_diseño_estudio.html') %>% normalizePath(),
  selfcontained = T,
  title = 'Diseño Estudio PFIM'
)

#-------------------------------------------------------------------------------#
# 3. Gráfico convergencia-----------------------------------------------------
#-------------------------------------------------------------------------------#

fig1 <- df_ls1b1 %>%
  plot_ly(x = ~A, y=~crit, frame = ~frame,
          type = 'scatter', mode = 'lines',  name = 'Conv',
          line = list(simplyfy = F)) %>% 
  layout(
    xaxis = list(title = 'Iteración', range = c(0, 300)),
    yaxis = list(title = 'Criterio', range = c(6E-3, 0.0111))
  ) %>% 
  animation_opts(
    frame = 10, 
    transition = 0, 
    redraw = FALSE
  )

htmlwidgets::saveWidget(
  fig1,
  file.path('figures', '002_diseño_convergencia.html') %>% normalizePath(),
  selfcontained = T,
  title = 'Convergencia Diseño Estudio PFIM'
)

#-------------------------------------------------------------------------------#
# 4. Figura completa -----------------------------------------------------
#-------------------------------------------------------------------------------#
axx <- list(
  gridcolor='rgb(255, 255, 255)',
  zerolinecolor='rgb(255, 255, 255)',
  showbackground=TRUE,
  backgroundcolor='rgb(230, 230,230)'
)

figTotal <- subplot(fig, fig1, nrows = 1, margin = 0.1, titleX = T, titleY = T) %>% 
  animation_opts(
    frame = 20, 
    transition = 0, 
    redraw = FALSE
  ) %>% 
  layout(title = "<b>Optimización Diseño Muestral Simplex</b><br>",
         scene = list(
           # domain=list(x=c(0,0.5),y=c(0.0,1)),
                      xaxis=append(axx, c(title = 'Tiempo (hrs)', range = c(0, 8))), 
                      yaxis=append(axx, c(title = 'Conc. plasmática (mg/L)')),
                      aspectmode='cube'),
         scene2 = list(
           # domain=list(x=c(0.5,1),y=c(0.0,1)),
                       xaxis=append(axx, c(title = 'Iteración', range = c(0, 300))), 
                       yaxis=append(axx, c(title = 'Criterio', range = c(6E-3, 0.0111))),
                       aspectmode='cube')
  )

htmlwidgets::saveWidget(
  figTotal,
  file.path('figures', '003_diseño_muestral.html') %>% normalizePath(),
  selfcontained = T,
  title = 'Convergencia Diseño Estudio PFIM'
)
