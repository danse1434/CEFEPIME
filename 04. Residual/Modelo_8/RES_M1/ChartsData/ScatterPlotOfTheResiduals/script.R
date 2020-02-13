# Carga de paquetes
require(tidyverse)
require(rlang)
require(grid)
require(gridExtra)
require(readxl)
##########################################################################-
# Selección de directorio de trabajo
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '04. Residual', 'Modelo_8', 'RES_M1', 'ChartsData', 
                'ScatterPlotOfTheResiduals'))



# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))

# Apertura de archivo de datos 
y_1_residuals <- read_csv("y_1_residuals.txt")


# Función con batería de pruebas de normalidad
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


# Resultados
normtest_batery(data = y_1_residuals,
                vector = c(4, 7, 13),
                alpha = 0.05)
