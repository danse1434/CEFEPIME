# Carga de paquetes 
require(tidyverse)
require(lixoftConnectors) # API Monolix
initializeLixoftConnectors(software="monolix")


# Creación de malla de vector
lamb_vec <- c(seq(-3, -0.1, 0.1),
              seq(+0.1, +3, 0.1)) 

# Asignación de lambda y ejecución
for (i in 1:length(lamb_vec)) {

# Asignación de lambda
lambda <- lamb_vec[[i]]

# Ruta a directorio de trabajo
dis <- file.path(
  'tbs_folder_1',
  paste0('lambda_', format(round(lambda, 2), nsmall = 2)))
# Ruta al archivo de control
filename <- file.path(dis, 'final_model.mlxtran')
# Lectura y asignación de archivo de control en envir.
A = readChar(filename, file.info(filename)$size) 

# Modificar el archivo de control
# Se adicionan gráficos con observaciones
B <-A %>% 
  str_replace_all(pattern = "residualsscatter", 
              replacement = "residualsscatter, outputplot, indfits, obspred")

# Sobrescribir el archivo de control
write_lines(B, filename)

# Carga proyecto con API
loadProject(projectFile = filename)
# Computar datos de gráficos
computeChartsData()
}