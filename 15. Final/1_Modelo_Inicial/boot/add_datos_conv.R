# Carga de paquetes 
require(glue)
require(tidyverse)
require(lixoftConnectors) # API Monolix
initializeLixoftConnectors(software="monolix")

# Asignación de lambda y ejecución
for (i in 1:1000) {
# Ruta a directorio de trabajo
dis <- glue('boot/B{i}')
# Ruta al archivo de control
filename <- file.path(dis, 'final_model.mlxtran')
# Lectura y asignación de archivo de control en envir.
A = readChar(filename, file.info(filename)$size) 
# Modificar el archivo de control
# Se adicionan gráficos con observaciones
B <-A %>% 
  str_replace_all(pattern = "populationParameters\\(\\)", 
              replacement = 'populationParameters()
plotResult(run = false, method = {saemresults})')

# Sobrescribir el archivo de control
write_lines(B, filename)

# Carga proyecto con API
loadProject(projectFile = filename)
# Computar datos de gráficos
computeChartsData()
}