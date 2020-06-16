# Carga de paquetes 
require(glue)
require(tidyverse)
require(lixoftConnectors) # API Monolix
initializeLixoftConnectors(software="monolix")

setwd(file.path(getwd(), 'BASE_MODEL', 'bootstrap', 'nonParametric'))

# Asignación de lambda y ejecución
for (i in 383:385) {
dis = glue('M_Error_bootstrap_bootstrap_{i}.mlxtran')
# Ruta al archivo de control
filename <- file.path(dis)
# Lectura y asignación de archivo de control en envir.
A = readChar(filename, file.info(filename)$size) 
# Modificar el archivo de control
# Se adicionan gráficos con observaciones
B <- A %>% 
  str_replace_all(pattern = "plotResult\\(run\\s\\=\\sfalse,method\\s\\=\\snone\\s\\)", 
              replacement = "plotResult(run = false, method = {saemresults})")

# Sobrescribir el archivo de control
write_lines(B, filename)

# Inicio carpeta
print(paste0('Proyecto ', i,'/1000', ' => Iniciado'))

# Carga proyecto con API
loadProject(projectFile = filename)
# Computar datos de gráficos
computeChartsData()
}

