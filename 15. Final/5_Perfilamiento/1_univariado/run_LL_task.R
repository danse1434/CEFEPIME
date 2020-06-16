# Carga de paquetes 
require(tidyverse)
require(lixoftConnectors) # API Monolix
initializeLixoftConnectors(software="monolix")

main_dir <- file.path('C:', 'Users', 'Daniel', 'OneDrive', 'Documents', 
                      '(Proyecto)_Estudio_PKPD','CEFEPIME','15. Final',
                      '5_Perfilamiento', '1_univariado')


# Asignación de lambda y ejecución
for (i in 1:100) {
  
  print(paste('Proyecto', i, '=>','Iniciado'))
  # Ruta al archivo de control
  filename <- file.path(main_dir, 'Q_pop', paste0('A',i), 'MODELO_FINAL.mlxtran')
  
  # Carga proyecto con API
  loadProject(projectFile = filename)
  
  # Computar registros de convergencia
  
  scenario = getScenario()
  scenario$plotList = c('saemresults')
  setScenario(scenario)
  
  computeChartsData()
}
