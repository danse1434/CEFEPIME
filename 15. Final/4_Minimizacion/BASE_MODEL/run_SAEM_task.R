# Carga de paquetes 
require(tidyverse)
require(lixoftConnectors) # API Monolix
initializeLixoftConnectors(software="monolix")


# Asignación de lambda y ejecución
for (i in 1:27) {
  
  # Ruta al archivo de control
  filename <- file.path(getwd(), 'BASE_MODEL', 'Assessment', paste0('Run',i))
  
  # Carga proyecto con API
  loadProject(projectFile = filename)
  
  # Computar estimación de máxima verosimilitud
  # runLogLikelihoodEstimation(linearization = FALSE, wait = TRUE)
  
  # Computar registros de convergencia
  runConditionalDistributionSampling()
  runConditionalModeEstimation()
  runStandardErrorEstimation()
  
  scenario = getScenario()
  scenario$plotList = c('saemresults')
  setScenario(scenario)
  
  computeChartsData()
}
