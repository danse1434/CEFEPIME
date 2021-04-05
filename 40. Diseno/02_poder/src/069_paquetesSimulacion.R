# Carga de paquetes
require(tidyverse)
theme_set(theme_bw())

require(gt)
monolix2019R2.path <-  "C:/ProgramData/Lixoft/MonolixSuite2019R2"
require(lixoftConnectors, lib.loc = monolix2019R2.path )
require(mlxR)
initMlxR(path = monolix2019R2.path)

# Número de poblaciones
nPob   <- 50