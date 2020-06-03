require(tidyverse)

# Creación de malla de vector
lamb_vec <- c(seq(-3, -0.1, 0.1),
              seq(+0.1, +3, 0.1)) 

# lamb_vec %>% format(., nsmall = 2)

# Creación iterativa de variable lambda y ejecución de source
for (i in 1:length(lamb_vec)) {
  lambda <- lamb_vec[[i]]
  source('generador/gen.R', encoding = 'UTF-8')
}

#-------------------------------------------------------------------------------#

