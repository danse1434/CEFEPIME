
#-------------------------------------------------------------------------------#
# Correr de forma externa, creando la variable lambda

# Directorio de lambda en texto
lambda_dir = format(round(lambda,2), nsmall=2)

# Leer los datos originales
data <- read_delim('generador/data/1_data_TSFD.csv', ';', na = '.')
# Leer archivo de control de Monolix
fileName <- 'generador/final_model.mlxtran'
A = readChar(fileName, file.info(fileName)$size)

# Crear carpeta de trabajo
dir.create(glue::glue('tbs_folder_1/lambda_{lambda_dir}'))
# Realizar la transformación Box-Cox de DV
data1 <- data %>% 
  mutate(DV_bc = ((DV^lambda)-1)/lambda)
# Escribir nuevo archivo de datos modificados
write_delim(data1, 
            glue::glue('tbs_folder_1/lambda_{lambda_dir}/data_TSFD.csv'),
            na = '.', delim = ';')
# Realizar la transformación de Box-Cox de IV
B <- str_replace_all(A, "lambda", as.character(lambda))
# Escribir nuevo archivo de control de Monolix
write_lines(B, 
            glue::glue('tbs_folder_1/lambda_{lambda_dir}/final_model.mlxtran'), 
                       sep = '\n')
# plotResult(run = false,method = {residualsscatter, residualsdistribution })