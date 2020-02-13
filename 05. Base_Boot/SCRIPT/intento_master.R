# Selección de directorio de trabajo
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '05. Base_Boot', 'SCRIPT'))


source('runMonolix.R', encoding = 'UTF-8')

direct <- file.path('..', 'BOOT', paste0('Data', 2))
direct2 <- normalizePath(direct, winslash = '/')


system('HELP')


\\n
       HELP')

C:/ProgramData/Lixoft/MonolixSuite2019R1/bin/monolix.bat
shell('C:/ProgramData/Lixoft/MonolixSuite2019R1/bin/monolix.bat')

C:\ProgramData\Lixoft\MonolixSuite2019R1\bin


runMonolix(
  name = c(direct2),
  mlxInstallDir = 'C:/ProgramData/Lixoft/MonolixSuite2019R1/bin/',
  saveGraphics = TRUE,
  display = TRUE
)






for (i in 1:10) {
  print(paste0('Modelo ', i, ' terminado'))
}

