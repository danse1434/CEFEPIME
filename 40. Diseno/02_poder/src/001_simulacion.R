##------------------------------------------------------------------------------#
## Nombre del Script: Simulación de regímenes de dosificación para Vancomicina
##  
## Propósito del Script: Realizar simulaciones para CLCR no definido
##  
## Autor: Daniel S. Parra González 
## Fecha de creación:  14-03-2021
##  
## Copyright (c) Daniel S. Parra, 2021 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------------#

require(data.table)
require(progress)

source(file.path('src', '069_paquetesSimulacion.R'), encoding = 'UTF-8')

adm <- list(time=0, amount=2000, tinf=0.5)
y   <- list(name='y', time=c(0.5, 1, 2, 4, 6, 8))
p   <- c(Cl_pop = 7.66, V_pop = 17.52, 
         omega_Cl = (8.04 - 7.66) * sqrt(100)/qnorm(0.975), 
         omega_V = (18.60 - 17.52) * sqrt(100)/qnorm(0.975),
         a = 2)

res <- simulx(model = file.path('models', 'reference_model.mlxtran'), 
              parameter = p,
              output = y,
              treatment = adm,
              group = list(size = 1e3, level='individual'))

res$y %>% 
  tibble() %>% 
  ggplot(aes(x=time, y=y, group=time)) +
  geom_boxplot() + coord_cartesian(ylim=c(0,2000))

fwrite(res$y, file.path('results', '001_DF_simulados.csv'))



# source(file.path('src', '080_fun_ConversionLog.R'), encoding = 'UTF-8')
# source(file.path('src', '082_fun_Simulacion.R'), encoding = 'UTF-8')
# source(file.path('src', '002_creacionDataFrame.R'), encoding = 'UTF-8')

# Directorio de modelo
mdir <- file.path('model', 'FINAL.mlxtran')
rdir <- file.path('model', 'reference_model.mlxtran')

# Formato nombre archivos de resultado
outexposure <- '001_exp_desconocido'
outputfile  <- '002_out_desconocido'

#-------------------------------------------------------------------------------#
# 1. Definición parámetros --------------------------
#-------------------------------------------------------------------------------#
set.seed(123456)

# Carga de parámetros del modelo (sin incertidumbre)
p_DF <- read_csv(file.path('model', 'FINAL', 'populationParameters.txt'))
p <- setNames(p_DF$value, p_DF$parameter)

# N. individuos a simular por régimen
N <- 1e3

# Simulación de parámetros y covariables para cada individuo
p_DataFrame <- creacionDF(p, N) %>% 
  add_column(
    # Peso variable
    WTKG = wangScenario3(60, 55, 67, 14) %>%
    {paramMoments(.$x, .$s^2)} %>%
    {rlnorm(N, .$mu, sqrt(.$sigma2))},
    # Aclaramiento de creatinina variable
    CLCRMLMIN = paramMoments(120.33, desvDeInt(120.33, 109.94, 130.72, 14)^2) %>%
    {rlnorm(N, .$mu, sqrt(.$sigma2))}
  )

# Liberar memoria de p_DF
rm(p_DF)

#-------------------------------------------------------------------------------#
# 2. Definición regímenes de dosificación -------------------
#-------------------------------------------------------------------------------#

#' 2.1. Ventana de determinación (output)
#' La exposición se determina en un periodo de 24 horas en estado estacionario 
#' con granularidad de 30 puntos en el intervalo de tiempo. El número de dosis a 
#' simular se ajusta para reflejar el estado estacionario en el tiempo de 
#' medición 72 - 96 h.
#' 
#' Como el intervalo de dosificación más pequeño es 96 se requieren máximo 96/6 
#' dosis en este régimen (y por extensión en los otros) para simular 
#' concentraciones de estado estacionario.
#'  
out <- list(name = 'Cc', time = seq(72, 96, length = 30))

#' 2.2. Regímenes de dosificación
#' En la ruta './data/adm_list.csv' se encuentra una tabla con las especificaciones 
#' de historial de administración, se prueban 21 configuraciones, con 3 dosis 
#' diarias, 4 intervalos de dosificación y 3 tiempos de infusión, eliminando 
#' historiales redundantes.
#' 
admDF <- fread(file.path('data', 'adm_list.csv'))

#' Se crea una lista con sublistas que tienen DD, ii, y tinf.
#' 
adm <- list()
# Crear una lista con información de dosis diaria, II, y tiempo de infusión
for (i in 1:dim(admDF)[1]) {
  row <- admDF[i,]
  adm[[row$ID]] <- list(DD = row$DD, ii = row$ii, tinf = row$tinf)
}

#-------------------------------------------------------------------------------#
# 3. Simulación -----------------------------------------------------
#-------------------------------------------------------------------------------#

# Seguimiento de tiempo
pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: (:eta)",
                       total = length(adm))

#' Se procesa cada regímen de dosificación por separado debido a limitaciones de 
#' RAM
#' 
for (k in 1:length(adm)) {
  
  pb$tick()
  
  # Lista vacía con parámetros y tratamientos para cada individuo a simular
  dP_ls <- list()
  
  #' 3.1. Creación de diseño de estudio
  #' 
  for (i in 1:N) {
    # 3.1.1. Selección de parámetros en la fila de cada individuo
    par <- as.vector(as.data.frame(p_DataFrame)[i, ])
    # 3.1.2. Creación de régimen de dosificación para cada individuo
    # Se multiplica la dosis diaria por el peso, ya que está original en mg/kg/d
    amountLS = listaTratamiento(adm[[k]]$DD * par$WTKG, adm[[k]]$ii, adm[[k]]$tinf, 16)
    # 3.1.3. Crear para un elemento en la lista de simulación
    dP_ls[[i]] <- list(
      parameter = par,
      treatment = list(amountLS),
      size      = 1,
      level     = 'individual'
    )
  }
  
  #' 3.2. Simulación y medición de exposición
  res <- exposure(model = rdir,
                  output = out,
                  group = dP_ls)
  
  #' 3.4. Separación de elementos de exposición *Cc* y perfiles plasmáticos *output$Cc*
  exposureDF <- as.data.table(res$Cc) 
  outputDF   <- as.data.table(res$output$Cc)
  # Liberación de memoria
  rm(res, dP_ls)
  
  #-------------------------------------------------------------------------------#
  #' 3.5. Manipulación de tabla de exposición ------------------------------
  #-------------------------------------------------------------------------------#
  #' Se eliminan cifras significativas que sólo ocupan memoria para *Cmax*, *Cmin*, 
  #' y *AUC*, se eliminan columnas innecesarias *t1*, *t2*, *step*, *tmin*, y *tmax* 
  
  exposureDF[, `:=`(
    t1   = NULL,
    t2   = NULL,
    step = NULL,
    tmin = NULL,
    tmax = NULL,
    auc  = round(auc, 3),
    cmax = round(cmax, 3),
    cmin = round(cmin, 3)
  )] 
  
  #-------------------------------------------------------------------------------#
  #' 3.6. Manipulación de perfiles plasmáticos ------------------------------
  #-------------------------------------------------------------------------------#
  #' Para los perfiles se reduce a dos cifras significativas a *time*, y *Cc*.
  #' 
  outputDF[, `:=`(time = round(time, 2),
                  Cc = round(Cc, 2))]
  #' Crear tabla de resumen con mediana e intervalo de predicción 95% de perfiles
  outputDF1 <- outputDF[, by = time,
                        .(
                          Q1 = quantile(Cc, probs = 0.025),
                          Q2 = median(Cc),
                          Q3 = quantile(Cc, probs = 0.975)
                        )]
  
  # Almacenamiento de simulaciones en archivos con permanencia
  fwrite(exposureDF, file.path('results', paste0(outexposure, k, '.csv')))
  fwrite(outputDF1,  file.path('results', paste0(outputfile, k, '.csv')))
  
  rm(exposureDF, outputDF, outputDF1)
}





