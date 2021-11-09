
setwd("C:/Users/Daniel/OneDrive/Documents/(Proyecto)_Estudio_PKPD/CEFEPIME/16. Simulaciones/1_Simulacion")

require(data.table)
require(progress)
require(Rcpp)
require(RcppArmadillo)

# Carga de archivo Rcpp
Rcpp::sourceCpp('src/90_verificacion_PTA.cpp')
# Carga de archivo con funciones R
source(file.path('src', '690_paquetesSimulacion.R'), encoding = 'UTF-8')
source(file.path('src', '80_funciones.R'), encoding = 'UTF-8')
source(file.path('src', '102_creacionDataFrame.R'), encoding = 'UTF-8')
source(file.path('src', '820_fun_Simulacion.R'), encoding = 'UTF-8')

# source(file.path('.', 'src', '080_fun_ConversionLog.R'), encoding = 'UTF-8')

set.seed(Sys.time())

# Cálculo de matriz
calculoMatriz <- function(df) {
  cvMat <- conv_matriz(df)  
  PTA1 <- pta_verificador(cvMat, mic_vec, crit = 1.0, treshold = 1)$PTA
  PTA2 <- pta_verificador(cvMat, mic_vec, crit = 0.6, treshold = 1)$PTA
  
  return(data.frame('PTA1' = PTA1, 'PTA2' = PTA2))
}

# N. regímenes
Nreg = 40 # 5000
# N. individuos a simular por régimen
N <- 1e3

regimenes <- tibble(
  DD = runif(Nreg, 2000, 6000),
  II = runif(Nreg, 4, 24),
  SCRMGDL = runif(Nreg, 0.26, 1.02)
) %>% mutate(tinf = map_dbl(II, ~ runif(1, 0.5, .x)))


listaIndicadores <- list()

# Directorio de modelo
# rdir <- file.path('model', 'reference_model.mlxtran')
mdir <- file.path('data', 'modeloFinal.txt')


# Carga de parámetros del modelo (sin incertidumbre)
p_DF <- read_csv(file.path('data', 'populationParameters.txt'))
p <- setNames(p_DF$value, p_DF$parameter)
p['f'] = 0.2

# Liberar memoria de p_DF
rm(p_DF)

# Lista de MIC
mic_vec <- c(1 * (2 ^ (seq(-4, 7, 0.5))))

# Seguimiento de tiempo
pb <- progress_bar$new(format = "[:bar] :current/:total (:percent) eta: (:eta)",
                       total = dim(regimenes)[1])

Sys.time() -> ptm

for (k in 1:dim(regimenes)[1]) {
  
  pb$tick()
  
  # Simulación de parámetros y covariables para cada individuo
  p_DataFrame <- creacionDF(p, N) %>% 
    add_column(SCRMGDL = regimenes$SCRMGDL[k])
  
  out <- list(name = 'ufCc', time = seq(72, 96, length = 30))

  par <- as.vector(as.data.frame(p_DataFrame)[1, ])

  dP_ls <- list()

  for (i in 1:N) {
    # 3.1.1. Selección de parámetros en la fila de cada individuo
    par <- as.vector(as.data.frame(p_DataFrame)[i, ])
    # 3.1.2. Creación de régimen de dosificación para cada individuo
    amountLS = listaTratamiento(regimenes$DD[k],
                                regimenes$II[k],
                                regimenes$tinf[k],
                                16)
    # 3.1.3. Crear para un elemento en la lista de simulación
    dP_ls[[i]] <- list(
      parameter = par,
      treatment = list(amountLS),
      size      = 1,
      level     = 'individual'
    )
  }
  # print(dP_ls)
  # print(mdir)
  # print(out)

  tryCatch(
    expr = {
      res <- simulx(model = mdir,
                    output = out,
                    group = dP_ls)
      
      group_by(res$ufCc, group) %>% 
        nest() %>% 
        mutate(PTA = map(data, ~calculoMatriz(.x))) %>% 
        # select(-data) %>% 
        unnest(PTA) %>% 
        group_by(PTA1.MIC) %>% 
        summarise(PTA1.PTA = mean(PTA1.PTA), 
                  PTA2.PTA = mean(PTA2.PTA)) -> exposureDF
      
      # Liberación de memoria
      # rm(res, dP_ls)
      
      listaIndicadores <- append(listaIndicadores, list(exposureDF))
      rm(exposureDF)
    },
    error = function(e){
      print(e)
      return(NA)
    },
    warning = function(w){
      print(w)
      return(NA)
    })
}

d <- Sys.time() - ptm
d

regimenesDF <- regimenes %>%
  add_column(ind = listaIndicadores) %>%
  unnest(ind) %>%
  rename(MIC = PTA1.MIC, PTA1 = PTA1.PTA, PTA2 = PTA2.PTA)
  

fwrite(regimenesDF, file.path('.', 'results', 'evaluacionRegimenes.csv'), append = T)
