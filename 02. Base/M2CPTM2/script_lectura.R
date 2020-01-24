##########################################################################-
# PROYECTO FARMACOCINÉTICA POBLACIONAL DE CEFEPIME EN NF ------------------
##########################################################################-
# Simulación de Concentraciones

# Carga de Paquetes ---------------------------------------------------------
require(tidyverse); theme_set(theme_classic())
require(ggExtra)
require(grid); 
vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
require(plyr)
require(rlang)
require(reshape2)
require(gridExtra)

# Selección de directorio ---------------------------------------------------
setwd(file.path('F:','Documentos','(Proyecto)_Estudio_PKPD','CEFEPIME',
                '02. Base','M2CPTMLLOQ'))
# dir <- dirname(parent.frame(2)$ofile) # frame(3) también funciona.
# setwd(dir)

##########################################################################-
# Lectura de archivo de Iteraciones ---------------------------------------
##########################################################################-
data <- read_csv("./2compclIT.csv",skip = 3)
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  (1) Seleción de nombres de columna
##  (2) Eliminar los carácteres "/"
##  (3) Eliminar los carácteres "'"
##  (4) Eliminar espacios
##  (5) Crear los archivos de datos
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
colnames(data) =
  colnames(data) %>%
  str_replace_all(., "/", "_") %>%
  str_replace_all(., "'", "") %>%
  str_replace_all(., "[:blank:]", "")
write.csv(x=data,file='./Tablas/data_IT.csv')
##########################################################################-
# Creación de mosaico con trayectorias de iteración
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  (1) Creación de expresiones por rlang, y objeto lista en blanco (ls)
##  (2) Hacer regresión con las últimas 10 filas y las variables (x,y)
##  (3) Pasar la función summary a cada objeto tipo lm
##  (4) Extraer la matriz de coeficientes
##  (5) Extraer la matriz completa
##  (6) Extraer el 7. objeto de la matriz que es el valor p de la pendiente
##  (7) Convertir a carácter y luego a número
##  (8) Deslistar y volver un vector
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
conv_asses =
  function(x,y,data,n) {
    if (is_missing(n)) {
      n = 10
    }
    xe = rlang::ensym(x); ye = rlang::ensym(y); ls = list()
    for (i in 1:(dim(data)[1])) {
      if (i < n) {
        ls[[i]] = NA
        
      } else {
        a = rlang::expr(!!ye ~ !!xe + 1)
        ls[[i]] = lm(a, data = data[(i - (n-1)):i, ])
      }
    }
    
    Z = ls %>%
      map( ~ summary(.x))  %>%
      map(., `[`, 'coefficients') %>%
      map(., `[[`, 1) %>%
      map( ~ pluck(.x, 8, .default = 1E-9)) %>%
      unlist(., recursive = F, use.names = F) 
    return(Z)
  }


v = data %>% 
    add_column(new = conv_asses(x = "Iteration",y = `-2LogLikelihood`,data,30)) %>% 
    mutate(new2 = if_else(new>=0.05,T,F)) %>% 
    filter(new2 == T) %>%
    summarise(A = first(Iteration)) %>% 
    map(., ~ pluck(.x,1)) %>% unlist()  


IT_list = list()
for (i in 2:18) {
  y = data %>% as.data.frame(.)
  plot = ggplot(data = y, mapping = aes_string(x = y$Iteration, y =
                                                 y[, i])) +
    geom_line(col = 'blue4') +
    xlab('Iteración') +
    ylab(paste(names(y)[i])) +
    scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
    scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL)) +
    theme(axis.line = element_line(colour = 'black')) + theme_classic() +
    geom_vline(xintercept = v)
  IT_list[[i - 1]] = plot
}

do.call("grid.arrange", c(IT_list, ncol=4))


# Data 1 --------------------------------------------------------------------

data1 <-
  read_csv(
    paste0(dir,"/2compclRSD.csv"),
    col_types = cols(X12 = col_skip())
  )

data1 <-
  data1 %>%
  mutate(IndividID = str_replace_all(IndividID, "[:blank: \']", "")) %>% 
  mutate(IndividID = str_replace_all(IndividID, "(?<=\\_)(?=\\d$)", '0'))

ggplot(data1, aes(x=Obser.Time-72, y = Data, group=IndividID,colour=IndividID))+
  geom_point() +
  geom_line(aes(y=ModelPred.), lty='dashed')  +
  geom_line(aes(y=PopModelPred.), lty='solid') +
  facet_wrap(~IndividID) +
  xlab('TAD')

##########################################################################-
## REVISAR
sigma = sum(data1$Residual^2)/(length(data1$Residual)-(4+2+4)) %>% sqrt()

epsilon_shrik <-
  data1 %>%
  mutate(IWRES = (Data - ModelPred.) / sigma) %>%
  summarise(1 - sd(IWRES))



# Archivo de datos ----------------------------------------------------------
conect = readLines(paste0(dir,"/2compclIND.csv"))
data3 <- read_csv(conect[-6], 
                  skip = 4)

colnames(data3) =
  colnames(data3) %>%
  str_replace_all(., "[/\\-]", "_") %>%
  str_replace_all(., "'", "") %>%
  str_replace_all(., "[:blank:]", "") %>% 
  str_replace_all(., "[\\(\\)]", "") 

data3 <-
  data3  %>%
  select(-X41) %>%
  mutate(IndividID = str_replace_all(IndividID, "[:blank: \']", "")) %>%
  mutate(IndividID = str_replace_all(IndividID, "(?<=\\_)(?=\\d$)", '0'))

data3 %>% 
  ggplot(aes(CLt, Vc)) +
  geom_point() +
  stat_smooth(method = 'lm')
