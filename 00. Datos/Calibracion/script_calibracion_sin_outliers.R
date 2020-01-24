# Documentos de Calibración -----------------------------------------------
#####################################################################################################-
# En este nuevo script se retira la curva de calibración 'P-2-4_F-051115_A-JCAR' que se considera co-
# mo anómala dentro del análisis de influencia realizado.
#####################################################################################################-
require(tidyverse)
require(lubridate)
require(tibble)
require(chemCal)

setwd(file.path('F:','Documentos','Calibracion'))

data = read_csv(file = './Datos/Calibracion.csv')
source("./Script/funciones.R")

# Se coloca en formato fecha a la columna con fechas
data$Fecha = dmy(data$Fecha) %>% 
  format(., format = "%d%Om%y")
  
# Se arreglan los carácteres, se les quita unas comillas extra al inicio y final
data$Pacientes = data$Pacientes %>% str_sub(., start = 2, end = -2)

# Tabla con las concentraciones por nivel
Conc_data = 
tribble(~"Nivel", ~"Conc.",
  "C1", 50,   "C2", 25,
  "C3", 12.5,   "C4", 6.25,
  "C5", 3.125)
# Se colocan las concentraciones en la tabla data
data1 = inner_join(data, Conc_data, by="Nivel")

# Identificación de posibles datos anómalos -------------------------------------------
  #####################################################################################################-
  # Se realiza una función que permite saber si un dato es anómalo por estar afuera de la distribución
  # en un rango Q1-1.5*IQR a Q3+1.5*IQR (con valores 1.5 veces IQR).
  #####################################################################################################-

out_det <- function(vec, val){
  a1 = quantile(x=vec, probs=0.75)[[1]] + IQR(vec)*1.5
  a2 = quantile(x=vec, probs=0.25)[[1]] - IQR(vec)*1.5
  if(val >= a2 & val <= a1){
    b = TRUE } else { b = FALSE }
  return(b)
}
out_det = Vectorize(FUN=out_det, vectorize.args = 'val') # Vectorizar función
  
  #####################################################################################################-
  # Se crea una nueva base de datos con cambios como: (1) adición de columna *ID* que tiene un código de
  # identificación para cada registro, (2) cambio de la variable *Nivel* de caractéres número, (3) agrupa-
  # ción de variables *Nivel* < *Autor* < *Fecha* < *Pacientes*, (4) adición de columna *Norm_Diam* que 
  # contiene datos de diámetro de halo de inhibición, (5) adición de variable tipo carácter que tiene inf-
  # ormación de *Pacientes*, *Fecha*, y *Autor* en *Grupo*, (6) adición de variable lógica (*out_status*) 
  # que indica si valor se encuentra dentro de rango especificado de 1.5*IQR, para cada grupo formado.
  #####################################################################################################-
  data2 = data1 %>% 
  rownames_to_column(.data = ., var = 'ID') %>% 
  mutate(Nivel = as.factor(Nivel)) %>% 
  group_by(Pacientes, Fecha, Autor, Nivel) %>% 
  mutate(Norm_Diam = (Diametro-mean(Diametro))/sd(Diametro)) %>% 
  mutate(Grupo = paste0('P-',Pacientes,"-F-",Fecha,"-A-",Autor)) %>% 
  mutate(out_status = out_det(vec=Diametro, val=Diametro))  %>% 
  filter(Grupo != 'P-2-4-F-051115-A-JCAR')
  
  #####################################################################################################-
  # Se elabora un gráfico *g1_plot* que contiene los datos outlier graficados junto a la distribución de 
  # los datos ilustrada mediante un boxplot. 
  # NÓTESE que el nivel de concentración C3 es el que posee la mayoría de determinaciones, porque se mi-
  # dió dos veces en cada ensayo.
  #####################################################################################################-
  g1_plot = 
  ggplot(data2, aes(x=Nivel, y = Norm_Diam)) +
  geom_boxplot(aes(col=Grupo, fill=Grupo), alpha=0.05,
               outlier.color = NULL, outlier.fill = NULL) +
  geom_text(data = filter(data2, out_status == FALSE),
            aes(x=Nivel, y = Norm_Diam, label=ID), 
            hjust=0, vjust=0,
            inherit.aes = FALSE) +
  facet_wrap(~Grupo) +
  xlab('Nivel de Concentración') +
  ylab('Diámetro de Halo, estandarizado por grupo') +
  coord_cartesian(ylim=c(-2.5,2.5)) +
  theme_bw() +
  geom_point(aes(col=Grupo), shape=1); print(g1_plot)
  # Almacenamiento de gráfico
  ggsave(filename = './Resultados/Figuras/21_Outliers_identificados.pdf', plot = g1_plot, 
         device = 'pdf', width = 8.27, height = 5.83)
  #####################################################################################################-

# Análisis de Exactitud ---------------------------------------------------
  #####################################################################################################-
  # Se realiza una comprobación del modelo de regresión lineal con mejor ajuste para la curva, se parte
  # del data frame *data2* que contiene toda la información relevante sin agrupación. El modelo de regre-
  # sión debe hacerse con la relación entre *Diametro* ~ log(*Conc.*). Dentro de cada grupo se debe hacer
  # un resumen de los datos de repetición, y hacerse la regresión con los datos de todos los grupos.
  # 
  # Se elabora el data frame *D3_COUT* que tiene los datos promedio (*m_Diam*) de las repeticiones en cada 
  # grupo, así como desviación estándar de Diámetro intra-grupo (*s_Diam*), y concentración promedio 
  # (*m_Conc*). Este nuevo data frame *D3_COUT* ya contiene a la variable concentración transformada a log-
  # aritmo, y esta agrupada por nivel.
  #####################################################################################################-
  D3_COUT = data2 %>% 
    summarise(m_Conc = mean(Conc.),
              m_Diam = mean(Diametro),
              s_Diam = sd(Diametro)) %>% 
    ungroup() %>% 
    mutate(log_conc = log(m_Conc)) %>% 
    group_by(Nivel) 
  #####################################################################################################-
  # Se elabora un gráfico con los datos en el modelo que se pretende obtener, se obtiene una amplia varia-
  # ción en los valores de respuesta con los datos concentración respuesta. Se tiene que P-2-4-F-201115-A-
  # JCAR se diferencia de los otros grupos por tener una menor sensibilidad a cambios en la concentración
  # que el resto. 
  #####################################################################################################-
  # A continuación, se calculan intervalos de predicción y confianza para una curva de regresión lineal 
  # con ponderación de tipo 1/log(conc).
  #
  bands = data.frame('log_conc' = seq(0,4,0.01))
  
  bands[,'log_conc'] = seq(0,4,0.01)
  bands[,'conf'] = predict(lm(m_Diam ~ log_conc, D3_COUT, weights=1/m_Diam), 
                           newdata = data.frame(log_conc = seq(0,4,0.01)),
                           interval = 'confidence')[,1]
  bands[,'conf_LI'] = predict(lm(m_Diam ~ log_conc, D3_COUT, weights=1/m_Diam), 
                               newdata = data.frame(log_conc = seq(0,4,0.01)),
                               interval = 'confidence', weights = 1/(bands$conf))[,2]
  bands[,'conf_LS'] = predict(lm(m_Diam ~ log_conc, D3_COUT, weights=1/m_Diam), 
                              newdata = data.frame(log_conc = seq(0,4,0.01)),
                              interval = 'confidence', weights = 1/(bands$conf))[,3]
  bands[,'pred_LI'] = predict(lm(m_Diam ~ log_conc, D3_COUT, weights=1/m_Diam), 
                               newdata = data.frame(log_conc = seq(0,4,0.01)),
                               interval = 'prediction', weights = 1/(bands$conf))[,2]
  bands[,'pred_LS'] = predict(lm(m_Diam ~ log_conc, D3_COUT, weights=1/m_Diam), 
                              newdata = data.frame(log_conc = seq(0,4,0.01)),
                              interval = 'prediction', weights = 1/(bands$conf))[,3]
  #####################################################################################################-
  # Estos intervalos se adicionan a una gráfica en conjunto con los datos.
  #####################################################################################################-
  g2_plot = 
  D3_COUT %>% 
    ungroup() %>% 
    mutate(Grupo = paste0('P-',Pacientes,"-F-",Fecha,"-A-",Autor)) %>% 
    group_by(Nivel) %>% 
    ggplot(., aes(x=log_conc, y=m_Diam, group=Grupo, 
                  col=Grupo, shape=Grupo)) +
    # Bandas de confianza y predicción
    geom_ribbon(data = bands, aes(x = log_conc, ymin = pred_LI, ymax=pred_LS), inherit.aes=F,
                fill='blue1',alpha=0.1)+
    geom_ribbon(data = bands, aes(x = log_conc, ymin = conf_LI, ymax=conf_LS), inherit.aes=F,
                fill='blue2',alpha=0.5)+
    geom_line(data = bands, aes(x = log_conc, y=conf), inherit.aes=F,col='blue4',lty='solid')+
    # Puntos
    geom_point() +
    theme_bw() +
    coord_cartesian(ylim = c(15,28), xlim=c(1.13,3.95))+
    ylab(expression(paste('Diámetro ',phi, ' de halo (mm)'))) +
    xlab(expression(paste('ln(',C[nom],')'))) +
    scale_x_continuous(breaks = seq(1,4,0.5)) +
    scale_y_continuous(breaks = seq(15,30,1)) +
  theme(panel.grid = element_line(colour = 'gray98'), 
          legend.spacing.y = unit(0.05,'cm'),
          legend.key.size = unit(0.25,'cm'), legend.key.width = unit(1,'cm'),
          legend.text =  element_text(size=6),
          legend.title = element_text(size=8, face = 'bold'),
          legend.position = c(0.8,0.15)); g2_plot
  ggsave(filename = './Resultados/Figuras/22_curva_calibracion_1.pdf', plot = g2_plot, 
         device = 'pdf', width = 5.13, height = 4.29)
  #####################################################################################################-
  # Se realiza una evaluación de diversos  modelos de regresión lineal que se diferencian por la función
  # de ponderación. No se realizan con ajuste del intercepto por el origen ya que los datos no fueron nor-
  # malizados frente a un blanco. 
  #####################################################################################################-
  # Se crea un objeto de tipo lista
  list_reg = list()
  # Se escribe dentro de la lista un grupo de modelos de regresión que varían en la ponderación
  list_reg[[1]] <- lm(formula = m_Diam ~ log_conc,data = D3_COUT)
  list_reg[[2]] <- lm(formula = m_Diam ~ log_conc,data = D3_COUT, weights = 1/m_Diam)
  list_reg[[3]] <- lm(formula = m_Diam ~ log_conc,data = D3_COUT, weights = 1/(m_Diam^2))
  list_reg[[4]] <- lm(formula = m_Diam ~ log_conc,data = D3_COUT, weights = 1/log_conc)
  list_reg[[5]] <- lm(formula = m_Diam ~ log_conc,data = D3_COUT, weights = 1/(log_conc^2))
  # Se crea una nueva lista con los resumenes de regresión
  list_reg_summ = lapply(list_reg, function(x) summary(x))
  # Se le mide la extensión a la lista creada
  l = length(list_reg_summ)
  #####################################################################################################-
  # Se hace un resumen de los datos, primero con los coeficientes de los modelos de regresión lineal es-
  # tudiados. 
  reg_sum =
  map(.x = list_reg_summ, .f = "coefficients") %>% 
    unlist(.) %>% 
    matrix(., byrow = T, ncol = 8)
  
  colnames(reg_sum) = rep(c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)'),each = 2)
  row.names(reg_sum) = c('OLS', 'WLS (x^-1)','WLS (x^-2)','WLS (y^-1)','WLS (y^-2)')
  # 
  lmp <- function (modelo) {
    # Función para obtener el valor p de los estadísticos
    if (class(modelo) != "lm") stop("No es un objeto de clase 'lm' ")
    f <- summary(modelo)$fstatistic
    p <- pf(f[1],f[2],f[3],lower.tail=F)
    attributes(p) <- NULL
    return(p)  }
  # Se hace un resumen de los parámetros de bondad de ajuste del modelo - - - - - 
  REG = matrix(ncol = 8, nrow = l)
  for (i in 1:l) {
    REG[i,2] = list_reg_summ[[i]]$sigma # Sigma
    REG[i,3] = list_reg_summ[[i]]$df[2] # Grados de Libertad
    REG[i,4] = list_reg_summ[[i]]$r.squared
    REG[i,5] = list_reg_summ[[i]]$adj.r.squared
    REG[i,6] = list_reg_summ[[i]]$fstatistic[["value"]]
    REG[i,7] = qf(.95,list_reg_summ[[i]]$df[2],1)
    REG[i,8] = lmp(list_reg[[i]])
  }
  REG = data.frame(REG)
  colnames(REG) <- c("Tipo","RSE","gL","r2","r2.adj","f_calc","f_tab","p")
  # Se encuentra que todos los modelos describen de manera adecuada el modelo y las pruebas de bondad de 
  # ajuste cumplen. El valor de r2 es algo bajo para las expectativas en una validación de método bio-
  # analítico.
  # En base a los datos de los modelos de regresión se debería escoger el modelo N.º 4 que tiene pondera-
  # ción de tipo (1/log_conc). 
  #####################################################################################################-
  # Se realizan estimaciones de valores retro-predichos para cada modelo de regresión utilizando la li-
  # brería chemCal de R. Las estimaciones con valores en el dominio de concentración, se almacenan los 
  # datos en la lista *L.df*.
  L.df = list_reg %>% 
    lapply(X = ., FUN = function(x){retro.predic(model = x,
                                                 data = D3_COUT)[[2]]})
  # Se unen los datos de todos los modelos en el data frame *L1.df*
  L1.df = plyr::rbind.fill(L.df[[1]],L.df[[2]],L.df[[3]],L.df[[4]],L.df[[5]])
  # Se crea la variable modelo que indica el modelo de origen
  L1.df["Modelo"] = as.factor(rep(1:5,each = dim(D3_COUT)[1]))
  # Se traen las variables contenidas en *D3_COUT*
  L1.df["Pacientes"] = as.factor(rep(D3_COUT$Pacientes,5))
  L1.df["Fecha"] = as.factor(rep(D3_COUT$Fecha,5))
  L1.df["Autor"] = as.factor(rep(D3_COUT$Autor,5))
  L1.df["Nivel"] = as.factor(rep(D3_COUT$Nivel,5))
  L1.df["m_Conc"] = as.double(rep(D3_COUT$m_Conc,5))
  L1.df["m_Diam"] = as.double(rep(D3_COUT$m_Diam,5))
  # Se reorganiza la tabla
  L1.df =
  L1.df %>% select("Pacientes","Fecha","Autor","Nivel","Modelo","m_Conc",
                   "Pred.","SE","Conf.","LI","LS","Modelo","m_Diam") %>% 
    mutate(EA = Pred.-m_Conc) %>% 
    mutate(ER = EA/m_Conc) %>% 
    mutate(Grupo = paste0('P-',Pacientes,"-F-",Fecha,"-A-",Autor)) %>% 
    rownames_to_column(., var = 'ID')
  #####################################################################################################-
  # GRÁFICOS DE LINEALIDAD - MODELOS DE REGRESIÓN
  #####################################################################################################-
  # Se elabora un gráfico que compara los datos de cada corrida de muestras frente a la predicción del 
  # modelo escogido N.4.
  g3_plot = 
  ggplot(data=filter(L1.df, Modelo==2), aes(x = m_Conc, y = Pred., 
                                            group=Grupo, col=Grupo, shape=Grupo)) + 
    # facet_wrap(~Grupo, labeller = label_value) + 
    geom_abline(slope = 1, intercept = 0, lty='dashed') +
    theme_bw() +
    stat_smooth(method = 'lm', se = FALSE, size=0.25, aes(lty=Grupo)) +
    geom_point(size=2) + 
    xlab('Concentración nominal (mg/L)') +
    ylab('Concentración predicha (mg/L)') +
    coord_cartesian(xlim=c(0,60), ylim=c(0,60)) +
    theme(panel.grid = element_line(colour = 'gray98'), 
          legend.spacing.y = unit(0.05,'cm'),
          legend.key.size = unit(0.25,'cm'), legend.key.width = unit(1,'cm'),
          legend.text =  element_text(size=6),
          legend.title = element_text(size=8, face = 'bold'),
          legend.position = c(0.2,0.7)); g3_plot
  ggsave(filename = './Resultados/Figuras/23_Comparacion_CC_1.pdf', plot = g3_plot, 
         device = 'pdf', width = 5.13, height = 4.29)  
  # Se elabora un gráfico que compara los modelos frente a los datos de todas las corridas PRED vs CONC.
  g4_plot=
  ggplot(data=L1.df, aes(x = m_Conc, y = Pred.)) + 
    geom_point(size=2) + 
    facet_wrap(~Modelo, labeller = label_both) + 
    geom_abline(slope = 1, intercept = 0, lty='solid') +
    stat_smooth(method = 'loess', se = FALSE, lty='dashed') +
    xlab('Concentración nominal (mg/L)') +
    ylab('Concentración predicha (mg/L)') +
    theme_bw() +
    theme(panel.grid = element_line(colour = 'gray98')); g4_plot
  ggsave(filename = './Resultados/Figuras/24_Comparacion_CC_2.pdf', plot = g4_plot, 
         device = 'pdf', width = 5.13*1.2, height = 4.29)  
  # Se elabora un gráfico de errores relativos vs concentración nominal de los datos 
  g5_plot = 
  ggplot(data=filter(L1.df, Modelo==2), 
         aes(x = m_Conc, y = ER)) + 
    geom_point(shape=1,col='blue4') +
    geom_hline(yintercept = 0.0, size=0.5)+
    geom_hline(yintercept = c(-0.20,+0.20), size=0.5, lty='dotted')+
    geom_hline(yintercept = c(-0.15,+0.15), size=0.5, lty='dashed')+
    stat_smooth(method = 'loess', se = TRUE, lty='dashed', 
                fill='blue1', alpha=0.1) +
    theme_bw() +
    coord_cartesian(ylim=c(-0.45,0.45)) +
    xlab('Concentración nominal (mg/L)') +
    ylab('Error relativo') + 
    geom_text(data = filter(.data = L1.df, ((ER>=0.15)|(ER<=-0.15))&(Modelo==4)), 
              aes(x=m_Conc, y=ER, label = ID),
              position = position_jitter(width = 0.1, height = 0.1, seed = 1323),
              inherit.aes=F, col='black',size=4) +
    scale_y_continuous(breaks = seq(-0.5,0.5,0.1),
                       minor_breaks = seq(-0.5,0.5,0.1)) +
    theme(panel.grid = element_line(colour = 'gray98'))
  ggsave(filename = './Resultados/Figuras/25_Grafico_error_concentracion.pdf', plot = g5_plot, 
         device = 'pdf', width = 5.13*1.2, height = 4.29)  
  # Se elabora un gráfico de errores relativos vs diámetro de halo de inhibición de los datos
  g6_plot =
  ggplot(data=filter(L1.df, Modelo==2), 
         aes(x = m_Diam, y = ER)) + 
    geom_point(shape=1,col='red4') +
    geom_hline(yintercept = 0.0, size=0.5)+
    geom_hline(yintercept = c(-0.20,+0.20), size=0.5, lty='dotted')+
    geom_hline(yintercept = c(-0.15,+0.15), size=0.5, lty='dashed')+
    stat_smooth(method = 'loess', se = TRUE, lty='dashed', col='red4', 
                fill='red1', alpha=0.1) +
    theme_bw() +
    geom_text(data = filter(.data = L1.df, ((ER>=0.15)|(ER<=-0.15))&(Modelo==4)), 
              aes(x=m_Diam, y=ER, label = ID),
              position = position_jitter(width = 0.1, height = 0.05, seed = 1323),
              inherit.aes=F, col='red',size=4) +
    coord_cartesian(ylim=c(-0.45,0.45)) +
    xlab('Diám. halo de inhibición (mm)') +
    ylab('Error relativo') +
    scale_y_continuous(breaks = seq(-0.5,0.5,0.1),
                       minor_breaks = seq(-0.5,0.5,0.1)) +
    theme(panel.grid = element_line(colour = 'gray98'))
  
  ggsave(filename = './Resultados/Figuras/26_Grafico_error_diametro_halo.pdf', plot=g6_plot, 
         device = 'pdf', width = 5.13*1.2, height = 4.29)  
  # Gráfico de exactitud por nivel de concentración
  g7_plot =   
  L1.df %>% 
    filter(Modelo == 2) %>% 
    group_by(Nivel) %>% 
    summarise(mn=mean(Pred.), nom=mean(m_Conc)) %>% 
    mutate(Error = (nom-mn)/nom) %>% 
    ggplot(data = ., mapping = aes(x=Nivel, y=Error)) +
    geom_bar(col='blue4', fill='blue4', alpha=0.5,
             stat = 'identity') +
    geom_hline(yintercept = c(-0.20,+0.20), size=0.5, lty='dotted')+
    geom_hline(yintercept = c(-0.15,+0.15), size=0.5, lty='dashed')+
    geom_hline(yintercept = c(0), size=0.5, lty='solid')+
    theme_bw() + xlab('Nivel de concentración') + 
    geom_text(aes(y=Error*1.3,label=round(Error,2)))+
    ylab(expression(paste('Error Relativo ',(C[pred]-C[nom])/C[nom])))  +
    theme(panel.grid = element_line(colour = 'gray98')); g7_plot
  ggsave(filename = './Resultados/Figuras/27_Exactitud_metodo.pdf', plot=g7_plot, 
         device = 'pdf', width = 4.29, height = 4.29)    
    # Gráfico de precisión por nivel de concentración
  g8_plot = 
  L1.df %>% 
    filter(Modelo == 2) %>% 
    group_by(Nivel) %>% 
    summarise(mn=mean(Pred.), sd=sd(Pred.)) %>% 
    mutate(RSD = sd/mn) %>% 
    # Gráfico GGPLOT
    ggplot(data = ., mapping = aes(x=Nivel, y=RSD)) +
    geom_bar(col='green4', fill='green4', alpha=0.5,
             stat = 'identity') +
    geom_hline(yintercept = c(+0.20), size=0.5, lty='dotted')+
    geom_hline(yintercept = c(+0.15), size=0.5, lty='dashed')+
    geom_text(aes(y=RSD*1.3,label=round(RSD,2)))+
    theme_bw() + xlab('Nivel de concentración') + 
    ylab('RSD')  +
    theme(panel.grid = element_line(colour = 'gray98'))
  ggsave(filename = './Resultados/Figuras/28_Repetibilidad_metodo.pdf', plot=g8_plot, 
         device = 'pdf', width = 4.29, height = 4.29)  
#####################################################################################################-
# Pruebas de normalidad
# 
shapiro.test(x = residuals(object = list_reg[[2]]))
require(car)


require(olsrr)
ols_test_breusch_pagan(list_reg[[2]])


ols_plot_cooksd_bar(list_reg[[2]])
ols_plot_cooksd_chart(list_reg[[2]])
ols_plot_dffits(list_reg[[2]])
ols_plot_hadi(list_reg[[2]])
  #####################################################################################################-
  # Repetibilidad #
  # Se realiza una estimación de la precisión en condiciones de repetibilidad teniendo en cuenta los da-
  # tos provenientes de repeticiones de datos. 
  data1 = 
  data1 %>% 
    filter((Pacientes != "2-4")|(Fecha!="051115")|(Autor!="JCAR")) %>% 
    mutate(log_conc = log(Conc.)) %>% 
    rename(m_Diam = Diametro) %>% 
    mutate(Grupo = paste0('P-',Pacientes,"-F-",Fecha,"-A-",Autor)) %>% 
    mutate(Grupo=as.factor(Grupo)) 
  # Predicción de Concentración
  data1[,"Pred."] = retro.predic(model = list_reg[[2]], data = data1, weights = 1/(data1$m_Diam))[[2]][,1]
  data1[,"Pred_LI"] = retro.predic(model = list_reg[[2]], data = data1, weights = 1/(data1$m_Diam))[[2]][,4]
  data1[,"Pred_LS"] = retro.predic(model = list_reg[[2]], data = data1, weights = 1/(data1$m_Diam))[[2]][,5]
  # A continuación, se muestran los resultados de un análisis de los datos de forma longitudinal, con un 
  # ANOVA de medidas repetidas. 
  # El grupo está en el estrato 1, y las repeticiones en el estrato 2
  AOV_list = list()
  AOV_list[[1]] = aov(Pred. ~ Grupo, 
                      data = subset(data1, subset = (Nivel=="C5")))
  AOV_list[[2]] = aov(Pred. ~ Grupo, 
                      data = subset(data1, subset = (Nivel=="C4")))
  AOV_list[[3]] = aov(Pred. ~ Grupo, 
                      data = subset(data1, subset = (Nivel=="C3")))
  AOV_list[[4]] = aov(Pred. ~ Grupo,
                      data = subset(data1, subset = (Nivel=="C2")))
  AOV_list[[5]] = aov(Pred. ~ Grupo, 
                      data = subset(data1, subset = (Nivel=="C1")))
  # 
  AOV_reg_summ = lapply(AOV_list, function(x) summary(x))
  
  AOV_reg = 
    lapply(AOV_reg_summ, `[[`, 1) %>% 
    unlist() %>% 
    matrix(., byrow = T, ncol = 5)
  
  colnames(reg_sum) = rep(c('Estimate', 'Std. Error', 't value', 'Pr(>|t|)'),each = 2)
  row.names(reg_sum) = c('OLS', 'WLS (x^-1)','WLS (x^-2)','WLS (y^-1)','WLS (y^-2)')
  
  data1 %>% group_by(Nivel) %>% summarise(mn = mean(Pred.))
  #####################################################################################################-
  # Límite de Cuantificación
  LOQ = loq(object = list_reg[[2]], alpha=0.05, n = 6, w.loq=log(3.125))
  exp(LOQ[[1]])
  
  exp(lod(object = list_reg[[1]])[[1]])
  exp(loq(object = list_reg[[1]])[[1]])
  
  list_reg[[2]] %>% summary(.)
  LOQ_2 = exp(10*(0.18333/3.82758))
  #####################################################################################################-
  #####################################################################################################-
  #####################################################################################################-
  #####################################################################################################-