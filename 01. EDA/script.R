require('tidyverse')
require('lubridate')
require('rlang')
require('grid')
require('PerformanceAnalytics')

setwd(dir = file.path('F:','Documentos','(Proyecto)_Estudio_PKPD','CEFEPIME','01. EDA'))
theme_set(theme_bw())

data = read_csv(file = 'DATA/NONMEM_Cefepime.csv', 
                col_names = T, 
                na = c('.'),
                col_types = cols(EVID = col_factor(),
                                 MDV = col_factor(),
                                 SS = col_factor(),
                                 SEXF = col_factor(), 
                                 ANTU = col_factor(),
                                 LLP = col_factor(),
                                 LMP = col_factor()))
# data[,1:30][data[,1:30]=='.'] <- NA # Para reemplazar todos los datos por un valor específico

eGFR = function(Scr,edad,sexo,raza){
  if (sexo == 1) {k = 0.7; b = 1; a = -0.329 # Mujeres
  } else {k = 0.9; b = 0; a = -0.411} # Hombres
  if (raza == 1){c = 1
  } else {c=0}
  # Tasa de filtración glomerular
  GFR = 141*(min(Scr/k,1))^(a)*(max(Scr/k,1))^-1.209*(0.993^edad)*(1.018^b)*(1.159^c)
  return(GFR)
}
eGFR = Vectorize(FUN=eGFR, vectorize.args=c('Scr','edad','sexo'))

eGFR(Scr = 0.46, edad = 43, sexo = 1, raza = 0)

data = data %>% 
  mutate(DATE = dmy(DATE)) %>% 
  mutate(SCM2 = 0.007184*(WTKG^0.425)*(HCM^0.725)) %>% 
  mutate(DATEHOUR = dmy_hms(paste(data$DATE, data$TIME))) %>% 
  mutate(Dx = case_when(LLP == 1 ~ 'LLA',
                        LMP == 1 ~ 'LMA/LMC',
                        (LLP==0 & LMP==0) ~ 'Otros')) %>% 
  mutate(Dx = as.factor(Dx)) %>% 
  mutate(eGFR = eGFR(SCRMGDL, AGEA,SEXF,0))

LOQ = 25/2/2/2/2 # Concentración que se puede cuantificar como mínino

# Pruebas de normalidad ---------------------------------------------------
normtest_batery = function(data, vector, alpha){
  df = matrix(nrow = length(unique(vector)), ncol = 7)
  for (j in vector) {
    X = dplyr::pull(data,j) # Selecciona como un vector atómico a una columna 
    i = match(j,vector) # Encuentra la posición en el vectror
    df[i,1] = colnames(data[,j])
    df[i,2] = ifelse(shapiro.test(X)$p.value < alpha, '+', '-') # Shapiro-Wilks
    df[i,3] = ifelse(nortest::ad.test(X)$p.value < alpha, '+', '-') #Anderson-Darling
    df[i,4] = ifelse(nortest::cvm.test(X)$p.value < alpha, '+', '-') #Cramer von Mises
    df[i,5] = ifelse(nortest::lillie.test(X)$p.value < alpha, '+', '-') #Liliefors
    df[i,6] = ifelse(nortest::pearson.test(X)$p.value < alpha, '+', '-') #Pearson
    df[i,7] = ifelse(nortest::sf.test(X)$p.value < alpha, '+', '-') #Shapiro Francia
  }
  colnames(df) = c('Variable','Shapiro', 'Anderson_Darling', 'Cramer_von_Mises',
                   'Liliefors','Pearson','Shapiro_Francia')
  return(df)
}

#colnames(data) %>% View(.) #Para ver todos los nombres de las columnas

  data %>% 
    filter(TAD == 0) %>% # Se seleccionan 15 datos correspondientes a los pacientes
    normtest_batery(data=., vector = c(16:27,31,34), alpha = 0.05) %>% 
    as.data.frame() %>% View()

#  ##########################################################################################-
# Estadística descriptiva -------------------------------------------------
#  ##########################################################################################-  
  # Se realiza un análisis desriptivo de las covariable de los pacientes y se obtiene de las 
  # variables continuas media, desviación estándar, límite inferior (IC95%), límite superior 
  # (IC95%), mínimo, máximo, Q1, mediana, Q3, y rango intercuartílico (IQR).
  
    descriptiva = function(data, vector){
    df = matrix(nrow = length(unique(vector)), ncol = 10)
    for (j in vector) {
      X = dplyr::pull(data,j) # Selecciona como un vector atómico a una columna 
      i = match(j,vector) # Encuentra la posición en el vectror
      df[i,1] = mean(X) %>% round(.,3)
      df[i,2] = sd(X) %>% round(.,3)
      df[i,3] = mean(X)-((sd(X)*qt(1-0.1/2, length(X)-1))/sqrt(length(X))) %>% round(.,3)
      df[i,4] = mean(X)+((sd(X)*qt(1-0.1/2, length(X)-1))/sqrt(length(X))) %>% round(.,3)
      df[i,5] = fivenum(X)[1] %>% round(.,3)
      df[i,6] = fivenum(X)[2] %>% round(.,3)
      df[i,7] = fivenum(X)[3] %>% round(.,3)
      df[i,8] = fivenum(X)[4] %>% round(.,3)
      df[i,9] = fivenum(X)[5] %>% round(.,3)
      df[i,10] = IQR(X) %>% round(.,3)
    }
    colnames(df) = c('Media', 'SD', 'LI', 'LS','MIN','Q1','Mediana','Q3','MAX','IQR')
    df = as.data.frame(df)
    
    for (j in vector) {
      i = match(j,vector) # Encuentra la posición en el vectror
      row.names(df)[i] = colnames(data[,j])
    }
    return(df)
    }
  
#  ##########################################################################################-
  data %>% 
    filter(TAD == 0) %>% # Se seleccionan 15 datos correspondientes a los pacientes
    descriptiva(data=., vector = c(16:27,31,34)) %>% 
    # select(-c('SD','Q1','Mediana','Q3','IQR')) %>% 
    as.data.frame() %>% View()
# Se realizaron análisis descriptivos de los datos.
#  ##########################################################################################-  

#  ##########################################################################################-  
# Estudio de correlaciones entre variables del paciente -------------------
#  ##########################################################################################-  
# Se realizó un análisis de correlación entre variables de los pacientes por un método gráfico.
  pdf(file = './RESULTADOS/Correlacion_Continua_1.pdf', width = 11.0,height = 8.50);{
  chart.Correlation(data[data[,'TAD']==0,
                         c(16,17,18,31,19,20,21,23,24,27)], 
                    histogram=TRUE, pch=19)}; dev.off()
  
  chart.Correlation(data[data[,'TAD']==0,
                         c(20:21,22:23,24:25,26:27)], 
                    histogram=TRUE, pch=19)

# Se discute en el texto las implicaciones del análisis de correlación.
#####################################################################################################-  
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

###########################################################################-
# Resumen de valores dentro del intervalo definido
###########################################################################-  
  data %>% 
    .[.[,'TAD']==0,
      c(16,17,18,31,19,20,21,23,24,27)] %>% 
  mutate_all(funs(out_det(vec=., val=.)))
  
##########################################################################-
# Determinación de coeficientes de correlación de Pearson -----------------
##########################################################################-
# Cálculo de coeficiente de correlación r2
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  Función modificada con Indices
##  1 Seleccionar los datos por medio de índices
##  2 Se calcula la correlación entre dos vectores seleccionados desde el
##  data.frame, por sus nombres x y y.
##  3 Retorna esta correlación
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
corr <- function(data, indices, x, y) {
  df <- data[indices,]
  Z = cor(df[,x], df[,y])
  return(Z)
}

Confint_Boot <- function(data, x, y) {
    B <- boot::boot(
        data = data,
        statistic = corr,
        R = 1000,
        x = x,
        y = y)
    C <- boot::boot.ci(B, type = "basic")
    return(list(B, C))
  }

  
# Confint_Boot(data, "AGEA", "HCM")
# Confint_Boot(data, "AGEA", "SCM2")
# Confint_Boot(data, "AGEA", "CLCRMLMIN")
# Confint_Boot(data, "AGEA", "ALBGDL")
# Confint_Boot(data, "WTKG", "SCM2")
# Confint_Boot(data, "WTKG", "HCM")
# Confint_Boot(data, "WTKG", "IMCKGM2")
# Confint_Boot(data, "WTKG", "ALBGDL")
# 
# Confint_Boot(data, "ALBGDL", "SCM2")
# Confint_Boot(data, "ALBGDL", "HCM")
# Confint_Boot(data, "ALBGDL", "SCRMGDL")
# 
# Confint_Boot(data, "HCM", "SCM2")


#  ##########################################################################################-   
#  ##########################################################################################-    
# Se crea una función que realiza el test de Kendall entre un par de variables al set de datos  
  kendall = function(xcol, ycol){
    xcol_qu = rlang::enquo(xcol); ycol_qu = rlang::enquo(ycol)
      A= data %>% 
        .[.[,'TAD']==0,] %>% 
        select(.,!!xcol_qu,!!ycol_qu) %>% 
        mutate(!!xcol_qu := as.numeric(!!xcol_qu)) %>%
        mutate(!!ycol_qu := as.numeric(!!ycol_qu)) %>% 
        as.data.frame(.)
      B = cor.test(x=A[,1], y=A[,2], method = 'kendall', 
                   alternative = 'greater')
      return(list(B$p.value, B$estimate[[1]]))
  }
  
  
  kendall(SEXF, AGEA)[[1]]; kendall(SEXF, WTKG)[[1]]
  kendall(SEXF, HCM)[[1]];   kendall(SEXF, IMCKGM2)[[1]]
  kendall(SEXF, SCRMGDL)[[1]];   kendall(SEXF, ALBGDL)[[1]]
  kendall(SEXF, SCM2)[[1]];   kendall(LMP, AGEA)[[1]]
  kendall(LMP, WTKG)[[1]];   kendall(LMP, HCM)[[1]]
  kendall(LMP, IMCKGM2)[[1]];   kendall(LMP, SCM2)[[1]]
  kendall(LLP, SCM2)[[1]]
  
  plot1 <- function(df, xcol, ycol){
    xcol_qu = rlang::enquo(xcol); ycol_qu = rlang::enquo(ycol)
    dplyr::filter(df, EVID==1 & TAD==0) %>% 
      ggplot(aes(x = !!xcol_qu, y = !!ycol_qu)) +
      theme_void() +
      theme(plot.title = element_text(size=14, hjust = 0.5), 
            panel.grid = element_blank(), 
            axis.line = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            plot.background = element_rect(colour = 'black'),
            plot.margin = margin(0,0,0,0,unit='cm'),
            axis.text = element_blank()) +
      geom_smooth(method = 'loess', formula = y~x, se = F, lty='solid',color='red2',size=0.4) +
      geom_smooth(method = 'lm', formula = y~x, se = F, lty='solid', color='blue')+
      geom_point(shape=1,size=1)
  }
  plot2 <- function(df, xcol, ycol, h = 0) {
    xcol_qu = rlang::enquo(xcol); ycol_qu = rlang::enquo(ycol)
    dplyr::filter(df, EVID==1&TAD==0) %>% 
      ggplot(aes(x = !!xcol_qu, y = !!ycol_qu, col = !!xcol_qu)) +
      theme_bw() +
      theme(legend.position = "none",
            panel.grid = element_blank()) +
      geom_boxplot(shape=1, size=0.4) +
      geom_point() +
      scale_color_hue(h.start = h)
      }

  plot2(data,SEXF,WTKG, 20)
  
################################################################################################-
# Gráficos de relaciones de covariables significativas
################################################################################################-
  
  dir_res = c("./RESULTADOS/kendal_")
  
  pdf(paste0(dir_res,"SEXF_1.pdf"), width=3.5, height=3.0);{
  plot2(data,SEXF,WTKG,100)}; dev.off()
  pdf(paste0(dir_res,"SEXF_2.pdf"), width=3.5, height=3.0);{
    plot2(data,SEXF,HCM,100)}; dev.off()
  pdf(paste0(dir_res,"SEXF_3.pdf"), width=3.5, height=3.0);{
    plot2(data,SEXF,SCRMGDL,100)}; dev.off()
  pdf(paste0(dir_res,"SEXF_4.pdf"), width=3.5, height=3.0);{
    plot2(data,SEXF,ALBGDL,100)}; dev.off()
  pdf(paste0(dir_res,"SEXF_5.pdf"), width=3.5, height=3.0);{
    plot2(data,SEXF,SCM2,100)}; dev.off()
  pdf(paste0(dir_res,"LLP_1.pdf"), width=3.5, height=3.0);{
    plot2(data,LLP,ALBGDL,20)}; dev.off()
  pdf(paste0(dir_res,"LLP_2.pdf"), width=3.5, height=3.0);{
    plot2(data,LLP,CLCRMLMIN,20)}; dev.off()
  pdf(paste0(dir_res,"LLP_3.pdf"), width=3.5, height=3.0);{
    plot2(data,LLP,SCM2,20)}; dev.off()
  pdf(paste0(dir_res,"LMP_1.pdf"), width=3.5, height=3.0);{
    plot2(data,LMP,HCM,20)}; dev.off()

###################################################################################################-
###################################################################################################-
# REVISIÓN DE DIFERENCIAS ENTRE GRUPOS POR VARIABLE DICOTÓMICA  
###################################################################################################-
###################################################################################################-
  
  comp_func_1 = function(varg, var1, arg){
    varg_q = rlang::enquo(varg); var1_q = rlang::enquo(var1)
    
    if (arg == 1) {
      A = data %>% 
        filter(TAD == 0) %>% 
        group_by(!!varg_q) %>% 
        summarise(mn = mean(!!var1_q), 
                  Q1 = fivenum(!!var1_q)[2],
                  Q3 = fivenum(!!var1_q)[4])
      return(A)  
    } else {
      B = data %>% 
        filter(TAD == 0) %>% 
        group_by(!!varg_q) %>% 
        summarise(mn = mean(!!var1_q), 
                  iqr = IQR(!!var1_q))
      return(B)
    }
  }
    
  
  comp_func_1(varg = SEXF, var1 = WTKG, 2)
  comp_func_1(varg = SEXF, var1 = HCM, 2)
  comp_func_1(varg = SEXF, var1 = SCRMGDL, 2)
  comp_func_1(varg = SEXF, var1 = ALBGDL, 2)
  
  
  
  comp_func_1(varg = LLP, var1 = ALBGDL,2)
  comp_func_1(varg = LLP, var1 = CLCRMLMIN,2)
  comp_func_1(varg = LLP, var1 = SCM2,2)
  
  comp_func_1(varg = LMP, var1 = HCM,2)
  
  
  
  # vplayout <- function(x, y) viewport(layout.pos.row = y, layout.pos.col = x)
  # grid.newpage();{
  # pushViewport(viewport(layout = grid.layout(16, 16)))
  # 
  plotly = function(x,y,X,Y,z){
    X_qu = rlang::enquo(X); Y_qu = rlang::enquo(Y)
    if (z==1) {print(plot1(data, !!X_qu, !!Y_qu), vp=vplayout(x, y))}
    if (z==2) {print(plot2(data, !!X_qu, !!Y_qu), vp=vplayout(x, y))}
  }
  
  plotly(6,1,SEXF,AGEA,2)
  # 
  # plotly(2,1,WTKG,AGEA,1);   plotly(3,1,HCM, AGEA,1)
  # plotly(4,1,IMCKGM2, AGEA,1);   plotly(5,1,SCM2, AGEA,1)
  # plotly(6,1,SEXF,AGEA,2);   plotly(7,1,SCRMGDL,AGEA,1)
  # plotly(8,1,CLCRMLMIN,AGEA,1);   plotly(9,1,PROGDL,AGEA,1)
  # plotly(10,1,ALBGDL,AGEA,1);   plotly(11,1,DIND,AGEA,1)
  # plotly(12,1,DNFD,AGEA,1);   plotly(13,1,RAL,AGEA,1)
  # plotly(14,1,RAN,AGEA,1);   plotly(15,1,ANTU,AGEA,2)
  # plotly(16,1,Dx,AGEA,2)
  # 
  # plotly(3,2,HCM, WTKG,1);   plotly(4,2,IMCKGM2, WTKG,1)
  # plotly(5,2,SCM2, WTKG,1);   plotly(6,2,SEXF,WTKG,2)
  # plotly(7,2,SCRMGDL,WTKG,1);  plotly(8,2,CLCRMLMIN,WTKG,1)
  # plotly(9,2,PROGDL,WTKG,1);  plotly(10,2,ALBGDL,WTKG,1)
  # plotly(11,2,DIND,WTKG,1);  plotly(12,2,DNFD,WTKG,1)
  # plotly(13,2,RAL,WTKG,1);  plotly(14,2,RAN,WTKG,1)
  # plotly(15,2,ANTU,WTKG,2);  plotly(16,2,Dx,WTKG,2)  
  # 
  # plotly(4,3,IMCKGM2, HCM,1);   plotly(5,3,SCM2, HCM,1)
  # plotly(6,3,SEXF,HCM,2);  plotly(7,3,SCRMGDL,HCM,1)
  # plotly(8,3,CLCRMLMIN,HCM,1);  plotly(9,3,PROGDL,HCM,1)
  # plotly(10,3,ALBGDL,HCM,1);  plotly(11,3,DIND,HCM,1)
  # plotly(12,3,DNFD,HCM,1);  plotly(13,3,RAL,HCM,1)
  # plotly(14,3,RAN,HCM,1);  plotly(15,3,ANTU,HCM,2)  
  # plotly(16,3,Dx,HCM,2)  ;plotly(5,4,SCM2, IMCKGM2,1)
  # plotly(6,4,SEXF,IMCKGM2,2); plotly(7,4,SCRMGDL,IMCKGM2,1)
  # plotly(8,4,CLCRMLMIN,IMCKGM2,1);   plotly(9,4,PROGDL,IMCKGM2,1)
  # plotly(10,4,ALBGDL,IMCKGM2,1);  plotly(11,4,DIND,IMCKGM2,1)
  # plotly(12,4,DNFD,IMCKGM2,1);  plotly(13,4,RAL,IMCKGM2,1)
  # plotly(14,4,RAN,IMCKGM2,1);  plotly(15,4,ANTU,IMCKGM2,2)
  # plotly(16,4,Dx,IMCKGM2,2)  
  # 
  # plotly(6,5,SEXF,SCM2,2);   plotly(7,5,SCRMGDL,SCM2,1)
  # plotly(8,5,CLCRMLMIN,SCM2,1);  plotly(9,5,PROGDL,SCM2,1)
  # plotly(10,5,ALBGDL,SCM2,1);  plotly(11,5,DIND,SCM2,1)
  # plotly(12,5,DNFD,SCM2,1);  plotly(13,5,RAL,SCM2,1)
  # plotly(14,5,RAN,SCM2,1);  plotly(15,5,ANTU,SCM2,2)
  # plotly(16,5,Dx,SCM2,2)
  # 
  # plotly(7,6,SCRMGDL,SEXF,2);  plotly(8,6,CLCRMLMIN,SEXF,2)
  # plotly(9,6,PROGDL,SEXF,2);  plotly(10,6,ALBGDL,SEXF,2)
  # plotly(11,6,DIND,SEXF,2);  plotly(12,6,DNFD,SEXF,2)
  # plotly(13,6,RAL,SEXF,2);  plotly(14,6,RAN,SEXF,2)
  # plotly(15,6,ANTU,SEXF,2);  plotly(16,6,Dx,SEXF,2) 
  # 
  # plotly(8,7,CLCRMLMIN,SCRMGDL,1);  plotly(9,7,PROGDL,SCRMGDL,1)
  # plotly(10,7,ALBGDL,SCRMGDL,1);  plotly(11,7,DIND,SCRMGDL,1)
  # plotly(12,7,DNFD,SCRMGDL,1);  plotly(13,7,RAL,SCRMGDL,1)
  # plotly(14,7,RAN,SCRMGDL,1);  plotly(15,7,ANTU,SCRMGDL,2)
  # plotly(16,7,Dx,SCRMGDL,2) 
  # }
  
  
###################################################################################################-
###################################################################################################-
# PRIMERA REVISIÓN DE DATOS FARMACOCINÉTICOS
###################################################################################################-
###################################################################################################-

theme_set(theme_bw() + 
            theme(legend.position = c(0.8, 0.85), 
                  legend.background = element_rect(fill = NULL, colour = 'gray50'), 
                  legend.key = element_rect(fill = NULL), 
                  legend.title = element_text(face = 'bold'), 
                  legend.spacing.y = unit(0.01, 'cm'),
                  legend.key.height = unit(0.30, 'cm')))
            
##########################################################################-
# Gráfico de concentración plasmática en escala normal 
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1. Data
##  2. Convertir la variable *ID* en un factor discreto
##  3. Abrir el ambien GGplot
##  4. Adicionar puntos, y adicionar líneas
##  5. Configurar etiqueta eje X
##  6. Configurar etiqueta eje Y
##  7. Configurar etiquetas para mostrar 5 filas
##  8. Adicionar línea en LLOQ
##  9. Configurar colores en escala discreta como colores Viridis
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

gperfil1 <- data %>% 
  mutate(ID = factor(ID)) %>% 
  filter(EVID == 0) %>% 
  ggplot(., aes(x = TAD, y = DV, group = ID, col = ID)) + 
  geom_point() + geom_line() + 
  xlab('TAD: tiempo tras dosis (hr)') +
  ylab('Concentración plasmática (mg/L)') + 
  guides(col = guide_legend(nrow=5)) + 
  geom_hline(yintercept = LOQ, lty = 'dashed') +
  scale_color_viridis_d()

##########################################################################-
# Resultados de 
ggsave(filename = './RESULTADOS/TAD_DV.pdf', plot = gperfil1, 
       device = 'pdf', width = 6, height = 4, units = 'in')

##########################################################################-
# Gráficos de perfiles plasmáticos

data %>% 
  filter(EVID == 0) %>% 
  mutate(DATEHOUR = date(DATEHOUR)) %>% 
  ggplot(., aes(x=DATEHOUR, y=DV, group=ID, col=ID)) + 
  geom_point() + 
  geom_line()  +
  xlab('Fecha') +
  ylab('Concentración plasmática (mg/L)') + 
  scale_x_date(date_breaks = '1 month',
               date_minor_breaks = '1 weeks',
               date_labels = "%b-%y") +
  guides(col = guide_legend(nrow=5)) +
  theme(panel.grid.minor.x = element_line(colour='gray96'))

data %>% 
  filter(EVID ==0) %>% 
  mutate(LOQ = if_else(condition = (DV <= 3.15), 
                       true = TRUE, 
                       false = FALSE)) %>% 
  count(LOQ)

##########################################################################-
# Gráfico de concentración plasmática en escala logarítmica
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1. Data
##  2. Convertir la variable *ID* en un factor discreto
##  3. Abrir el ambien GGplot
##  4. Adicionar puntos, y adicionar líneas
##  5. Configurar etiqueta eje X
##  6. Configurar etiqueta eje Y
##  7. Configurar etiquetas para mostrar 5 filas
##  8. Adicionar línea en LLOQ
##  9. Convertir la escala de Y en logarítmo
##  10. Configurar colores en escala discreta como colores Viridis
##  11. Adicionar puntos de quiebre en la escala
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

gperfil2 <- 
  data %>% 
  mutate(ID = factor(ID)) %>% 
  filter(EVID == 0) %>% 
  ggplot(., aes(x = TAD, y = DV, group = ID, col = ID)) + 
  geom_point() + geom_line() + 
  xlab('TAD: tiempo tras dosis (hr)') +
  ylab(expression(log(C[p]))) + 
  guides(col = guide_legend(nrow=5)) + 
  geom_hline(yintercept = LOQ, lty = 'dashed') +
  scale_y_continuous(trans = 'log10') +
  scale_color_viridis_d() +
  annotation_logticks(sides = 'lr')

##########################################################################-
# Resultados de 
ggsave(filename = './RESULTADOS/log_TAD_DV.pdf', plot = gperfil2, 
       device = 'pdf', width = 6, height = 4, units = 'in')


###############################################################################################-
###############################################################################################-
#### REVISIÓN DE ESTACIONALIDAD DE LOS DATOS
###############################################################################################-
###############################################################################################-

data1 = data %>% 
  filter(EVID ==0) %>% 
  group_by(ID) %>% 
  summarise(valle = min(DV)) %>% 
  left_join(.,data,by=c("ID","valle"="DV")) %>% 
  rename("DV" = "valle") %>% 
  select(ID,DV,TSFD,TAD)


data %>% 
  filter(EVID == 0) %>% 
  ggplot(., aes(x=TSFD, y=DV, group=ID, col=ID)) + 
  geom_point() + 
  geom_line() +
  xlab('TSFD: tiempo tras la primera dosis (hr)') +
  ylab('Concentración plasmática (mg/L)') + 
  geom_smooth(method = 'loess', data = data1, aes(y=DV, x = TSFD),
              inherit.aes = F,
              linetype='dashed') +
  # guides(col = guide_legend(nrow=5)) +
  coord_cartesian(ylim = c(0,100))+
  scale_color_viridis_c()
