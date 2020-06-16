# PROYECTO FARMACOCINÉTICA POBLACIONAL DE CEFEPIME EN NF
# Simulación de Concentraciones
library(ggplot2); theme_set(theme_classic())
require(ggExtra)
require(grid); vplayout <- function(x, y) viewport(layout.pos.row = x, layout.pos.col = y)
require(magrittr)
library(plyr)
library(reshape2)

setwd(file.path("F:","Documentos","Personal","Jorge",
                "simulaciones","Paciente Típico"))


# PK-PD -------------------------------------------------------------------
# Leer archivo tipo CSV
data <- read.csv("bolusPOP.csv",header = T)

data1 <- data[,-c(2:9)]
data2 <- melt(data1, id.vars = 'Indiv..')
data2['Tiempo'] <- rep(seq(0, 23.9, 0.1), each = 10000)

# df1 <- subset(data2, data2$Indiv.. == c(1:10))
# 
# ggplot(data = df1, aes(x = Tiempo, y = value, group = Indiv..)) +
#   geom_line(alpha = 0.05) + 
#   labs(x = 'Tiempo (min)', y = 'Concentración (mg/L)') +
#   stat_summary(fun.y='mean', colour="green4", geom="line", inherit.aes = F, aes(x = Tiempo, y = value)) +
#   theme_classic()+
#   theme(panel.background = element_rect(colour = 'black'))

G1 <- ggplot(data = data2, aes(x = Tiempo, y = value, group = Indiv..)) +
  geom_line(alpha = 0.05) + 
  labs(x = 'Tiempo (min)', y = 'Concentración (mg/L)')+
  stat_summary(fun.y='mean', colour="green4", geom="line", inherit.aes = F, aes(x = Tiempo, y = value)) +
  theme_classic()+
  theme(panel.background = element_rect(colour = 'black'))
  
# ggsave(filename = 'Infusion_5_min.png', device = 'png', plot = G1, 
#        width = 5.83, height = 4.13, units = 'in',dpi = 320)



# Evaluación de Objetivo PD -----------------------------------------------
#MIC = c(0.02,0.04,0.08,0.016,0.032,0.064,0.125,0.25,
#        0.50,1,2,4,8,16,32,64,128,256,512)

MIC = c(1*(2^(seq(-10,10,0.1))))

PTA = vector(length = length(MIC))

df1 <- as.matrix(data1[,-1])

for (j in 1: length(MIC)) {
  y = vector(length=dim(df1)[1])
  y1 = for(i in 1:(dim(df1)[1])) {
    y[i] = (length(subset(df1[i,], subset = df1[i,]>=MIC[j]))/length(df1[i,]))
  }
  PTA[j] = length(y[y >= 0.5])/length(y) }

############PUNTOS DE GRAFICO##########################
MIC1 = c(1*(2^(seq(-10,10,1))))
PTA1 = vector(length = length(MIC1))

for (j in 1: length(MIC1)) {
  y = vector(length=dim(df1)[1])
  y1 = for(i in 1:(dim(df1)[1])) {
    y[i] = (length(subset(df1[i,], subset = df1[i,]>=MIC1[j]))/length(df1[i,]))
  }
  PTA1[j] = length(y[y >= 0.5])/length(y) }

# PTA vs MIC
PTA.1 = data.frame(PTA, MIC)
PTA.2 = data.frame(PTA1, MIC1)

ggplot(PTA.1, aes(x = MIC, y = PTA)) + 
  geom_line() +
  geom_point(data = PTA.2, aes(x = MIC1, y = PTA1), inherit.aes = FALSE) +
  scale_x_continuous(trans='log2',breaks = c(0.002,0.004,0.008,0.016,0.032,0.064,0.125,0.25,
                                             0.50,1,2,4,8,16,32,64,128,256,512)) +
  scale_y_continuous(breaks = c(seq(0,1, by = 0.1)))+ 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) +
  labs(title="Modelo PK-PD 100% f T > MIC", subtitle="Cefepime 2000mg Bolos - 10,000 individuos",
       x="Concentración Mínima Inhibitoria (MIC)", y="Probabilidad de Cumplimiento Objetivo (PTA)")

write.csv(x = PTA.1, file = 'PTA_Inf_Ext.csv')






  
  geom_line(data=PTA.2, aes(x=MIC1, y=PTA_1), inherit.aes = F) + 
  geom_point() + geom_line()+
  scale_x_continuous(trans='log2',breaks = c(0.02,0.04,0.08,0.016,0.032,0.064,0.125,0.25,
                                             0.50,1,2,4,8,16,32,64,128,256,512)) +
  scale_y_continuous(breaks = c(seq(0,100, by = 10))) + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1)) + 
  geom_line(aes(y = c(rep(90,19))), col = "red4") + 
  labs(title="Modelo PK-PD 100% f T > MIC", subtitle="Cefepime 2000mg Bolos - 10,000 individuos",
       x="Concentración Mínima Inhibitoria (MIC)", y="Probabilidad de Cumplimiento Objetivo (PTA)")

# Gráfico Spaguetti ------------------------------------------------------
graph1=ggplot(C_r, aes(x = V4, y = Conc.))+
  geom_line(size=0.01) + 
  labs(title="Farmacocinética - Set de datos Öbrink-Hansen K. AAC - 2015", subtitle="Modelo Poblacional de Piperacilina en Bolos",
       x="Tiempo (horas)", y="Concentración Plasmática (mg/L)") +
  scale_x_continuous(sec.axis=dup_axis(name=NULL,labels = NULL))+ 
  scale_y_continuous(sec.axis=dup_axis(name=NULL,labels = NULL)) + theme(legend.position="none")

graph2= graph1 + aes(colour = factor(Sujeto))
graph2

{graph3=ggplot(C_r, aes(x = V4, y = V5))+
    geom_point(size=0.01) + 
    labs(title="Farmacocinética - Set de datos Öbrink-Hansen K. AAC - 2015", subtitle="Modelo Poblacional PIPC Bolos - error proporcional",
         x="Tiempo (horas)", y="Concentración Plasmática (mg/L)") +
    scale_x_continuous(sec.axis=dup_axis(name=NULL,labels = NULL))+ 
    scale_y_continuous(sec.axis=dup_axis(name=NULL,labels = NULL))
  
  graph4= graph3 + aes(colour = factor(Sujeto))
  graph4}



# Superposición de Gráficos -----------------------------------------------
VPC2 <- read.csv(file = 'Datos_VPC2.csv')
data_OBS <- read.csv(file = 'Datos_Observados.csv')

vpc_colours <- c('P99' <- 'green1', 'P95' <- 'green1', 'P90' <- 'green1')
vpc_alphas <- c('P99' <- 1, 'P95' <- 0.4, 'P90' <- 0.1)

VPC_final <- ggplot(data_OBS, x = t, y = C) +
  
  geom_ribbon(data = VPC2, aes(x = VPC2$t, ymax = VPC2$P99, ymin = VPC2$P01, fill = 'P99', alpha = 'P99'), 
              inherit.aes = F) +
  geom_ribbon(data = VPC2, aes(x = VPC2$t, ymax = VPC2$P95, ymin = VPC2$P05, fill = 'P95', alpha = 'P95'), 
              inherit.aes = F) +
  geom_ribbon(data = VPC2, aes(x = VPC2$t, ymax = VPC2$P90, ymin = VPC2$P10, fill = 'P90', alpha = 'P90'), 
              inherit.aes = F) +
  geom_line(data = data2, aes(x = Tiempo, y = value, group = Indiv..), alpha = 0.05, inherit.aes = FALSE) +
  geom_point(aes(x = t, y = C)) +
  coord_cartesian(xlim=c(0,8))+
  scale_x_continuous(sec.axis=dup_axis(name=NULL,labels = NULL))+ 
  scale_y_continuous(sec.axis=dup_axis(name=NULL,labels = NULL))+
  labs(x="Tiempo (horas)", y="Concentración Plasmática (mg/L)") +
  scale_fill_manual(values = vpc_colours, name = 'Percentiles')+
  scale_alpha_manual(values = vpc_alphas, name = 'Percentiles')


 ggsave(filename = 'Comparativa.png', device = 'png', plot = VPC_final, 
        width = 5.83, height = 4.13, units = 'in',dpi = 320)

