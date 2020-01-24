require(mlxR)
require(tidyverse)

model1 = inlineModel("
DESCRIPTION:
The administration is via an infusion (requires INFUSION RATE or INFUSION DURATION column-type in the data set).
The PK model has one compartment (volume V) and a linear elimination (clearance Cl).

[LONGITUDINAL]
input = {V, Cl}

EQUATION:
; PK model definition
Cc = pkmodel(V, Cl)

OUTPUT:
output = Cc

[INDIVIDUAL]
input = {V_pop, Cl_pop, omega_V, omega_Cl}

DEFINITION: 
V = {distribution = lognormal, reference = V_pop, sd = omega_V}
Cl = {distribution = lognormal, reference = Cl_pop, sd = omega_Cl}
")


pPK = c(V_pop = 33.20, Cl_pop = 13.1,
        omega_V = sqrt(log(((9.56/100)^2)+1)),
        omega_Cl = sqrt(log(((0.362/100)^2)+1)))
  
Cc = list(name='Cc',time=seq(from=0, to=24*8.5, by=0.1)) # U: mg/L  
tto = list(time=seq(0,24*8,8), amount = 2000, tinf = 0.5)
g1 = list(size=2000, level='individual')

res1 = simulx(model = model1, 
              parameter = pPK,
              output = list(name=c("Cc","PCA"), time='steady.state'), 
              group = g1,
              treatment = tto,
              settings = list(seed=102010))



DATA <- read_csv("(Proyecto)_Estudio_PKPD/CEFEPIME/02. Base/NONMEM/DATA/NONMEM_Cefepime.csv", 
                 col_types = cols(DV = col_double(),
                                  ID = col_factor()), 
                 locale = locale())





DATA1 = DATA %>% 
  filter(EVID == 0)


G1 = prctilemlx(res1$Cc,
           band = list(number=2, level=95)) + 
  geom_point(data = DATA1, aes(x=TSFD, y=DV, group=ID, colour=ID), inherit.aes = FALSE) +
  theme_bw() +
  ylab('Concentración plasmática (mg/L)') + xlab('Tiempo (horas)') +
  theme(legend.position = 'none')
  
G1 %>% ggsave(filename = 'Rplot.pdf', device = 'pdf', width = 10, height = 4, units = 'in')  


# Datos de Exposición -----------------------------------------------------

res2 = exposure(
  model = model1,
  parameter = pPK,
  output = list(name=c("Cc"), time='steady.state'),
  group = g1,
  treatment = list(amount = 2000, tinf = 0.5, ii = 8),
  settings = list(seed = 102010)
)


a = (range(res2$Cc$cmin)[2]-range(res2$Cc$cmin)[1])/38
ggplot(data = res2$Cc, mapping = aes(x=cmin))+
  geom_histogram(binwidth = a)+
  geom_density(aes(y=a * ..count..))+
  geom_vline(xintercept = 25/2/2/2/2, lty='dashed')

a = (range(res2$Cc$cmax)[2]-range(res2$Cc$cmax)[1])/38
ggplot(data = res2$Cc, mapping = aes(x=cmax)) +
  geom_histogram(binwidth = a)+
  geom_density(aes(y=a * ..count..))+
  geom_vline(xintercept = 25*2, lty='dashed')

a = (range(res2$Cc$auc)[2]-range(res2$Cc$auc)[1])/38
ggplot(data = res2$Cc, mapping = aes(x=auc)) +
  geom_histogram(binwidth = a, fill='blue1', colour='blue4')+
  geom_density(aes(y=a * ..count..)) +
  theme_bw()



