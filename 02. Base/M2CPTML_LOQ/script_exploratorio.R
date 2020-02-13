require(mlxR)
require(tidyverse)

model1 = inlineModel("
[LONGITUDINAL]
input = {Cl, V1, Q, V2, a, b}

EQUATION: ;Parameter transformations 
V = V1 
k = Cl/V1 
k12 = Q/V1 
k21 = Q/V2

; PK model definition
Cc = pkmodel(V, k, k12, k21)
g = a+b*Cc

DEFINITION:
y = {distribution = normal, reference = Cc, sd = g}

OUTPUT:
output = {Cc, y}

[INDIVIDUAL]
input = {V1_pop, Cl_pop, V2_pop, Q_pop,
        omega_V1, omega_Cl, omega_V2, omega_Q}

DEFINITION: 
V1 = {distribution = lognormal, reference = V1_pop, sd = omega_V1}
V2 = {distribution = lognormal, reference = V2_pop, sd = omega_V2}
Cl = {distribution = lognormal, reference = Cl_pop, sd = omega_Cl}
Q = {distribution = lognormal, reference = Q_pop, sd = omega_Q}
")


pPK = c(V1_pop = 29.38942, V2_pop = 22.74502, 
        Cl_pop = 13.61520, Q_pop = 3.257956,
        omega_V1 = sqrt(log(((45.36368/100)^2)+1)),
        omega_V2 = sqrt(log(((0.9445930/100)^2)+1)),
        omega_Cl = sqrt(log(((12.72771/100)^2)+1)),
        omega_Q = sqrt(log(((37.57154/100)^2)+1)),
        a = 0.695, b = 0.225)

y = list(name='y',time=seq(from=0, to=24*8.5, by=0.1)) # U: mg/L  
tto = list(time=seq(0,24*8,8), amount = 2000, tinf = 0.5)
g1 = list(size=1000, level='individual')

res1 = simulx(model = model1, 
              parameter = pPK,
              output = y, 
              group = g1,
              treatment = tto,
              settings = list(seed=102010))



DATA <- read_csv("./(Proyecto)_Estudio_PKPD/CEFEPIME/02. Base/NONMEM/DATA/NONMEM_Cefepime.csv", 
                 col_types = cols(DV = col_double(),
                                  ID = col_factor()), 
                 locale = locale())

DATA1 = DATA %>% 
  filter(EVID == 0)

G1 = prctilemlx(res1$y) + 
  geom_point(data = DATA1, aes(x=TAD+72, y=DV, group=ID, colour=ID), inherit.aes = FALSE) +
  theme_bw() +
  ylab('Concentración plasmática (mg/L)') + xlab('Tiempo (horas)') +
  coord_cartesian(xlim=c(72,72+9))+
  theme(legend.position = 'none')
  
G1 %>% ggsave(filename = 'Perfil_2cptm.pdf', device = 'pdf', width = 10, height = 4, units = 'in')  

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
  geom_histogram(binwidth = a, fill='green1', colour='green4')+
  geom_density(aes(y=a * ..count..))+
  geom_vline(xintercept = 25/2/2/2/2, lty='dashed') +
  xlab(expression(C[min]~(SS))) + ylab("Conteo")+
  theme_bw()

a = (range(res2$Cc$cmax)[2]-range(res2$Cc$cmax)[1])/38
ggplot(data = res2$Cc, mapping = aes(x=cmax)) +
  geom_histogram(binwidth = a, fill='red1', colour='red4')+
  geom_density(aes(y=a * ..count..))+
  geom_vline(xintercept = 25*2, lty='dashed')+
  xlab(expression(C[max]~(SS))) + ylab("Conteo")+
  theme_bw()

a = (range(res2$Cc$auc)[2]-range(res2$Cc$auc)[1])/38
ggplot(data = res2$Cc, mapping = aes(x=auc)) +
  geom_histogram(binwidth = a, fill='blue1', colour='blue4')+
  geom_density(aes(y=a * ..count..)) +
  xlab(expression(AUC[a]^b~(mg*h/L))) + ylab("Conteo")+
  theme_bw()


res2$Cc %>% 
  select(auc,cmin, cmax) %>% 
  summarise_all(list(
    Q1 = ~fivenum(.)[2],
    Q2 = ~fivenum(.)[3],
    Q3 = ~fivenum(.)[4],
    IQR = ~IQR(.)
  ))


