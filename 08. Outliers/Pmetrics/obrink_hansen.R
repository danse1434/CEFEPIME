#Use PMmanual() for help
require(tidyverse)
require(Pmetrics)

#Run 1 - add your run description here
##########################################################################-
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '08. Outliers', 'Pmetrics', 'J12'))

data <- read_delim("data.csv", ";", escape_double = FALSE, na = ".",
                   trim_ws = TRUE)
data1 <- data %>% 
  mutate(DUR = AMT/RATE) %>% 
  mutate(INPUT = case_when(EVID == 1 ~ 1,
                           TRUE ~ NA_real_)) %>%
  mutate(OUTEQ = case_when(EVID == 0 ~ 1,
                           TRUE ~ NA_real_)) %>% 
  mutate(C0 = case_when(EVID == 0 ~ 1, #0.081
                        TRUE ~ NA_real_)) %>% 
  mutate(C1 = case_when(EVID == 0 ~ 0.1, #0.1152
                        TRUE ~ NA_real_)) %>% 
  mutate(C2 = case_when(EVID == 0 ~ 0,
                        TRUE ~ NA_real_)) %>% 
  mutate(C3 = case_when(EVID == 0 ~ 0,
                        TRUE ~ NA_real_)) %>% 
  rename(OUT = DV) %>% 
  rename(DOSE = AMT) %>%
  rename("#ID" = ID) %>% 
  select(`#ID`,'EVID', 'TIME', 'DUR', 'DOSE', 'ADDL', 'II', 'INPUT', 
         'OUT', 'OUTEQ', 'C0', 'C1', 'C2', 'C3')
  
# data2 <- data1 %>%
#   filter(EVID == 1) %>%
#   mutate(TIME = TIME + DUR) %>%
#   mutate(DUR = 0) %>%
#   mutate(DOSE = 0)

data3 <- data1 %>%
  # bind_rows(data2) %>%
  arrange(`#ID`, TIME) #%>%
  # mutate(`#ID` = case_when(
  #   `#ID` <= 9 ~ paste0('ID0', `#ID`),
  #   `#ID` >= 10 ~ paste0('ID', `#ID`),
  #   TRUE ~ NA_character_
  # ))

write_csv(x = data3, path = 'Base/LASTdat.csv', na = '.')

# Añadir manualmente al archivo la frase 
# POPDATA DEC_11

# # Añadir archivos
# fileName <- 'data3.csv'
# Z = readChar(fileName, file.info(fileName)$size)
# # Z = substr(Z, 4, nchar(Z))
# # 
# write_lines(x = c("POPDATA DEC_11", Z),
#             path = 'data3.csv',
#             sep = '\n')

##########################################################################-
# ITrun(model = "modini.txt", data = "data3.csv")

##########################################################################-
# Carga última -----------------------------------------------------
##########################################################################-
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '08. Outliers', 'Pmetrics', 'J12'))

list_file = list.files(file.path(getwd(), "Base"),
                       "modini\\.txt$|LASTdat\\.csv$", full.names = T)
file.copy(list_file, to = file.path(getwd()),
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

A = PMreadMatrix("LASTdat.csv") 
plot(A, pch = 1, cex = 1.2, lwd = 0.1, col = "blue1") 
plot(A, pch = 1, cex = 1.2, lwd = 0.1, col = "blue1", tad = TRUE) 

NPrun(model = "modini.txt",
      data = "LASTdat.csv",
      cycles = 1000)



##########################################################################-
PMload(run = 1)

plot(mdata.1, pred = post.1, overlay = T, pch = 1, cex = 0.5, lwd = 0.1, 
     col = "blue3", col.pred = "blue1", xlab = "Tiempo (horas)", 
     ylab = "Conc (mg/L)")


post.1 %>% 
  ggplot(aes(time, pred, col = factor(id))) +
  geom_line() +
  geom_point(data = mdata.1, mapping = aes(x = time, y = out)) +
  facet_wrap(~ factor(id)) +
  theme_bw() +
  scale_color_viridis_d()


op.1 %>% 
  ggplot(aes(y = obs, x = pred, colour = pred.type)) +
  geom_point() + 
  stat_smooth() +
  theme_bw() +
  facet_wrap(~ pred.type)



summary(post.1)
summary(op.1)


x = final.1$popPoints[,1]
y = final.1$popPoints[,2]
z = final.1$popPoints[,3]

rgl::plot3d(x = final.1$popPoints[, 1], 
            y = final.1$popPoints[, 3], 
            z = final.1$popPoints[, 5], 
            size = 4, radius = 2, 
            xlab = 'V', ylab = 'Cl', zlab = 'P', col = 'red')




##########################################################################-
# Polinomio ---------------------------------------------------------------
##########################################################################-

C <-
  read_csv(
    "F:/Documentos/(Proyecto)_Estudio_PKPD/CEFEPIME/00. Datos/Calibracion/Datos/Calibracion.csv"
  )

C1 <- C %>% 
  mutate(Conc = exp((Diametro - 11.850) / 3.828)) %>% 
  group_by(Nivel) %>% 
  summarise(mn = mean(Conc),
         sd = sd(Conc))

xEP <- C1$mn
yEP <- C1$sd

makeErrorPoly(obs = xEP, sd = yEP)


##########################################################################-
# Lectura de -----------------------------------------------------
##########################################################################-
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '08. Outliers', 'Pmetrics', 'J12'))

data <- read_delim("Monolix_data_TAD.csv", ",", escape_double = FALSE, na = ".",
                   trim_ws = TRUE)
data1 <- data %>% 
  mutate(DUR = AMT/RATE) %>% 
  mutate(INPUT = case_when(EVID == 1 ~ 1,
                           TRUE ~ NA_real_)) %>%
  mutate(OUTEQ = case_when(EVID == 0 ~ 1,
                           TRUE ~ NA_real_)) %>% 
  mutate(C0 = case_when(EVID == 0 ~ 0.046759110, #0.081
                        TRUE ~ NA_real_)) %>% 
  mutate(C1 = case_when(EVID == 0 ~ 0.082352685, #0.1152
                        TRUE ~ NA_real_)) %>% 
  mutate(C2 = case_when(EVID == 0 ~ 0.002449326,
                        TRUE ~ NA_real_)) %>% 
  mutate(C3 = case_when(EVID == 0 ~ 0,
                        TRUE ~ NA_real_)) %>% 
  rename(OUT = DV) %>% 
  rename(DOSE = AMT) %>%
  rename("#ID" = ID) %>% 
  select(`#ID`,'EVID', 'TIME', 'DUR', 'DOSE', 'ADDL', 'II', 'INPUT', 
         'OUT', 'OUTEQ', 'C0', 'C1', 'C2', 'C3')

data3 <- data1 %>%
  arrange(`#ID`, TIME)

write_csv(x = data3, path = 'Base/LASTda1.csv', na = '.')

list_file = list.files(file.path(getwd(), "Base"),
                       "modini\\.txt$|LASTda1\\.csv$", full.names = T)
file.copy(list_file, to = file.path(getwd()),
          overwrite = TRUE, recursive = FALSE, copy.mode = TRUE)

A = PMreadMatrix("LASTda1.csv") 
plot(A, pch = 1, cex = 1.2, lwd = 0.1, col = "blue1") 
plot(A, pch = 1, cex = 1.2, lwd = 0.1, col = "blue1", tad = TRUE) 


NPrun(model = "modini.txt", data = "LASTda1.csv",
      cycles = 200)












