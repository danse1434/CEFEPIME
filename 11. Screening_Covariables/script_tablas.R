##------------------------------------------------------------------------#
## Nombre del Script: análisis estadístico de tamizaje covariables de estudio
##  
## Proposito del Script:  
##  
## Autor: Daniel S. Parra Gonzalez 
## Fecha de creacion:  21-03-2020
##  
## Copyright (c) Daniel S. Parra, 2020 
##  
## Email: dsparrag@unal.edu.co 
##------------------------------------------------------------------------#
# Especificar directorio de trabajo
setwd(file.path('F:', 'Documentos', '(Proyecto)_Estudio_PKPD', 'CEFEPIME', 
                '11. Screening_Covariables'))

# Función de lectura de archivo en líneas específicas
source_lines <- function(file, lines) {
  source(textConnection(readLines(file)[lines]))
}


# Lectura de archivo sólo parte inicial
source_lines('script.R', 1:89)

##########################################################################-
# Introducción ------------------------------------------------------------
##########################################################################-
#' Función para la creación de fórmulas en lm() con evaluación tidy
#'
#' @param x variable dependiente (respuesta eta)
#' @param y variable independiente (covariable)
#' @param flatten condicional si se debe aplanar el quosure
#'
#' @return objeto de tipo fórmula que utiliza las variables con environment
#' @export
#' @examples f(disp, am + (!! sym(var)))
#' 
f <- function(x, y, flatten = TRUE) {
  # Creación de quosures
  x <- enquo(x); y <- enquo(y)
  
  # Environments de quosures debe ser el mismo en \code{x} y \code{y}
  # Podrían ser diferentes si son forwarded mediante puntos
  env <- get_env(x) # Environment de \code{x} 
  stopifnot(identical(env, get_env(y)))
  
  # Aplanar los quosures. Esto avisa al usuario, si los quosures son anidados. 
  # Estos no son soporados por funciones como lm().
  if (flatten) {
    x <- quo_squash(x, warn = TRUE)
    y <- quo_squash(y, warn = TRUE)
  }
  
  new_formula(x, y, env = env)
}

# Transformación de variables categóricas en factor
data_ori <- data_ori %>%
  mutate(SEXF = factor(SEXF),
         ANTU = factor(ANTU),
         LLP = factor(LLP),
         LMP = factor(LMP))

# Definición de covariables
covariates <- c('SEXF', 'AGEA', 'WTKG', 'HCM', 'IMCKGM2', 'SCM2', 'SCRMGDL', 
                'CLCRMLMIN', 'PROGDL', 'ALBGDL', 'DIND', 'RAL', 'RAN', 
                'ANTU', 'LLP', 'LMP')

##########################################################################-
#' Creación de lista con regresiones GLM
#'
#' @param eta desviación \code{eta} de interés
#' @param cov_vec vector de covariables de interés
#' @param intercepto condicional indica si la regresión tiene intercepto
#' @param df tabla de datos 
#' @return
#' @export 
#' Permite obtener una lista con regresión glm para cada covariable y 
#' desviación eta de interés
#' @examples GLM_ls_1("eta_Cl_SAEM", covariates, data_ori)
#' 
GLM_ls_1 <- function(eta, cov_vec, df, intercepto = TRUE, tipo = 'glm') {
  glm_list <- vector(mode = "list", length = length(cov_vec))
  if (tipo == 'lm') {
    for (i in 1:length(cov_vec)) {
      glm_list[[i]] <-
        lm(f(!!sym(eta), !!sym(cov_vec[i])), data = df)
    }
  } else if (tipo == 'glm') {
    if (intercepto == TRUE) {
      for (i in 1:length(cov_vec)) {
        glm_list[[i]] <-
          glm(f(!!sym(eta), !!sym(cov_vec[i])), data = df)
      }
    } else {
      for (i in 1:length(cov_vec)) {
        glm_list[[i]] <-
          glm(f(!!sym(eta), 0+!!sym(cov_vec[i])), data = df)
      }
    }
  }
  return(glm_list %>% set_names(cov_vec))
}

##########################################################################-
#' Procesamiento de lista con regresiones GLM para obtener una tabla con 
#' valores de coeficientes de cada modelo.
#'
#' @param eta tipo de desviación necesitada
#' @param inter condicional si los modelos deben tener intercepto o no
#'
#' @return
#' Esta función llama a GLM_ls_1(), y la procesa para la obtención de un 
#' data.frame. El data.frame contiene sólo los resultados de los coeficie-
#' ntes del modelo, y tiene opciones para obtener regresiones GLM() con y 
#' sin intercepto; así como regresiones lm(). 
#' @export
#' @examples
#' glm_df("eta_Cl_SAEM", TRUE, 'lm')
#' 
glm_df <- function(eta, inter, tipo = 'glm') {
  stopifnot(is_character(eta))
  stopifnot(is_bool(inter))
  
  GLM_ls_1(eta, covariates, data_ori, intercepto = inter, tipo = tipo) %>%
    map(~ summary(.x)) %>%
    map(~ magrittr::use_series(.x, "coefficients")) %>%
    map( ~ as.data.frame(.x)) %>% 
    map( ~ rownames_to_column(.x, var = 'Parametro')) %>% 
    map_dfr(~ as.data.frame(.x), .id = 'Tipo') %>%
    as_tibble() %>% 
    rename(SE = `Std. Error`, tval = `t value`, p = `Pr(>|t|)`) %>% 
    rename_at(vars('Estimate', 'SE', 'tval', 'p'), list( ~ paste0(eta, '_', .)))
}

##########################################################################-
# Tablas con resultados de regresión por eta ------------------------------
##########################################################################-
# Tabla para regresión glm() con pendiente e intercepto
DF1 <- glm_df("eta_Cl_SAEM", TRUE) %>% 
  left_join(glm_df("eta_Q_SAEM", TRUE), by = c("Tipo", "Parametro")) %>% 
  left_join(glm_df("eta_V1_SAEM", TRUE), by = c("Tipo", "Parametro")) %>% 
  left_join(glm_df("eta_V2_SAEM", TRUE), by = c("Tipo", "Parametro")) %>% 
  select(matches("Tipo|Parametro|p$")) %>% 
  gather(matches("eta"), key = 'ETA', value = 'p') %>% 
  filter(Parametro != "(Intercept)") %>% 
  mutate(Interp. = ifelse(p < 0.05, 'Signif.', 'NS'))

# Tabla para regresión glm() con pendiente 
DF2 <- glm_df("eta_Cl_SAEM", FALSE) %>% 
  left_join(glm_df("eta_Q_SAEM", FALSE), by = c("Tipo", "Parametro")) %>% 
  left_join(glm_df("eta_V1_SAEM", FALSE), by = c("Tipo", "Parametro")) %>% 
  left_join(glm_df("eta_V2_SAEM", FALSE), by = c("Tipo", "Parametro")) %>% 
  select(matches("Tipo|Parametro|p$")) %>% 
  gather(matches("eta"), key = 'ETA', value = 'p') %>% 
  filter(Parametro != "(Intercept)") %>% 
  mutate(Interp. = ifelse(p < 0.05, 'Signif.', 'NS'))

# Tabla para regresión lm() con pendiente e intercepto
DF3 <- glm_df("eta_Cl_SAEM", TRUE, 'lm') %>%
  left_join(glm_df("eta_Q_SAEM", TRUE, 'lm'), by = c("Tipo", "Parametro")) %>%
  left_join(glm_df("eta_V1_SAEM", TRUE, 'lm'), by = c("Tipo", "Parametro")) %>%
  left_join(glm_df("eta_V2_SAEM", TRUE, 'lm'), by = c("Tipo", "Parametro")) %>%
  select(matches("Tipo|Parametro|p$")) %>%
  gather(matches("eta"), key = 'ETA', value = 'p') %>%
  filter(Parametro != "(Intercept)") %>% 
  mutate(Interp. = ifelse(p < 0.05, 'Signif.', 'NS'))


##########################################################################-
# Correlación de Pearson --------------------------------------------------
##########################################################################-
#' Función de correlación adaptada
#'
#' @param data Tabla de datos 
#' @param x variable independiente (covariable)
#' @param y variable dependiente (desviación eta)
#'
#' @return coeficiente de correlación lineal de Pearson entre \code{x} y 
#' \code{y}
#' @examples
#' corr_test(data_ori, "eta_Cl_SAEM", "AGEA")
#' 
corr_test <- function(data, x, y) {
  a <- pull(data, x) %>% as.numeric(.)
  b <- pull(data, y) %>% as.numeric(.)
  cor.test(a, b, method = 'spearm', conf.level = 0.95)
}

##########################################################################-
#' Lista de correlaciones entre eta y covariables
#'
#' @param eta (carácter) correlación 
#' @param cov_vec (carácter) vector de correlación
#' @param df (data.frame) tabla con datos
#'
#' @return lista con correlación entre covariables y eta
#' @examples
#' corr_ls("eta_Q_SAEM", covariates, data_ori)
#' 
corr_ls <- function(eta, cov_vec, df) {
  c_ls <- vector(mode = "list", length = length(cov_vec))
  
    for (i in 1:length(cov_vec)) {
      c_ls[[i]] <- corr_test(df, eta, cov_vec[i])
    }
  return(c_ls %>% set_names(cov_vec))
}

##########################################################################-
#' Procesamiento de lista con correlaciones
#' Esta función llama a corr_ls() y realiza unas modificaciones para generar 
#' una tabla desde una lista.
#'
#' @param eta (carácter) correlación seleccionada
#'
#' @return tabla de datos con correlaciones entre eta seleccionado vs 
#' covariables.
#' @examples
#' corr_table("eta_Cl_SAEM")
#' 
corr_table <- function(eta) {
  corr_ls(eta, covariates, data_ori) %>%
    map(~ magrittr::use_series(.x, "p.value")) %>%
    map_dfr(~ as.data.frame(.x), .id = 'Tipo') %>%
    as_tibble() %>%
    rename(p = .x) %>%
    mutate(Interp. = ifelse(p < 0.05, 'Signif.', 'NS')) %>%
    add_column(ETA = paste0(eta, "_p"), .before = "Tipo")
}
  
##########################################################################-
# Creación de una tabla con las correlaciones para todos los etas
DF4 <- bind_rows(corr_table("eta_Cl_SAEM"),
                 corr_table("eta_V1_SAEM"),
                 corr_table("eta_Q_SAEM"),
                 corr_table("eta_V2_SAEM"))


##########################################################################-
# Modelo Aditivos Generalizados -------------------------------------------
##########################################################################-
#' Función de lista con modelos aditivos generalizados (GAM)
#'
#' @param eta (carácter) variable eta
#' @param cov_vec (carácter) vector de covariables
#' @param df (data.frame) tabla de datos
#'
#' @return lista con resultados de regresión GAM entre etas y covariables 
#' continuas
#' @examples
#' GAM_ls('eta_Cl_SAEM', cov_vec = covariates[2], df = data_ori)
#' GAM_ls('eta_Cl_SAEM', cov_vec = covariates[c(2:12)], df = data_ori)[[1]]
#' 
GAM_ls <- function(eta, cov_vec, df) {
  gam_ls <- vector(mode = "list", length = length(cov_vec))

  for (i in 1:length(cov_vec)) {
    gam_ls[[i]] <-
      mgcv::gam(f(!!sym(eta), !!expr(s(!!sym(cov_vec[i])))), 
                data = df)
  }

  return(gam_ls %>% set_names(cov_vec))
}

#' Extracción en tabla de significancia de parámetros de alisamiento GAM
#' @param eta (carácter) desviación eta de interés
#' @return tabla con significancia del parámetro de alisamiento en 
#' regresión GAM vs covariable eta
#' @examples
#' gam_df("eta_Q_SAEM")
#' 
gam_df <- function(eta) {
  stopifnot(is_character(eta))
  
  GAM_ls(eta, cov_vec = covariates[c(2:12)], df = data_ori) %>%
    map( ~ summary(.x)) %>%
    map( ~ magrittr::use_series(.x, "s.table")) %>%
    map(~ as.data.frame(.x)) %>%
    map(~ rownames_to_column(.x, var = 'Parametro')) %>%
    map_dfr( ~ as.data.frame(.x), .id = 'Tipo') %>%
    as_tibble(.) %>%
    rename(p = `p-value`) %>% 
    add_column(ETA = eta, .before = "Tipo")
}

##########################################################################-
# Función LRT para modelos GAM --------------------------------------------
##########################################################################-
#' Función conversión lista de GAM para cada modelo de covariable vs modelo 
#' base
#' @param eta (carácter) desviación eta de interés
#' @param m0 (objeto GAM) modelo sin inclusión de las covariables sólo 
#' intercepto
#'
#' @return tabla de datos con valores de AIC vs grados de libertad de modelo 
#' base (reducida) y modelo full.
#' @examples
#' m0 <- mgcv::gam(formula = eta_Cl_SAEM ~ 1, data = data_ori)
#' lrt_GAM_df('eta_Cl_SAEM', m0)
#' 
lrt_GAM_df <- function(eta, m0) {
  a <- GAM_ls(eta, covariates[c(2:12)], df = data_ori) %>%
    map(~ AIC(.x, m0)) %>%
    map(~ as.data.frame(.x)) %>%
    map(~ rownames_to_column(.x, var = 'Parametro')) %>%
    map_dfr( ~ as.data.frame(.x), .id = 'Covariable') %>%
    mutate(Parametro = case_when(Parametro == '.x' ~ "full",
                                 Parametro == 'm0' ~ "reduced"))
  
  a1 <- a %>%
    select(-df) %>%
    spread(key = 'Parametro', value = 'AIC')
  
  a2 <- a %>%
    select(-AIC) %>%
    spread(key = 'Parametro', value = 'df')
  
  left_join(a1,
            a2,
            by = c('Covariable'),
            suffix = c('_AIC', '_df')) %>%
    mutate(LRT = reduced_AIC - full_AIC + 2 * (full_df - reduced_df)) %>%
    mutate(p = pchisq(LRT, 1, lower.tail = FALSE),
           interp = if_else(p <= 0.05, '< 0.05', 'NS')) %>%
    add_column(ETA = eta, .before = "Covariable")
}


##########################################################################-
# Función LRT para modelos GLM --------------------------------------------
##########################################################################-
#' Tabla de comparación de modelos glm ('lm') para cada coviarable vs modelo 
#' base sin inclusión de covariable
#'
#' @param eta (carácter) desviación eta de interés
#' @param m0 (objeto GLM) modelo base (reducido) con sólo intercepto
#' @return tabla de datos con LRT entre modelo completo (covariable + 
#' intercepto) vs modelo reducido (sólo intercepto). Sólo se acepta la reg-
#' resión 'lm' con pendiente e intercepto
#'
#' @examples
#' m0 <- glm(eta_Cl_SAEM ~ 1, data = data_ori)
#' lrt_GLM_df('eta_Cl_SAEM', m0)
#' 
lrt_GLM_df <- function(eta, m0) {
  a <- GLM_ls_1(eta, covariates, data_ori) %>%
    map( ~ AIC(.x, m0)) %>%
    map( ~ as.data.frame(.x)) %>%
    map( ~ rownames_to_column(.x, var = 'Parametro')) %>%
    map_dfr(~ as.data.frame(.x), .id = 'Covariable') %>%
    mutate(Parametro = case_when(Parametro == '.x' ~ "full",
                                 Parametro == 'm0' ~ "reduced"))
  
  a1 <- a %>%
    select(-df) %>%
    spread(key = 'Parametro', value = 'AIC')
  
  a2 <- a %>%
    select(-AIC) %>%
    spread(key = 'Parametro', value = 'df')
    
  left_join(a1, a2, by = c('Covariable'), suffix = c('_AIC', '_df')) %>%
      mutate(LRT = reduced_AIC - full_AIC + 2 * (full_df - reduced_df)) %>%
      mutate(p = pchisq(LRT, 1, lower.tail = FALSE),
             interp = if_else(p <= 0.05, '< 0.05', 'NS')) %>%
      add_column(ETA = eta, .before = "Covariable")
}

##########################################################################-
# Unión de todas las tablas -----------------------------------------------
##########################################################################-
# Tabla únicas y exportación a CSV
DF1 %>% 
  left_join(DF2, by = c('Tipo', 'Parametro', 'ETA'), suffix = c("CO", "SI")) %>% 
  left_join(DF3, by = c('Tipo', 'Parametro', 'ETA')) %>% 
  left_join(DF4, by = c('Tipo', 'ETA'), suffix = c("lm", "pearson")) %>% 
  write_csv(path = "Figuras/tabla.csv")

##########################################################################-
# Cálculo de comparativas entre modelos full y reduced --------------------
##########################################################################-
# Modelos base (BASE)
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
##  1 Lista para con datos de LRT en modelos GLM completos y reducidos
##  2 Incluir el objeto df a cada una de los elementos de las listas
##:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
lrt_glm_ls <- vector(mode = 'list', length = 4L)

# Eta Cl
m0 <- glm(eta_Cl_SAEM ~ 1, data = data_ori)
lrt_glm_ls[[1]] <- lrt_GLM_df('eta_Cl_SAEM', m0)

# Eta Q
m0 <- glm(eta_V2_SAEM ~ 1, data = data_ori)
lrt_glm_ls[[2]] <- lrt_GLM_df('eta_Q_SAEM', m0)

# Eta V1
m0 <- glm(eta_V1_SAEM ~ 1, data = data_ori)
lrt_glm_ls[[3]] <- lrt_GLM_df('eta_V1_SAEM', m0)

# Eta V2
m0 <- glm(eta_V2_SAEM ~ 1, data = data_ori)
lrt_glm_ls[[4]] <- lrt_GLM_df('eta_V2_SAEM', m0)


# Modelos GAM
lrt_gam_ls <- vector(mode = 'list', length = 4L)

# Eta Cl
m0 <- mgcv::gam(eta_Cl_SAEM ~ 1, data = data_ori)
lrt_gam_ls[[1]] <- lrt_GAM_df('eta_Cl_SAEM', m0)

# Eta Q
m0 <- mgcv::gam(eta_Q_SAEM ~ 1, data = data_ori)
lrt_gam_ls[[2]] <- lrt_GAM_df('eta_Q_SAEM', m0)

# Eta V1
m0 <- mgcv::gam(eta_V1_SAEM ~ 1, data = data_ori)
lrt_gam_ls[[3]] <- lrt_GAM_df('eta_V1_SAEM', m0)

# Eta V2
m0 <- mgcv::gam(eta_V2_SAEM ~ 1, data = data_ori)
lrt_gam_ls[[4]] <- lrt_GAM_df('eta_V2_SAEM', m0)


# Creación de tablas con resultados LRT
lrt_glm_ls %>% 
  map_dfr(~ as.data.frame(.x)) %>% 
  write.csv(file = 'Figuras/tabla_lrt_glm.csv')

lrt_gam_ls %>% 
  map_dfr(~ as.data.frame(.x)) %>% 
  write.csv(file = 'Figuras/tabla_lrt_gam.csv')



  
  
  
  
  
    
