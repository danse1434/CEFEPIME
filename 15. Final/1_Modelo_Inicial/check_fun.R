# Carga de paquetes
# require(magrittr)
# require(glue)
require(rlang)
require(tidyverse)
require(Rcpp)
require(RcppArmadillo)

Rcpp::sourceCpp('boot/funcion_definitiva.cpp')

A <- data.frame(a = c(rnorm(30, 100, 10), 
                      rchisq(30, df = 2), 
                      rt(30, 12),
                      rgamma(30, shape = 2, rate = 2), 
                      rlnorm(30, 20, 2),
                      rweibull(30, 2, 2)),
                b = rep(1:6, each = 30)) %>% 
  as_tibble() 

B <- A %>%
  nest_by(b) %>%
  mutate(mean = map_dbl(data, ~ mean(.x))) %>% 
  mutate(B = pmap(list(param=data, t0=mean, alpha=0.05), confints)) %>% 
  unnest(B) %>% 
  pivot_longer(cols = contains("."), names_to = c('metodo','tipo'), 
               values_to = 'limit', names_sep = '\\.') %>% 
  pivot_wider(names_from = 'tipo', values_from = 'limit') %>% 
  add_column(Tipo = rep(c('Normal', 'Chisq', 't-Student', 'Gamma',
                            'Log-normal', 'Weibull' ), each = 5))



# B %>%
#   unnest(data) %>%
#   ggplot() +
#   geom_point(aes(x = metodo, y = mean), col = 'red') +
#   geom_point(aes(x = metodo, y = a), position = 'jitter', size = 1) +
#   geom_errorbar(aes(x = metodo, ymin = LI, ymax = LS)) +
#   facet_wrap(. ~ Tipo, scales = 'free_x') +
#   coord_flip()


# Prueba de funciones en el set de datos de interés
test <- read_csv('boot/parametros_df.csv') %>%
  nest_by(parameter, mean)

test %>%
  mutate(B = pmap(list(param=data, t0=mean, alpha=0.05), confints)) %>%
  unnest(B) %>%
  view()
