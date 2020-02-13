
require(tidyverse)

# Selección de tema
theme_set(theme_classic() +
            theme(panel.border = element_rect(fill = NA, colour = 'black')))


X <- data.frame(x = seq(0, 5, 1e-5))

X <- X %>% 
  mutate(y = dchisq(x, df = 1)) 


ggplot(X, mapping = aes(x = x, y = y)) +
  geom_line() +
  stat_function(fun = dchisq, args = list(df = 1), xlim = c(0, 5),geom = 'area') +
  geom_vline(xintercept = qchisq(p = 0.05, df = 1, lower.tail = F), lty = 'dashed')


 

qchisq(p = 0.05, df = 1, lower.tail = F)

qchisq(p = 0.01, df = 1, lower.tail = F)


pchisq(q = 3.841459, df = 1, lower.tail = F)









##########################################################################-
# Comparación de modelo de 1CPTM vs 2CPTM
# Dominio original
LRT_1 = 580.121 - 534.092 + 2 * (9 - 5)
LRT_t1 = qchisq(p = 0.05, df = 9 - 5, lower.tail = F)
pval1 = pchisq(q = 54.029, df = 9 - 5, lower.tail = F)


# Dominio logarítmico
LRT_1.ln = 62.740 - 40.726 + 2 * (9 - 5)
LRT_t1.ln = qchisq(p = 0.05, df = 9 - 5, lower.tail = F)
pval1.ln = pchisq(q = LRT_1.ln, df = 9 - 5, lower.tail = F)


##########################################################################-
# Comparación de modelo de 2CPTM vs 3CPTM
# Dominio original
LRT_2 = 534.092 - 542.394 + 2 * (13 - 9)
LRT_t2 = qchisq(p = 0.05, df = 13 - 9, lower.tail = F)
pval2 = pchisq(q = LRT_2, df = 13 - 9, lower.tail = F)


# Dominio logarítmico
LRT_2.ln = 40.726 - 47.778 + 2 * (13 - 9)
LRT_t2.ln = qchisq(p = 0.05, df = 13 - 9, lower.tail = F)
pval2.ln = pchisq(q = LRT_2.ln, df = 13 - 9, lower.tail = F)
