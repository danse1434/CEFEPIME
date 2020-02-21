A = convergence_list %>%
  magrittr::extract2(1) %>%
  select(convergenceIndicator) %>% 
  pull()


slop_sign <- function(A, i, l) {
  stopifnot(is.vector(A))
  
  if (i > l) {
    
    x = 1:(l+1)
    y = A[(i - l):i]
    
    B = data.frame(x, y)

    Zsl = lm(y ~ x, B) %>%
      summary(.) %>%
      magrittr::use_series('coefficients') %>%
      magrittr::extract2(2, 1)
    
    Zin = lm(y ~ x, B) %>%
      summary(.) %>%
      magrittr::use_series('coefficients') %>%
      magrittr::extract2(2, 2)
    
    Z_pen = 
      lm(y ~ x, B)
    
    T = Zsl/Zin
    
    return(Z_pen)
  } else {
    return(NA_real_)
  }
}

slop_sign(A, 101, 100) %>% 
  summary()

##########################################################################-
# 

# qt(p = 0.0229, df = 99, lower.tail = F)
# 
# pt(q = 2.311, df = 100, lower.tail = F)
# 
# data.frame(x = 12:14) %>%
#   mutate(y = 0.0113 - 
#            round(pt(q = 2.582, df = x, lower.tail = F),4)) %>%
#   ggplot(aes(x, y)) +
#   geom_line() +
#   geom_hline(yintercept = 0)
# 
# 2*pt(q = 2.582, df = 99, lower.tail = F)


require(Rcpp)
require(RcppArmadillo)

Rcpp::sourceCpp('funcion_CPP.cpp')

x <- func4(A = A, j = 101, l = 100)
y <- func5(A = A, j = 101, l = 100)


func6(x,y)
func7(x,y)
func8(x,y)

slop_signC(A = A, j = 101, l = 100)
  
microbenchmark::microbenchmark(slop_sign(A, 101, 100),
                               slop_signC(A, 101, 100))








