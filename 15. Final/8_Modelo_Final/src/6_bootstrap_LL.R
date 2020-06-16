require(Rsmlx)

project <- '1_M_Error_1/M_Error_bootstrap.mlxtran'

r <- bootmlx(project, nboot=1000, parametric = FALSE,
             tasks = c(populationParameterEstimation = TRUE) )

bootmlx(project, dataFolder='1_M_Error_1/BASE_MODEL/bootstrap/', parametric = FALSE,
        tasks = c(plot=TRUE) )


r.prl <- confintmlx(project, method="proflike", 
                    parameters = c("V1_pop", "Cl_pop", "V2_pop", 
                                   "beta_Cl_tSCRMGDL", "omega_Q"))

# /**********************************************************************/ 
# parameter  omega_Q 
# Value  0.845 
# CI =  [0.089 , 1.236]  
# diff. =  [-0.756 , 0.391] 
# rel. diff. =  [-89.406 , 46.301] 
# /**********************************************************************/ 
#   parameter  V1_pop 
# Value  23.814 
# CI =  [20.789 , 30.071]  
# diff. =  [-3.025 , 6.257] 
# rel. diff. =  [-12.702 , 26.274] 
# /**********************************************************************/ 
#   parameter  Cl_pop 
# Value  20.594 
# CI =  [14.386 , 28.689]  
# diff. =  [-6.208 , 8.095] 
# rel. diff. =  [-30.144 , 39.307] 
# /**********************************************************************/ 
#   parameter  V2_pop 
# Value  12.891 
# CI =  [6.517 , 18.418]  
# diff. =  [-6.374 , 5.526] 
# rel. diff. =  [-49.441 , 42.872] 
# /**********************************************************************/ 
#   parameter  beta_Cl_tSCRMGDL 
# Value  -0.415 
# CI =  [-0.732 , -0.081]  
# diff. =  [-0.318 , 0.333] 
# rel. diff. =  [76.72 , -80.496] 
# /**********************************************************************/ 
#   parameter  omega_Q 
# Value  0.845 
# CI =  [0.089 , 1.236]  
# diff. =  [-0.756 , 0.391] 
# rel. diff. =  [-89.406 , 46.301] 
