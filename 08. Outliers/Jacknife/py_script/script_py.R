source('script.R', encoding='UTF-8')

require(reticulate)
use_python('C:/Users/ACER/AppData/Local/Programs/Python/Python38-32/python')


source_python('py_script/script_py.py')
sumarNumeros(3,4)


py_run_file('py_script/script_py.py')
py_run_string("x = sumarNumeros(4,5.3)")
py$x


pd <- import("pandas", convert = FALSE)
np <- import("numpy", convert = FALSE)
plt <- import("matplotlib.pyplot", convert = FALSE)




a <- np$array(c(1:4))
