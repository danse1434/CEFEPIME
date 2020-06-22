@ECHO OFF
REM CREACIÓN DE VINCULO SIMBÓLICO PARA FUNCIONAMIENTO CON R

REM DIRECTORIO PADRE
SET parent=F:\Documentos\(Proyecto)_Estudio_PKPD\CEFEPIME\

REM CARPETA ORIGINAL (OBJETIVO)
SET original="%parent%15. Final\8_Modelo_Final\"

REM CARPETA NUEVA
SET target="%parent%16. Simulaciones\1_Simulacion\model\"

REM VINCULO
MKLINK /D %target% %original%
PAUSE