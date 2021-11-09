ECHO OFF
for /l %%x in (1, 1, 50) do (
   echo Lote %%x
   "C:\Program Files\R\R-3.6.3\bin\x64\R.exe" -e "source('104_regimenesDosificacionRandom.R', encoding = 'UTF-8')"
)

PAUSE