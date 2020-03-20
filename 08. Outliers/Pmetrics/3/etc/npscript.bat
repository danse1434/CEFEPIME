cd "F:\Documentos\(Proyecto)_Estudio_PKPD\CEFEPIME\08. Outliers\Pmetrics\J12\3"
echo Windows>time.txt
echo %time%>>time.txt
np_prep DOS < PMcontrol
echo 1 > extnum
echo go > go
gfortran -O3 -w -fopenmp -fmax-stack-var-size=32768 -o np_run "C:\Users\ACER\AppData\Roaming\Pmetrics\compiledFortran\pNPeng.o" npagdriv.f
np_run <go
echo. & echo Cleaning up.... & echo.
echo off
mkdir inputs
mkdir outputs
mkdir wrkcopy
mkdir etc
echo LASTda1.csv >> NP_RF0001.TXT
if not exist NP_RF0001.TXT (set error=1) ELSE (set error=0)
if exist DEN* move DEN* outputs
if exist OUT0* move OUT0* outputs
if exist OUTT* move OUTT* outputs
if exist PRTB* move PRTB* outputs
if exist ILOG* move ILOG* outputs
if exist NP_RF* move NP_RF* outputs
if exist ERROR* move ERROR* outputs
move instr.inx etc
move PMcontrol etc
move modini.for etc\modini.for
move modini.txt inputs\modini.txt
move XQZPJ*.ZMQ wrkcopy
move extnum etc
move npag*.* etc
erase CHMAX*.*
if exist FROM0001 move FROM0001 inputs
erase fort.*
erase go
move np_prep* etc
move np_run* etc
move LASTda1.csv inputs
echo %time% >> time.txt
move time.txt outputs
if %error% == 0 (
"C:\Program Files\R\R-3.6.0\bin\x64\Rscript" "C:/Program Files/R/R-3.6.0/library/Pmetrics/report/NPrepScript.R" "F:\Documentos\(Proyecto)_Estudio_PKPD\CEFEPIME\08. Outliers\Pmetrics\J12\3\outputs" median TRUE
start "NPAG Report" "F:\Documentos\(Proyecto)_Estudio_PKPD\CEFEPIME\08. Outliers\Pmetrics\J12\3\outputs\NPAGreport.html")
copy npscript.bat etc
echo. & echo Press any key to complete run and close this window... & echo.
pause > nul
erase npscript.bat 
