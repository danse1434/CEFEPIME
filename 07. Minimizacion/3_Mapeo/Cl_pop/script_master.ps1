for ($i=1; $i -le 100; $i++) 
    {
  $StartTime = $(get-date)
  "Carpeta $i iniciada"
  $Var = "F:\Documentos\(Proyecto)_Estudio_PKPD\CEFEPIME\07. Minimizacion\3_Mapeo\Clearance\assessment\A$i\MODELO_BASE.mlxtran"
  
  C:\ProgramData\Lixoft\MonolixSuite2019R1\bin\monolix.bat --no-gui --thread 8 --mode "none" -p $Var
  $elapsedTime = $(get-date) - $StartTime
  $totalTime = "{0:HH:mm:ss}" -f ([datetime]$elapsedTime.Ticks)
  $totalTime
    }
