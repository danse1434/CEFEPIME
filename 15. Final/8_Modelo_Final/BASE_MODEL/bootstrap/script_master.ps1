$main_dir = "C:\Users\Daniel\OneDrive\Documents\(Proyecto)_Estudio_PKPD\CEFEPIME\15. Final\3_Transform_Both_Sides\tbs_folder_2\boot"

for ($i=1; $i -le 1000; $i++) 
    {
  $StartTime = $(get-date)
  
  "Carpeta $i iniciada"
  
  $Var = "$main_dir\B$i\final_model.mlxtran"

  C:\ProgramData\Lixoft\MonolixSuite2019R1\bin\monolix.bat --no-gui --thread 8 --mode "none" -p $Var
  
  $elapsedTime = $(get-date) - $StartTime
  $totalTime = "{0:HH:mm:ss}" -f ([datetime]$elapsedTime.Ticks)
  $totalTime
    }
