runMonolix <- function(name, mlxInstallDir, display=NULL, saveGraphics=NULL, scenario=NULL){
  #####################################################################################################
  # Function allowing the run a Mlxtran project
  #
  # -> name: Define the Mlxtran files to run in absolute path
  # --> If it is a list of Mlxtran files, run them
  # --> If it is a list of directory, run all the Mlxtran files in it
  # -> mlxInstallDir : Define where MonolixSuite is installed
  # -> display : Boolean to define if Monolix calculation is displayed or not (FALSE by default)
  # -> saveGraphics : Boolean to define if the graphics are saved or not (FALSE by default)
  # -> scenario : Monolix scenario could be a list of 
  #     - saem: estimate population parameters
  #     - indiv: estimate individual parameters
  #     - fim: estimate the standard errors of the estimates and Fisher information matrix
  #     - ll: estimate the log-likelihood
  #     - graphics: generate the result graphics
  #   By default, the run option is used to get the scenario from the project
  # Typically for windows mlxInstallDir <- 'C:/ProgramData/Lixoft/MonolixSuite2016R1/bin/'
  #
  # Lixoft - v1.0
  #####################################################################################################
  
  runCalculation <- TRUE
  
  # Define the operating system: the associated function is different between Windows and Linux
  if(Sys.info()[1]=="Windows"){ext<- 'bat'}else{ext <- 'sh'}
  
  # Check the existence of Monolix
  monolixPath <- file.path(mlxInstallDir,paste0('/monolix.',ext))
  if(!file.exists(monolixPath)){
    warning(paste0('Unable to find monolix.bat with your definition of mlxInstallDir: ',mlxInstallDir))
    runCalculation <- FALSE
  }
  
  # Check the existance of all the Mlxtran projects and/or directory
  if(!is.vector(name)){name<-c(name)}
  # Check if it is a list of directory or files
  files<- NULL
  for(index_name in 1:length(name)){
    if(file.exists(name[index_name])||dir.exists(name[index_name])){
      nameInfo <- file.info(name[index_name])
      if(nameInfo$isdir){
        files <- c(files,file.path(name[index_name],list.files(path = name[index_name], pattern = ".mlxtran$"))) 
      }else{
        files <- c(files,name[index_name])
      }
    }else{
      base::cat(paste0('File/folder :',name[index_name],' does not exist \n'))
    }
  }
  
  # initialization of the default values
  if(is.null(display)){display <- FALSE}
  if(display){displayCommand <- ''}else{displayCommand<-' -nowin '}
  if(is.null(saveGraphics)){saveGraphics <- FALSE}
  
  # Define the scenario
  feasibleTasks <- c('saem','SAEM','fim','FIM','ll','LL','graphics','GRAPHICS','indiv','INDIV')
  if(is.null(scenario)){
    displayScenario <- ' -f run'
  }else{
    displayScenario <- NULL
    for(indexScenario in 1:length(scenario)){
      task <- scenario[indexScenario]
      if(is.element(task,feasibleTasks)){
        displayScenario<-paste0(displayScenario,' -f ',scenario[indexScenario])}
      else{
        warning(paste0(task),' in not a valid task for the scenario')
      }
    }
    if(is.null(displayScenario)){
      warning('scenario does not contain any valid task')
      runCalculation <- FALSE
    }
  }
  
  # Change the graphicsToPrint
  if(saveGraphics){
    for(index_project in 1:length(files)){
      filesXmlx <- sub(x = files[index_project],pattern = ".mlxtran",replacement = '_graphics.xmlx')
      if(file.exists(filesXmlx)){
        for(index in 1:length(filesXmlx)){
          fileLines <- readLines(filesXmlx,-1)
          isInGraphics <- 0
          for(ind_Line in (1:1000)){
            # Check if we are in the <graphicsToPrint> Section
            if(length(grep(x = fileLines[ind_Line],pattern = "<graphicsToPrint>"))>0){isInGraphics = 1}
            if(length(grep(x = fileLines[ind_Line],pattern = "</graphicsToPrint>"))>0){isInGraphics = 0}
            if(isInGraphics==1){
              if(length(grep(x = fileLines[ind_Line],pattern = "value=")) >0){
                fileLines[ind_Line] <- sub("0","1",fileLines[ind_Line])
              }
            }
          }
          writeLines(fileLines,paste0(resultDirectory,filesXmlx))
        }
      }
    }
  }
  
  #####################################################################################################
  # Compute all the projects
  #####################################################################################################
  if(runCalculation){
    errorDisplay <- '\n/***************************************************/ \nSummary \n'
    for(index in 1:length(files)){
      base::cat(paste0('Running project ',basename(files[index]),' '))
      out <-system(paste0(monolixPath,' -p ',files[index],displayCommand,displayScenario),show.output.on.console=FALSE,intern = TRUE)
      # Check if there were an error
      indexError <- grep(x=out,pattern="Error message:")
      if(length(indexError)>0){
        errorDisplay <- paste0(errorDisplay,'Project: ',files[index],' NOK \n => ',out[indexError],'\n')
        base::cat(paste0('NOK => ',out[indexError],'\n'))
      }else{
        base::cat('OK\n')
        errorDisplay <- paste0(errorDisplay,'Project: ',files[index],' OK \n')}
    }
    errorDisplay <- paste0(errorDisplay,'/***************************************************/\n')
    base::cat(errorDisplay)
  }
}

