library(shiny)
library(shinythemes)
library(shinyalert)
library(proceduralnames)
library(shinyjs)
library(emayili)
library(shinyvalidate)

jscode <- "shinyjs.refresh = function() { history.go(0); }"

# shiny options
options(shiny.sanitize.errors = FALSE)
options(shiny.maxRequestSize=2*1024^2) 
options(shiny.sanitize.errors = FALSE)

server <- function(input, output, session){
	shinyOptions(shiny.maxRequestSize=1*1024^2)
  iv <- InputValidator$new()
  iv$add_rule("abundanceTable", sv_required())
  iv$add_rule("asvFasta", sv_required())
  iv$add_rule("emailid", sv_email())
  iv$enable()
  appDir <- '/opt/shiny-server-apps/MicFunPred'
  rvalues <- reactiveValues(abundanceTable=NULL,asvFasta=NULL,pident=NULL,gccutoff=NULL,emailid=NULL,jobID=NULL)
  # 1. run prediction
  observeEvent(input$runPrediction,{
    if(!(is.null(input$abundanceTable)|is.null(input$asvFasta))){
      # 1. Assign jobID and make directory
      jobID <- paste(make_english_names(n = 1,n_words = 2),paste(format(Sys.time(),'%d'),format(Sys.time(),'%m'),format(Sys.time(),'%y'),sep=''),sep='_')
      dir.create(path = file.path(appDir,'jobs',jobID),mode = '777')
      system(paste('chmod 777',file.path(appDir,'jobs',jobID)))
      # 2. copy input files
     file.copy(from = input$abundanceTable$datapath,to = file.path(appDir,'uploads',paste(jobID,'abundance.txt',sep = '_')))
     file.copy(from = input$asvFasta$datapath,to = file.path(appDir,'uploads',paste(jobID,'seq.fna',sep = '_')))
      # 3. Alert
      shinyalert(title = 'Success. Sit tight now. We will let you know.',text = paste('Once job is completed, we will mail you at ',input$emailid,'\n','Job ID:',jobID))
      shinyjs::reset('predictionForm')
      # 4. run MicFunPred
      system(paste('touch',file.path(appDir,'jobs',jobID,'log.txt')))
      system(paste('chmod 777',file.path(appDir,'jobs',jobID,'log.txt')))
      # change directory to jobDIr
      jobDir <- file.path(appDir,'jobs',jobID)
      setwd(jobDir)
      abundaceFilePath <- file.path(appDir,'uploads',paste(jobID,'abundance.txt',sep = '_'))
      asvFilePath <- file.path(appDir,'uploads',paste(jobID,'seq.fna',sep = '_'))
      cmd <- paste(file.path(appDir,'pyenv/bin/MicFunPred_run_pipeline.py'),'-i',abundaceFilePath,'-r',asvFilePath,'-p',input$pident,'-c',input$gccutoff,'-o',jobDir,'-t',4,'-v','-j',jobID,'-e',input$emailid)
      system(cmd,wait = FALSE)
    }
  })
  # 2. Download output
  output$downloadJob <- downloadHandler(
    filename <- function() {
      paste(input$jobID, "tar.gz", sep=".")
    },
    
    content <- function(file) {
      setwd(file.path(appDir,'jobs'))
      if(dir.exists(input$jobID)){
        tar(file, input$jobID,compression = 'gzip',compression_level = 9) 
      }
      else{
        shinyalert(title='Your job is still running or not found')
      }
    }
  ) 
} 

