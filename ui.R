library(shiny)
library(shinythemes)
library(shinyalert)
library(shinyvalidate)
library(shinyjs)

options(shiny.sanitize.errors = FALSE)
ui <- navbarPage(useShinyalert(),shinyjs::useShinyjs(),theme = shinytheme('flatly'),
                 tags$head(
                   tags$style(".para {margin-left:200px; margin-right:200px; margin-top:50px ;font-size:16px; text-align: justify}")
                 ),
                 title = a('MicFunPred',href="http://micfunpred.microdm.net.in/",style = "color:white;"),
                 collapsible = T,
                 # About page
                 tabPanel('About',includeMarkdown("about.md")),
                 # Run prediction
                 tabPanel('Run prediction',
                          fluidRow(
                            # test data
                            column(id='col1',3,style = "background-color:#EBF0F1;padding-bottom:20px;",align="center",h3('Test data'),br(),
                                   HTML("<font family='Lato', size=3><b>Abundance table:</b> <a href='test_counts.tsv'>test_counts.tsv</a></font>"),
                                   br(),br(),
                                   HTML("<font family='Lato', size=3><b>ASV/OTU sequence:</b> <a href='test.fasta'>test.fasta</a></font>"),
                                   br(),br(),
                                   HTML("<font family='Lato', size=3><b>Test output:</b> <a href='MicFunPred_test_out.zip'>MicFunPred_test_out.zip</a></font>")
                            ),
                            # prediction form
                            div(id='predictionForm',column(id='col2',3,offset = 1,style = "background-color:#EBF0F1; padding-bottom:20px;",align="center",h3('Run prediction'),fileInput(inputId = 'abundanceTable',label = 'Abundance table (tab separated)',accept = c('text/csv','.txt','.tsv','.csv')),
                                                           fileInput(inputId = 'asvFasta',label = 'ASV/OTU sequences (FASTA format)',accept = c('.fasta','.fna')),
                                                           sliderInput(inputId = 'pident',label = 'Percent identity cut-off',min = 97,max = 100,value = 97,step = 0.5),
                                                           sliderInput(inputId = 'gccutoff',label = 'Gene coverage cut-off',min = 0,max = 1,value = 0.5,step = 0.05),
                                                           textInput(inputId = 'emailid',label = 'Email ID:'),
                                                           actionButton(inputId = 'runPrediction',label = 'Submit job'))),
                            # Download form
                            div(id='downloadForm',column(id='col3',3,offset = 1,style = "background-color:#EBF0F1; padding-bottom:20px;",align="center",h3('Download output'),
                                                         textInput(inputId = 'jobID',label = 'Job ID'),
                                                         downloadButton(outputId = 'downloadJob',label = 'Download output')
                            )),
                            column(3,offset = 1,style = "background-color:#EBF0F1; padding-bottom:20px;",align="center",HTML("<font family='Lato', size=3, color='red'>Server based pradiction is available for input files <=2MB. For larger input files, download and run MicFunPred on local machine</font>
"),h3('Download/Install'),actionButton(inputId='ab1', label="Download/Install", onclick ="window.open('https://github.com/microDM/MicFunPred', '_blank')"))
                          )),
                 tabPanel('Help',includeHTML("help.html")),
                 tabPanel('MicFunPred Team',includeHTML("www/html/contact.html")),
                 tabPanel('Contact',includeMarkdown('contact.md'))
)
