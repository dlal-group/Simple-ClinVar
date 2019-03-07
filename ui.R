library(shiny) 
library(ggplot2) 
library(readr) 
library(stringr)
library(tidyr)
library(Hmisc)
library(ggrepel)
library(DT)
library(ggrepel)
library(dplyr)
library(tm)
library(shinycssloaders)
library("shinythemes")
# With Drawprotein
library(drawProteins)
library(magrittr)
#library(rlang)


# Define UI for application that draws a histogram
shinyUI(fluidPage(
  navbarPage(theme = shinytheme("yeti"), id = "inTabset",selected = "panel1",
             title = "SCV",
             # Main box search and description -----------------------------------------
             tabPanel("Home", value = "panel1",
                      br(),br(),
                      h2("Shiny ClinVar (SCV)", align = "center"),
                      br(),br(),
                      fluidRow(
                        column(width = 6, offset = 3, align = "center",
                               wellPanel("Type gene name or disease term",
                                         style = "background-color: #333333;
                                         color: white;
                                         border-top-color: #333333;
                                         border-left-color: #333333;
                                         border-right-color: #333333;
                                         box-shadow: 3px 3px 3px #d8d8d8;
                                         margin-bottom: 0px;
                                         padding:5px"), 
                               wellPanel(br(),br(),
                                         textInput(inputId = "genename", label = NULL),#, value = "clinvar"),
                                         br(),
                                         actionButton(inputId ="goButton", label = "Submit", class = "btn-primary"), 
                                         style = "background-color: #ffffff;
                                         border-bottom-color: #333333;
                                         border-left-color: #333333;
                                         border-right-color: #333333;
                                         box-shadow: 3px 3px 3px #d8d8d8;
                                         margin-top: 0px")
                               ) # WellPanel
                               ), #Fluid row
                      fluidRow(column(width = 6, offset = 3, br(), br(), p("How many missense variants are associated to heart disease?
                                                                           What are the top 10 genes mutated in Alzheimer? Does CDKL5 have pathogenic mutations? If so, where?
                                                                           SCV is able to answer these questions and more, in a matter of seconds. \n SCV was developed to 
                                                                           provide gene and disease wise summary statistic based on all available genetic variants from ClinVar.",
                                                                           align = "center"), p(""), style = "background-color: #ffffff")
                      )
                      ),
             
             # Filtering "left" side -------------------------------------------------------------
             tabPanel("Results", value = "panel2",
                      fluidRow(column(width = 12, h6(strong("Displaying results for:")))), 
                      fluidRow(column(width = 2, style = "padding-right: 0px",
                                      wellPanel(textOutput(outputId = "displayquery"),
                                                style = "background-color:#008CBA;
                                                color:white;
                                                border-color:#008CBA;
                                                box-shadow: 3px 3px 3px #d8d8d8;
                                                margin-bottom: 0px;
                                                padding:15px;
                                                width: 80%"),
                                      fluidRow(column(width = 12,h5(""))),
                                      fluidRow(column(width = 12,h5(""))),
                                      wellPanel(h6(strong("Choose your filter")),
                                                style = "background-color:white;
                                                border-bottom: 2px solid #EEEEEE;
                                                border-top-color: white;
                                                border-right-color: white;
                                                border-left-color: white;
                                                box-shadow: 0px 0px 0px white;
                                                padding:3px;
                                                width: 100%"),
                                      checkboxGroupInput(inputId = "clinicalinput",
                                                         label = "Clinical Significance",
                                                         choices = c("Pathogenic",
                                                                     "Likely pathogenic",
                                                                     "Risk factor/Association",
                                                                     "Uncertain/conflicting",
                                                                     "Protective/Likely benign",
                                                                     "Benign"),
                                                         selected = NULL ),
                                      checkboxGroupInput(inputId = "consequenceinput",
                                                         label = "Consequence",
                                                         choices= c("Missense",
                                                                    "Stop-gain",
                                                                    "In frame indel",
                                                                    "Frameshift",
                                                                    "Synonymous"),
                                                         selected = NULL ),
                                      checkboxGroupInput(inputId = "typeinput",
                                                         label = "Type",
                                                         choices= c("SNV",
                                                                    "Deletion",
                                                                    "Duplication",
                                                                    "Deletion",
                                                                    "Indel",
                                                                    "Inversion"),
                                                         selected = NULL ),
                                      checkboxGroupInput(inputId = "reviewinput",
                                                         label = "Review status",
                                                         choices= c("Criteria provided/ multiple submitters/ no conflicts",
                                                                    "Reviewed by expert panel",
                                                                    "Practice guideline",
                                                                    "Criteria provided/ single submitter",
                                                                    "Criteria provided/ conflicting interpretations",
                                                                    "No assertion criteria provided"),
                                                         selected = NULL ),
                                      actionButton(inputId = "filter", label = "Filter", width = "80%", 
                                                   style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                   margin-bottom: 0px;
                                                   padding:5px"),
                                      br(),br(),
                                      actionButton(inputId = "newsearch", label = "New Search", 
                                                   width = "80%",
                                                   style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                   margin-bottom: 0px;
                                                   padding:5px")
                                      ),
                               
                               # Results/Plots Panel -----------------------------------------------------
                               column(width = 10, style = "padding-left: 0px",
                                      fluidRow(
                                        column(2, offset = 1, 
                                               actionButton(inputId = "seevariants", label = textOutput(outputId = "displayobservation"),
                                                            width = "100%",
                                                            style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                            margin-bottom: 0px;
                                                            padding:15px",
                                                            class = "btn-success")),
                                        column(2, offset = 1,
                                               actionButton(inputId = "seegenes", label = textOutput(outputId = "displayngenes"),
                                                            width = "100%",
                                                            style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                            margin-bottom: 0px;
                                                            padding:15px",
                                                            class = "btn-danger")),
                                        column(2, offset = 1,
                                               actionButton(inputId = "seephenotypes", label = textOutput(outputId = "displayphenotypes"),
                                                            width = "100%",
                                                            style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                            margin-bottom: 0px;
                                                            padding:15px",
                                                            class = "btn-warning")),
                                        column(2, offset = 1,  
                                               actionButton(inputId = "seetable", label = "Show Table",
                                                            width = "100%",
                                                            style = "box-shadow: 3px 3px 3px #d8d8d8;
                                                            margin-bottom: 0px;
                                                            padding:15px"))
                                               ),
                                      br(),br(),
                                      fluidRow(column(width = 12,
                                                      tabsetPanel(type  = "tabs", id = "displaywhat", selected = "1",
                                                                  tabPanel(title = NULL, value = "1",#br(),
                                                                           fluidRow(
                                                                             column(3, h6("\tType"), align = "center",
                                                                                    withSpinner(plotOutput("Type", width="100%"))),
                                                                             column(3, h6("\tConsequence"), align = "center",
                                                                                    withSpinner(plotOutput("consequence", width="100%"))),
                                                                             column(3, h6("\tClinical Significance"), align = "center",
                                                                                    withSpinner(plotOutput("ClinicalSignificance", width="100%"))),
                                                                             column(3, h6("\tReview status"), align = "center",
                                                                                    withSpinner(plotOutput("review", width="100%")))
                                                                           )
                                                                  ), 
                                                                  tabPanel(title = NULL, value = "2", 
                                                                           uiOutput("genefield")
                                                                  ),
                                                                  tabPanel(title = NULL, value = "3",#br(),
                                                                           fluidRow(column(8, h6("Top 10 phenotypes associated"), align = "center", br(),
                                                                                           plotOutput("DiseaseCount")),
                                                                                    column(4, h6("Total phenotypes associations"), align = "center",
                                                                                           wellPanel(DT::dataTableOutput("DiseaseCountTotal"),
                                                                                                     style = "background-color: #ffffff;
                                                                                                     border-color: #ffffff;
                                                                                                     box-shadow: 0px 0px 0px #ffffff;
                                                                                                     margin-bottom: 5px")
                                                                                           )
                                                                                           ) #Fluid row
                                                                                           ),
                                                                  tabPanel(title = NULL, value = "4",
                                                                           fluidRow(column(12,br(), align = "center",
                                                                                           wellPanel(DT::dataTableOutput("clinvartable1"),
                                                                                                     style = "background-color: #ffffff;
                                                                                                     border-color: #ffffff;
                                                                                                     box-shadow: 0px 0px 0px #ffffff;
                                                                                                     margin-bottom: 5px")
                                                                                           ) # Column
                                                                                           ), 
                                                                           fluidRow(column(12,downloadLink("downloadData1", "Download"))
                                                                                           ), 
                                                                           fluidRow(column(12,br(), align = "center",
                                                                                          wellPanel(DT::dataTableOutput("clinvartable2"),
                                                                                                   style = "background-color: #ffffff;
                                                                                                   border-color: #ffffff;
                                                                                                   box-shadow: 0px 0px 0px #ffffff;
                                                                                                   margin-bottom: 5px")
                                                                                                              )# Column
                                                                                                              ),# Fluid row
                                                                           fluidRow(column(12,downloadLink("downloadData2", "Download"))
                                                                           )
                                                                                           ) #TabPanle
                                                                           ) #TabsetPanel
                                                                           ) # Column
                                                                  ) #Fluidrow aqui abajo solo de bonito
                                                                                    ) # column
                                                      ),
                      fluidRow(column(width = 12,
                                      wellPanel(style = "background-color:white;
                                                border-bottom: 2px solid #EEEEEE;
                                                border-top-color: white;
                                                border-right-color: white;
                                                border-left-color: white;
                                                box-shadow: 0px 0px 0px white;
                                                padding:0px;
                                                 width: 100%"))) #,
                     # verbatimTextOutput("test")
                                      ), #mainTabPanel,
             
             # about page --------------------------------------------------------------
             tabPanel("About", value = "panel3",                  
                      br(),br(),
                      fluidRow(column(width = 5, offset = 1, align = "center",
                                      
                                      wellPanel("Beta Release",
                                                style = "background-color: #333333;
                                                color: white;
                                                border-top-color: #333333;
                                                border-left-color: #333333;
                                                border-right-color: #333333;
                                                box-shadow: 3px 3px 3px #d8d8d8;
                                                margin-bottom: 0px;
                                                padding:5px"), 
                                      wellPanel(br(),p("The present release of Shiny ClinVar is a beta version."),
                                                p("We encourage you to contact us with bugs or suggestions."),
                                                br(),
                                                style = "background-color: #ffffff;
                                                border-color: #333333;
                                                box-shadow: 3px 3px 3px #d8d8d8;
                                                margin-top: 0px")
                                      ), # Column
                               column(width = 5, 
                                      wellPanel("Contact", align = "center",
                                                style = "background-color: #333333;
                                                color: white;
                                                border-top-color: #333333;
                                                border-left-color: #333333;
                                                border-right-color: #333333;
                                                box-shadow: 3px 3px 3px #d8d8d8;
                                                margin-bottom: 0px;
                                                padding:5px"), 
                                      wellPanel(align = "center",
                                                p(br(),icon("envelope", lib = "glyphicon"),"   e.perezpalma@uni-koeln.de"),
                                                p(icon("envelope", lib = "glyphicon"),"   mgramm1@smail.uni-koeln.de"),
                                                br(),
                                                style = "background-color: #ffffff;
                                                border-color: #333333;
                                                box-shadow: 3px 3px 3px #d8d8d8;
                                                margin-top: 0px")
                                      ) # Column
                                      ),br(),br() #Fluid row
                      # fluidRow(column(width = 6, offset = 3, align = "center",
                      #                 wellPanel("Tutorial",
                      #                           style = "background-color: #333333;
                      #                           color: white;
                      #                           border-top-color: #333333;
                      #                           border-left-color: #333333;
                      #                           border-right-color: #333333;
                      #                           box-shadow: 3px 3px 3px #d8d8d8;
                      #                           margin-bottom: 0px;
                      #                           padding:5px"), 
                      #                 wellPanel(br(),br(),br(),
                      #                           p(">youtube placeholder<"),br(),br(),br(),
                      #                           style = "background-color: #ffffff;
                      #                           border-color: #b2baba;
                      #                           box-shadow: 3px 3px 3px #d8d8d8;
                      #                           margin-top: 0px")
                      #                 ) #column
                      #                 ) #Fluid row
                                      ) # TabPanel
                                      ) # NavbarPage
                      ))
