#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
#library(shinyWidgets)
#library(shinydashboard)
library(ggplot2)
library(plotly)
library(DT)
library(shinybusy)

# Define UI for application that draws a histogram
shinyUI(
    tagList(
        #useShinydashboard(),
        busy_start_up(
            loader = spin_epic("orbit", color ="#FFB81C"), #"#FFF"
            #text = "Loading... \n aa",
            text = HTML(paste("Bethan Data", "Loading...", sep = '<br/>')),
            timeout = 3500,
            color = "#FFF",
            background = "#003B5C"   #"#112446"
        ),
        
        includeCSS("www/mycss.css"),
        
    #shinythemes::themeSelector(),
    navbarPage(fluid=T, windowTitle = "Bethan Data",
        theme = shinytheme("cosmo"),  # <--- To use a theme, uncomment this
        "Bethan Data",
        tabPanel("Navbar 1",#fluid=FALSE,
                 sidebarPanel(fluid=FALSE, width = 2,
                     #fileInput("spList", "Species List:", accept = c(".xlsx")),
                     #textInput("sp", "Species:", "LPC 16:1"),
                     
                     #tags$h5("Search Gene"),
                     #textInput("gene", "Search Gene", "Apoh"),
                     selectizeInput("gene", "Search Gene", 
                                 choices=c("Apoh"), 
                                 selected = "Apoh"),
                     actionButton("search_gene.lipid.corr", "Search"),
                     
                     #tags$h5("Download Data:"),
                     
                     
                     #sliderInput("slider", "Label input:", 1, 24, 12),
                     #sliderInput("slider.D", "D input:", 0,100,50),
                     #sliderInput("slider.S", "S input:", 0,100,50),
                     #tags$h5("Default actionButton:"),
                     #actionButton("action", "Search"),
                     
                     #tags$h5("actionButton with CSS class:"),
                     #actionButton("action2", "Action button", class = "btn-primary")
                 ),
                 mainPanel(width=9,
                     tabsetPanel(
                         tabPanel("Gene-Lipid Corr",
                                  
                                  h3("Gene Boxplot"),
                                  plotlyOutput("gene.boxplot",height = 600),
                                  hr(),
                                  
                                  h3("Correlation Table"),
                                  DT::dataTableOutput("corr.table", height = "450px"),
                                  hr(),
                                  
                                  h3("Scatter Plot(selected in table)"),
                                  plotlyOutput("genelip.scatplot",height = 500),
                                  hr(),
                     
                                  #h3("Scatter Plot(selected in table)"),
                                  #plotlyOutput("genelip.scatplot"),
                                  #fluidRow(class = "myRow1", 
                                  #    splitLayout(cellWidths = c("55%", "45%"), 
                                  #                plotlyOutput("gene.boxplot",height = 500),
                                  #                plotlyOutput("genelip.scatplot",height = 500)),
                                  #),
                                  
                                  
                                  #h4("row"),
                                  #verbatimTextOutput("spInput"),
                                  #h1("Header 1"),
                                  #h2("Header 2"),
                                  #h3("Header 3"),
                                  #h4("Header 4"),
                                  #h5("Header 5")
                         ),
                         
                         ##tab boxplot
                         #tabPanel("Gene Boxplot", 
                         #         #plotlyOutput("gene.boxplot")
                                  
                         #         ),
                        
                         tabPanel("RNA-Seq Data",
                                  h3("RNA-Seq Data"),
                                  DT::dataTableOutput("rnaseq")),
                         
                         tabPanel("Lipid Data", 
                                  h3("Lipid Data"),
                                  DT::dataTableOutput("lipidtable"),
                                  hr(),
                                  
                                  h3("Lipid Box Plot(selected in table)"),
                                  plotlyOutput("Lipid.boxplot")
                                  ),
                         
                         tabPanel("Tab 1", "This panel is intentionally left blank")
                     )
                 )
        ),
        tabPanel("Navbar 2", "This panel is intentionally left blank"),
        tabPanel("Navbar 3", "This panel is intentionally left blank")
    )
)
)
