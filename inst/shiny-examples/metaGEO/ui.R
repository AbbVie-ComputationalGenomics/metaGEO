library(shiny)
library(shinyjs)
library(shinyBS)
#source('global.R')

DisplayNames <- datasets
names(datasets) <- DisplayNames

# Define the overall UI
shinyUI(
  fluidPage(
    useShinyjs(),
    extendShinyjs(text=jsCode, functions='flip_checkbox'),
    # Give the page a title
    titlePanel("metaGEO: Plot gene expression across multiple GEO datasets"),
	  h2(metageo_title),
    # Define the sidebar with two inputs
    sidebarLayout(
      sidebarPanel(
        "Instructions:",
        tags$ol(
          tags$li("Enter HuGO gene symbol."),
          tags$li("Click 'Scan Datasets' to filter out sets with zero matching probes."),
          tags$li("Click 'Generate Plots'.")
        ),
        hr(),
        selectizeInput('gene', "Gene Symbol (exact HuGO nomenclature)", choices=NULL),
		    actionButton("scanDatasets", "Scan Datasets", icon("search"),
		                 style="color: #fff; background-color: #337ab7; border-color: #2e6da4; margin-bottom: 10px; margin-top: -3px;"),
        checkboxGroupInput("datasets", "Datasets containing expression data for the selected gene:",
                           choices=datasets, selected=datasets),
        # Copy the line below to make a date selector
        selectizeInput("method", "Display Method:",
                       choices=c(plotmethods, textmethods), list(placeholder = 'Select an output type')),
        actionButton("generatePlots", "Generate Plots", icon("bar-chart"), 
                     style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
		    hr(),
		    helpText("Configurable multi-plot dimensions. Re-rendering may be slow."),
		    fluidRow(
		      column(6, numericInput("plotHeight", label="Height", value=2400, step=100, min=600, max=10000, width='120px')),
		      column(6, numericInput("plotWidth", label="Width", value=800, step=50, min=400, max=1600, width='120px'))
		    )
      ),

      # Create a spot for the barplot
      mainPanel(
	    bsAlert("alert"),
        plotOutput("plot"),
        dataTableOutput('text')
      )
    )
  )
)

