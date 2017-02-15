library(shiny)
library(devtools)
library(logging)
library(shinyjs)
library(shinyBS)

# load supporting methods
source(paste0(packagedir, 'R/GetData.R'))
source(paste0(packagedir, 'R/PlotData.R'))

# load the instance's Rdata, if available
#load(paste0(packagedir, 'data/ObjectWithGEOdata_AS.Rdata'))

DisplayNames <- datasets
names(datasets) <- DisplayNames

genes_list <- read.table(paste0(packagedir, 'data/hgnc_complete_set.txt.gz'), as.is=T, sep='\t',
						                          quote='', comment.char = '', header=T)[,2]

plotHeight=1600

colors=c('red', 'blue', 'purple', 'black')

# Define a server for the Shiny app
shinyServer(function(input, output, session) {
  
  getPlotHeight <- function() {
    return(input$plotHeight)
  }
  getPlotWidth <- function() {
    return(input$plotWidth)
  }
  setPlotHeight <- function(newHeight) {
    updateNumericInput(session, 'plotHeight', value = newHeight)
    #plotHeight=newHeight
  }
  {myheight=1600}
  
  getInputGene <- eventReactive(input$generatePlots, {
    input$gene
  })
  
  observeEvent(input$generatePlots, {
	  closeAlert(session, 'fromURL')
  })

  observeEvent(TRUE, {
	
    query <- parseQueryString(session$clientData$url_search)
    qres <- query$gene

    updateSelectizeInput(session, 'gene', 
                       choices=genes_list,
                       selected=ifelse(is.null(qres), '', qres), 
                       server=TRUE)

    if (!is.null(qres)) {
      disable(id='generatePlots')
	  disable(id='scanDatasets')
      withProgress(message=paste0('!! Finding ', qres, ' in datasets...'),
                   detail='Please wait a moment.', value=0, {
        geneindex <- getSingleGeneIndex(qres, data, bioCdb)
        incProgress(1/3)
        print("updating checkboxes...")
        print(names(data)[sapply(geneindex$indices, length) > 0])
        sets_with_gene <- names(data)[sapply(geneindex$indices, length) > 0]
        incProgress(2/3)
        updateCheckboxGroupInput(session, "datasets", choices=datasets,
                                 selected=sets_with_gene)
        incProgress(3/3)
      })
	  alertText = paste0('Greetings, the datasets containing ', qres, ' have been selected.')
	  createAlert(session, 'alert', 'fromURL', title=alertText,
	              content='Please click "Generate Plots"', append=FALSE)
	  enable(id='generatePlots')
	  enable(id='scanDatasets')
    }
  })


  observeEvent(input$scanDatasets, {
    disable(id='generatePlots')
    withProgress(message=paste0('!! Finding ', input$gene, ' in datasets...'),
                 detail='Please wait a moment.', value=0, {
      geneindex <- getSingleGeneIndex(input$gene, data, bioCdb)
      incProgress(1/3)
      print("updating checkboxes...")
      print(names(data)[sapply(geneindex$indices, length) > 0])
      sets_with_gene <- names(data)[sapply(geneindex$indices, length) > 0]
      incProgress(2/3)
      updateCheckboxGroupInput(session, "datasets", choices=datasets, 
                               selected=sets_with_gene)
      incProgress(3/3)
    })
    enable(id='generatePlots')
  })
  
  
  output$plot <- renderPlot({
    withProgress(message=paste0('Building ', input$method, ' plots...'),
                 detail='Please wait a moment.', value=0, {
      gene <- getInputGene()
      incProgress(1/5)
      AtoB <- match(input$datasets, names(data))
      incProgress(2/5)
      geneindex <- getSingleGeneIndex(gene, data[AtoB], bioCdb)
      numPlots <- length(unlist(geneindex))
      print(paste0('numPlots: ', numPlots))
      setPlotHeight(ceiling(numPlots / 2) * 400)
      incProgress(3/5)
      # Render a boxplot
      if(input$method == 'Box') {
            pp <- (boxOneGene(data[AtoB], geneindex, groups[AtoB], 
                          main=paste('Boxplot with combined individuals and probes',gene, sep='\n'), 
                          colors=colors, ylab='(Raw data from GEO)', notch=T))
            incProgress(4/5)
      }
      if(input$method == 'Dot') {
            pp <- (stripOneGene(data[AtoB], geneindex, groups[AtoB], 
                            main=paste('Dotplot for each individual and probe',gene, sep='\n'), 
                            colors=colors, ylab='(Raw data from GEO)'))
            incProgress(4/5)
      }
      if(input$method == 'Violin') {
            pp<-violinOneGene(data[AtoB], geneindex, groups[AtoB], ylab='(Raw data from GEO)',  
                          main=paste('All Samples and Probes',gene, sep='\n'), colors=colors)
            incProgress(4/5)
      }
      incProgress(5/5)
      if(exists('pp')) return(pp)
    })    
  #}, width=800, height=2400)
  }, width=800, height=function(){getPlotHeight()})

  
    
  output$text <- renderDataTable({
      gene <- getInputGene()
      AtoB <- match(input$datasets, names(data))
      geneindex <- getSingleGeneIndex(gene, data[AtoB], bioCdb)
      
      if(input$method == 'Long') {
        xx<-tableOneGene(data[AtoB], geneindex, groups[AtoB])
      }
      if(input$method == 'Ragged') {
        xx<-raggedOneGene(data[AtoB], geneindex, groups[AtoB])
      }
      if(grepl('Grouping', input$method)) {
        xx<-groupingDataset(data, strsplit(input$method, '_')[[1]][2], groups)
      }
      if(exists('xx')) return(xx)
    },  
      options = list(
        pageLength=-1,
        lengthMenu = list(c(50, 100, -1), c('50', '100', 'All'))
      )
  )

  
})


