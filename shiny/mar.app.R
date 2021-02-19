#
# This is a Shiny web application. 
#

# Load libraries
library(shiny)
library(shinyjs)
library(ggplot2)
library(plyr)
library(dplyr)
library(data.table)

# # Load functions
# source("helpers.R")
# 
# # load file with codes
# human <- fread("data/listGenes/listHumanGenes.csv", header = FALSE)
# macaque <- fread("data/listGenes/listMacaqueGenes.csv", header = FALSE)
# mouse <- fread("data/listGenes/listMouseGenes.csv", header = FALSE)
# rat <- fread("data/listGenes/listRatGenes.csv", header = FALSE)
# rabbit <- fread("data/listGenes/listRabbitGenes.csv", header = FALSE)
# opossum <- fread("data/listGenes/listOpossumGenes.csv", header = FALSE)
# chicken <- fread("data/listGenes/listChickenGenes.csv", header = FALSE)
# 
# listGenes <- list(human, macaque, mouse, rat, rabbit, opossum, chicken)
# 
# ## Load data - for human
# humanLabels <- fread("data/humanLabels.csv", header = FALSE)
# #humanGenes <- fread("data/humanGenes.txt", header = FALSE)
# humanGenes <- fread("data/humanRounded.csv", header = FALSE)
# humanOrder <- fread("data/humanOrder.csv", header = FALSE)
# 
# ## Load data - for mouse
# mouseLabels <- fread("data/mouseLabels.csv", header = FALSE)
# #mouseGenes <- fread("data/mouseGenes.txt", header = FALSE)
# mouseGenes <- fread("data/mouseRounded.csv", header = FALSE)
# mouseOrder <- fread("data/mouseOrder.csv", header = FALSE)
# 
# ## Load data - for rat
# ratLabels <- fread("data/ratLabels.csv", header = FALSE)
# #ratGenes <- fread("data/ratGenes.txt", header = FALSE)
# ratGenes <- fread("data/ratRounded.csv", header = FALSE)
# ratOrder <- fread("data/ratOrder.csv", header = FALSE)
# 
# ## Load data - for opossum
# opossumLabels <- fread("data/opossumLabels.csv", header = FALSE)
# #opossumGenes <- fread("data/opossumGenes.txt", header = FALSE)
# opossumGenes <- fread("data/opossumRounded.csv", header = FALSE)
# opossumOrder <- fread("data/opossumOrder.csv", header = FALSE)
# 
# ## Load data - for rabbit
# rabbitLabels <- fread("data/rabbitLabels.csv", header = FALSE)
# #rabbitGenes <- fread("data/rabbitGenes.txt", header = FALSE)
# rabbitGenes <- fread("data/rabbitRounded.csv", header = FALSE)
# rabbitOrder <- fread("data/rabbitOrder.csv", header = FALSE)
# 
# ## Load data - for macaque
# macaqueLabels <- fread("data/macaqueLabels.csv", header = FALSE)
# #macaqueGenes <- fread("data/macaqueGenes.txt", header = FALSE)
# macaqueGenes <- fread("data/macaqueRounded.csv", header = FALSE)
# macaqueOrder <- fread("data/macaqueOrder.csv", header = FALSE)
# 
# ## Load data - for chicken
# chickenLabels <- fread("data/chickenLabels.csv", header = FALSE)
# #chickenGenes <- fread("data/chickenGenes.txt", header = FALSE)
# chickenGenes <- fread("data/chickenRounded.csv", header = FALSE)
# chickenOrder <- fread("data/chickenOrder.csv", header = FALSE)
# 
# # Define pallet
colTissue <- c(Brain = "#3399CC", Cerebellum = "#33CCFF", Heart = "#CC0000", Kidney = "#CC9900", Liver = "#339900", Ovary = "#CC3399", Testis = "#FF6600")

# Vector of species
species <- c("Human", "Macaque", "Mouse", "Rat", "Rabbit", "Opossum", "Chicken")

# Define UI for application
ui <- fluidPage(
  
  # link css styles
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
  # link font
  tags$head(tags$link(rel = "stylesheet", type = "text/css", href = "https://fonts.googleapis.com/css?family=Quicksand")),
  # link javascript file
  tags$head(tags$script(src = "scripts.js")),

  useShinyjs(),

  absolutePanel(
    top = 0, left= 0, right = 0,
    fixed = TRUE,
    tags$div(
      id = "topBar",
      style="padding: 8px 8px 0 26px; border-bottom: 1px solid #CCC; background: #FFF;",
      tags$img(id = "mainPicture", src = "Mammalian_evolutionTransparent.png",
               width = "150px", style = "margin-left: 10px"),
      titlePanel("Evo-devo mammalian organs"),
      ##
      tags$div(class = "pbusy",
              tags$img(id = "spinner", src="ajax-loader.gif", width = "80px"),
              p("Loading...")
      ),
      ##
      htmlOutput("gene_title")
    )),
  
  sidebarLayout(
    sidebarPanel(width = 3,
      tags$div(class = "side_panel",
        # Panel for the search options
        # "gene" - Create input of type textInput() and submittButton (in the future with autocomplete) 
        # for specifying a gene and selecting it
        fluidRow(textInput("gene", label = h4("Gene"), value = NULL, placeholder = "type gene or ensembl id"),
                 actionButton("click", label = "Search")),
        fluidRow(
          # "species" - Create input of type checkboxGroupInput() for selection of species,
          column(6, checkboxGroupInput("species", label = h4("Species:"),
                                       choices = species,
                                       selected = species)),
          # "tissues" - Create input of type checkboxGroupInput() for selection of tissues
          column(6, checkboxGroupInput("tissues", label = h4("Tissue:"),
                                       choices = names(colTissue),
                                       selected = names(colTissue))),
          column(12, actionButton("update", label = "Update"))
          ),
        # download button
        fluidRow(column(12, downloadButton(outputId = "report", label = "Download plots"))),
        fluidRow(column(12, p("Download expression tables (cpm)")),
                 column(9, selectInput(inputId = "speciesData", label = NULL, choices = species)),
                 column(3, downloadButton(outputId = "downloadData", label = "")))
        )
    ),
    
    # Panel for the plots
    mainPanel(width = 9,
      tags$div(id = "to_hide",
        # titlePanel("Evo-devo mammalian organs"),
        # htmlOutput("gene_title"),
        fluidRow(tags$div(uiOutput("human"),
                        plotOutput("humanPlot")),
                 tags$div(uiOutput("macaque"),
                        plotOutput("macaquePlot")),
                 tags$div(uiOutput("mouse"),
                        plotOutput("mousePlot")),
                 tags$div(uiOutput("rat"),
                        plotOutput("ratPlot")),
                 tags$div(uiOutput("rabbit"),
                        plotOutput("rabbitPlot")),
                 tags$div(uiOutput("opossum"),
                        plotOutput("opossumPlot")),
                 tags$div(uiOutput("chicken"),
                        plotOutput("chickenPlot"))
        )
      )
    )
  )
)

# Define server logic required to draw the plots
server <- function(input, output) {
  
  # Get information from the ensemble API
  geneTitle <- eventReactive(input$click, {
    geneT <- getEnsembl(input$gene, listGenes)
    return(geneT)
  })
  
  geneT <- eventReactive(input$click, {
    if (any(geneTitle() == 0)) {
      return("Sorry, that gene is invalid.")
    }
    if (geneTitle()["human","display_name"][1] != "CHARACTER(0)") {
      return(geneTitle()["human","display_name"][1])
    }
    if (geneTitle()["macaque","display_name"][1] != "CHARACTER(0)") {
      return(geneTitle()["macaque","display_name"][1])
    }
    if (geneTitle()["mouse","display_name"][1] != "CHARACTER(0)") {
      return(geneTitle()["mouse","display_name"][1])
    }
    if (geneTitle()["rat","display_name"][1] != "CHARACTER(0)") {
      return(geneTitle()["rat","display_name"][1])
    }
    if (geneTitle()["opossum","display_name"][1] != "CHARACTER(0)") {
      return(geneTitle()["opossum","display_name"][1])
    }
    if (geneTitle()["rabbit","display_name"][1] != "CHARACTER(0)") {
      return(geneTitle()["rabbit","display_name"][1])
    }
    if (geneTitle()["chicken","display_name"][1] != "CHARACTER(0)") {
      return(geneTitle()["chicken","display_name"][1])
    }
  })

  # Create a reactive expression to figure out the gene ensembl for all the species
  # Returns list with gene data for each species
  geneEns <- eventReactive(input$click, {
    ensembl <- geneTitle()
    if (any(ensembl == 0)) {
      return(3)
    }
    else {
      # make datatable for each species and bind them in a list
      human <- geneInfo(ensembl["human","id"], humanGenes, humanLabels)
      mouse <- geneInfo(ensembl["mouse","id"], mouseGenes, mouseLabels)
      rat <- geneInfo(ensembl["rat","id"], ratGenes, ratLabels)
      rabbit <- geneInfo(ensembl["rabbit","id"], rabbitGenes, rabbitLabels)
      opossum <- geneInfo(ensembl["opossum","id"], opossumGenes, opossumLabels)
      macaque <- geneInfo(ensembl["macaque","id"], macaqueGenes, macaqueLabels)
      chicken <- geneInfo(ensembl["chicken","id"], chickenGenes, chickenLabels)
      speciesList <- list(human, mouse, rat, rabbit, opossum, macaque, chicken)
      names(speciesList) <- c("human", "mouse", "rat", "rabbit", "opossum", "macaque", "chicken")

      # if HUMAN is selected
      toggle(id = "humanPlot", condition = ("Human" %in% input$species) & (is.data.frame(speciesList[["human"]])) & (any(input$species > 0)))
      toggle(id = "human", condition = "Human" %in% input$species)
      
      # if MOUSE is selected
      toggle(id = "mousePlot", condition = ("Mouse" %in% input$species) & (is.data.frame(speciesList[["mouse"]])) & (any(input$species > 0)))
      toggle(id = "mouse", condition = "Mouse" %in% input$species)
      
      # if RAT is selected
      toggle(id = "ratPlot", condition = ("Rat" %in% input$species) & (is.data.frame(speciesList[["rat"]])) & (any(input$species > 0)))
      toggle(id = "rat", condition = "Rat" %in% input$species)
      
      # if OPOSSUM is selected
      toggle(id = "opossumPlot", condition = ("Opossum" %in% input$species) & (is.data.frame(speciesList[["opossum"]])) & (any(input$species > 0)))
      toggle(id = "opossum", condition = "Opossum" %in% input$species)
      
      # if RABBIT is selected
      toggle(id = "rabbitPlot", condition = ("Rabbit" %in% input$species) & (is.data.frame(speciesList[["rabbit"]])) & (any(input$species > 0)))
      toggle(id = "rabbit", condition = "Rabbit" %in% input$species)
      
      # if MACAQUE is selected
      toggle(id = "macaquePlot", condition = ("Macaque" %in% input$species) & (is.data.frame(speciesList[["macaque"]])) & (any(input$species > 0)))
      toggle(id = "macaque", condition = "Macaque" %in% input$species)
      
      # if CHICKEN is selected
      toggle(id = "chickenPlot", condition = ("Chicken" %in% input$species) & (is.data.frame(speciesList[["chicken"]])) & (any(input$species > 0)))
      toggle(id = "chicken", condition = "Chicken" %in% input$species)

      return(speciesList)
    }
  })
  
  # Create a reactive expression that subsets the tissues from the species list
  geneTissues <- eventReactive(input$update, {
    geneData <- geneEns()
    selectedTissues <- tissues(geneData,input$tissues)
    return(selectedTissues)
  })
  
  # check selected species and hide the ones that aren't selected
  observeEvent(input$update, {
    if (any(geneTitle() != 0)) {
      geneE <- geneEns()
      # if HUMAN is selected
      toggle(id = "humanPlot", condition = ("Human" %in% input$species) & (is.data.frame(geneE[["human"]])) & (any(input$species > 0)))
      toggle(id = "human", condition = "Human" %in% input$species)
  
      # if MOUSE is selected
      toggle(id = "mousePlot", condition = ("Mouse" %in% input$species) & (is.data.frame(geneE[["mouse"]])) & (any(input$species > 0)))
      toggle(id = "mouse", condition = "Mouse" %in% input$species)
  
      # if RAT is selected
      toggle(id = "ratPlot", condition = ("Rat" %in% input$species) & (is.data.frame(geneE[["rat"]])) & (any(input$species > 0)))
      toggle(id = "rat", condition = "Rat" %in% input$species)
  
      # if OPOSSUM is selected
      toggle(id = "opossumPlot", condition = ("Opossum" %in% input$species) & (is.data.frame(geneE[["opossum"]])) & (any(input$species > 0)))
      toggle(id = "opossum", condition = "Opossum" %in% input$species)
  
      # if RABBIT is selected
      toggle(id = "rabbitPlot", condition = ("Rabbit" %in% input$species) & (is.data.frame(geneE[["rabbit"]])) & (any(input$species > 0)))
      toggle(id = "rabbit", condition = "Rabbit" %in% input$species)
  
      # if MACAQUE is selected
      toggle(id = "macaquePlot", condition = ("Macaque" %in% input$species) & (is.data.frame(geneE[["macaque"]])) & (any(input$species > 0)))
      toggle(id = "macaque", condition = "Macaque" %in% input$species)
  
      # if CHICKEN is selected
      toggle(id = "chickenPlot", condition = ("Chicken" %in% input$species) & (is.data.frame(geneE[["chicken"]])) & (any(input$species > 0)))
      toggle(id = "chicken", condition = "Chicken" %in% input$species)
    }
  })

  # render title for the gene, or write error message if the gene is invalid
  output$gene_title <- renderUI({ 
    if (any(geneTitle() == 0)) {
      h4(geneT())
    }
    else {
      tagList(p("Gene"),
              h3(geneT()))
    }
    
  })
 
  # render the title for each species' plot
  output$human <- renderUI({ 
    if (any(geneTitle() == 0)) { return()}
    else if (!is.data.frame(geneEns()[["human"]])) {
      if (geneTitle()["human","id"] != "CHARACTER(0)") {
        link <- paste("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", geneTitle()["human","id"], sep = "")
        tagList(h4("Human - no data available"),
                tags$a(href=link, target="_blank", "ensembl link")
        )
      } else {
        h4("Human - no data available")
      }
    }
    else { 
      link <- paste("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", geneTitle()["human","id"], sep = "")
      tagList(h4("Human - " , geneTitle()["human","id"]),
              tags$a(href=link, target="_blank", "ensembl link")
      )
    }
  })
  
  output$mouse <- renderUI({
    if (any(geneTitle() == 0)) { return()}
    else if (!is.data.frame(geneEns()[["mouse"]])) {
      if (geneTitle()["mouse","id"] != "CHARACTER(0)") {
        link <- paste("https://www.ensembl.org/Mus_musculus/Gene/Summary?g=", geneTitle()["mouse","id"], sep = "")
        tagList(h4("Mouse - no data available"),
              tags$a(href=link, target="_blank", "ensembl link")
        )
      } else {
        h4("Mouse - no data available")
      }
    }
    else { 
      link <- paste("https://www.ensembl.org/Mus_musculus/Gene/Summary?g=", geneTitle()["mouse","id"], sep = "")
      tagList(h4("Mouse - " , geneTitle()["mouse","id"]),
              tags$a(href=link, target="_blank", "ensembl link")
      )
    }
  })
  
  output$rat <- renderUI({ 
    if (any(geneTitle() == 0)) { return()}
    else if (!is.data.frame(geneEns()[["rat"]])) {
      if (geneTitle()["rat","id"] != "CHARACTER(0)") {
        link <- paste("https://www.ensembl.org/Rattus_norvegicus/Gene/Summary?g=", geneTitle()["rat","id"], sep = "")
        tagList(h4("Rat - no data available"),
                tags$a(href=link, target="_blank", "ensembl link")
        )
      } else {
        h4("Rat - no data available")
      }
    }
    else { 
      link <- paste("https://www.ensembl.org/Rattus_norvegicus/Gene/Summary?g=", geneTitle()["rat","id"], sep = "")
      tagList(h4("Rat - " , geneTitle()["rat","id"]),
              tags$a(href=link, target="_blank", "ensembl link")
      )
    }
  })
  
  output$opossum <- renderUI({ 
    if (any(geneTitle() == 0)) { return()}
    else if (!is.data.frame(geneEns()[["opossum"]])) {
      if (geneTitle()["opossum","id"] != "CHARACTER(0)") {
        link <- paste("https://www.ensembl.org/Monodelphis_domestica/Gene/Summary?g=", geneTitle()["opossum","id"], sep = "")
        tagList(h4("Opossum - no data available"),
                tags$a(href=link, target="_blank", "ensembl link")
        )
      } else {
        h4("Opossum - no data available")
      }
    }
    else { 
      link <- paste("https://www.ensembl.org/Monodelphis_domestica/Gene/Summary?g=", geneTitle()["opossum","id"], sep = "")
      tagList(h4("Opossum - " , geneTitle()["opossum","id"]),
              tags$a(href=link, target="_blank", "ensembl link")
      )
    }
  })
  
  output$rabbit <- renderUI({ 
    if (any(geneTitle() == 0)) { return()}
    else if (!is.data.frame(geneEns()[["rabbit"]])) {
      if (geneTitle()["rabbit","id"] != "CHARACTER(0)") {
        link <- paste("https://www.ensembl.org/Oryctolagus_cuniculus/Gene/Summary?g=", geneTitle()["rabbit","id"], sep = "")
        tagList(h4("Rabbit - no data available"),
                tags$a(href=link, target="_blank", "ensembl link")
        )
      } else {
        h4("Rabbit - no data available")
      }
    }
    else { 
      link <- paste("https://www.ensembl.org/Oryctolagus_cuniculus/Gene/Summary?g=", geneTitle()["rabbit","id"], sep = "")
      tagList(h4("Rabbit - " , geneTitle()["rabbit","id"]),
              tags$a(href=link, target="_blank", "ensembl link")
      )
    }
  })

  output$macaque <- renderUI({ 
    if (any(geneTitle() == 0)) { return()}
    else if (!is.data.frame(geneEns()[["macaque"]])) {
      if (geneTitle()["macaque","id"] != "CHARACTER(0)") {
        link <- paste("https://www.ensembl.org/Macaca_mulatta/Gene/Summary?g=", geneTitle()["macaque","id"], sep = "")
        tagList(h4("Macaque - no data available"),
                tags$a(href=link, target="_blank", "ensembl link")
        )
      } else {
        h4("Macaque - no data available")
      }
    }
    else { 
      link <- paste("https://www.ensembl.org/Macaca_mulatta/Gene/Summary?g=", geneTitle()["macaque","id"], sep = "")
      tagList(h4("Macaque - " , geneTitle()["macaque","id"]),
              tags$a(href=link, target="_blank", "ensembl link")
      )
    }
  })

  output$chicken <- renderUI({
    if (any(geneTitle() == 0)) { return()}
    else if (!is.data.frame(geneEns()[["chicken"]])) {
      if (geneTitle()["chicken","id"] != "CHARACTER(0)") {
        link <- paste("https://www.ensembl.org/Gallus_gallus/Gene/Summary?g=", geneTitle()["chicken","id"], sep = "")
        tagList(h4("Chicken - no data available"),
                tags$a(href=link, target="_blank", "ensembl link")
        )
      } else {
        h4("Chicken - no data available")
      }
    }
    else { 
      link <- paste("https://www.ensembl.org/Gallus_gallus/Gene/Summary?g=", geneTitle()["chicken","id"], sep = "")
      tagList(h4("Chicken - " , geneTitle()["chicken","id"]),
              tags$a(href=link, target="_blank", "ensembl link")
      )
    }
  })
  
  # render the graphic for HUMAN
  output$humanPlot <- renderPlot({
    if (input$update == 0) {
      return()
    }
    isolate({
      if (any(input$tissues > 0)) {
        human <- geneTissues()[["human"]]
        # only render graphic if there's data for it
        if (is.data.frame(human) & (any(input$species > 0))) {
          # adjust order
          selectedTissues <- as.character(t(human["colId"]))
          newOrder <- humanOrder
          newOrder <- newOrder[newOrder$V1 %in% selectedTissues]
          xLabels <- as.character(newOrder$V2)
          newOrder <- as.character(newOrder$V1)
          # plot graphic
          plotFun(human, human["colId"], human[,4], colTissue, newOrder, xLabels, "newborn")
        }
      }
    })
  })
  
  # render the graphic for MOUSE
  output$mousePlot <- renderPlot({
    if (input$update == 0) {
      return()
    }
    isolate({
      if (any(input$tissues > 0)) {
        mouse <- geneTissues()[["mouse"]]
        # only render graphic if there's data for it
        if (is.data.frame(mouse) & (any(input$species > 0))) {
          # adjust order
          selectedTissues <- as.character(t(mouse["colId"]))
          newOrder <- mouseOrder
          newOrder <- newOrder[newOrder$V1 %in% selectedTissues]
          xLabels <- as.character(newOrder$V2)
          newOrder <- as.character(newOrder$V1)
          # plot graphic
          plotFun(mouse, mouse["colId"], mouse[,4], colTissue, newOrder, xLabels, "P0")
        }
      }
    })
  })
  
  # render the graphic for RAT
  output$ratPlot <- renderPlot({
    if (input$update == 0) {
      return()
    }
    isolate({
      if (any(input$tissues > 0)) {
        rat <- geneTissues()[["rat"]]
        # only render graphic if there's data for it
        if (is.data.frame(rat) & (any(input$species > 0))) {
          # adjust order
          selectedTissues <- as.character(t(rat["colId"]))
          newOrder <- ratOrder
          newOrder <- newOrder[newOrder$V1 %in% selectedTissues]
          xLabels <- as.character(newOrder$V2)
          newOrder <- as.character(newOrder$V1)
          # plot graphic
          plotFun(rat, rat["colId"], rat[,4], colTissue, newOrder, xLabels, "P0")
        }
      }
    })
  })
  
  # render the graphic for OPOSSUM
  output$opossumPlot <- renderPlot({
    if (input$update == 0) {
      return()
    }
    isolate({
      if (any(input$tissues > 0)) {
        opossum <- geneTissues()[["opossum"]]
        # only render graphic if there's data for it
        if (is.data.frame(opossum) & (any(input$species > 0))) {
          # adjust order
          selectedTissues <- as.character(t(opossum["colId"]))
          newOrder <- opossumOrder
          newOrder <- newOrder[newOrder$V1 %in% selectedTissues]
          xLabels <- as.character(newOrder$V2)
          newOrder <- as.character(newOrder$V1)
          # plot graphic
          plotFun(opossum, opossum["colId"], opossum[,4], colTissue, newOrder, xLabels, "P0")
        }
      }
    })
  })
  
  # render the graphic for RABBIT
  output$rabbitPlot <- renderPlot({
    if (input$update == 0) {
      return()
    }
    isolate({
      if (any(input$tissues > 0)) {
        rabbit <- geneTissues()[["rabbit"]]
        # only render graphic if there's data for it
        if (is.data.frame(rabbit) & (any(input$species > 0))) {
          # adjust order
          selectedTissues <- as.character(t(rabbit["colId"]))
          newOrder <- rabbitOrder
          newOrder <- newOrder[newOrder$V1 %in% selectedTissues]
          xLabels <- as.character(newOrder$V2)
          newOrder <- as.character(newOrder$V1)
          # plot graphic
          plotFun(rabbit, rabbit["colId"], rabbit[,4], colTissue, newOrder, xLabels, "P0")
        }
      }
    })
  })
  
  # render the graphic for MACAQUE
  output$macaquePlot <- renderPlot({
    if (input$update == 0) {
      return()
    }
    isolate({
      if (any(input$tissues > 0)) {
        macaque <- geneTissues()[["macaque"]]
        # only render graphic if there's data for it
        if (is.data.frame(macaque) & (any(input$species > 0))) {
          # adjust order
          selectedTissues <- as.character(t(macaque["colId"]))
          newOrder <- macaqueOrder
          newOrder <- newOrder[newOrder$V1 %in% selectedTissues]
          xLabels <- as.character(newOrder$V2)
          newOrder <- as.character(newOrder$V1)
          # plot graphic
          plotFun(macaque, macaque["colId"], macaque[,4], colTissue, newOrder, xLabels, "P0")
        }
      }
    })
  })
  
  # render the graphic for CHICKEN
  output$chickenPlot <- renderPlot({
    if (input$update == 0) {
      return()
    }
    isolate({
      if (any(input$tissues > 0)) {
        chicken <- geneTissues()[["chicken"]]
        # only render graphic if there's data for it
        if (is.data.frame(chicken) & (any(input$species > 0))) {
          # adjust order
          selectedTissues <- as.character(t(chicken["colId"]))
          newOrder <- chickenOrder
          newOrder <- newOrder[newOrder$V1 %in% selectedTissues]
          xLabels <- as.character(newOrder$V2)
          newOrder <- as.character(newOrder$V1)
          # plot graphic
          plotFun(chicken, chicken["colId"], chicken[,4], colTissue, newOrder, xLabels, "P0")
        }
      }
    })
  })

  # downloadHandler contains 2 arguments as functions, namely filename, content
  output$report <- downloadHandler(
      # For PDF output - "report.pdf"
      filename = function() {
        paste(geneTitle()["human","display_name"], ".pdf", sep = "")
      },
      content = function(file) {
        # Copy the report file to a temporary directory before processing it, in
        # case we don't have write permissions to the current working dir (which
        # can happen when deployed).
        tempReport <- file.path(tempdir(), "report.Rmd")
        file.copy("report.Rmd", tempReport, overwrite = TRUE)
        
        # Set up parameters to pass to Rmd document
        params <- list(gene = geneTitle(), 
                       tissues = input$tissues, 
                       species = input$species, 
                       plotData = geneTissues(),
                       humanOrder = humanOrder,
                       mouseOrder = mouseOrder,
                       ratOrder = ratOrder,
                       opossumOrder = opossumOrder,
                       rabbitOrder = rabbitOrder,
                       macaqueOrder = macaqueOrder,
                       chickenOrder = chickenOrder
                      )
        
        # Knit the document, passing in the `params` list, and eval it in a
        # child of the global environment (this isolates the code in the document
        # from the code in this app).
        rmarkdown::render(tempReport, output_file = file,
                          params = params,
                          envir = new.env(parent = globalenv())
        )
      }
    )
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste(input$speciesData, "_cpm.txt", sep = "")
    },
    content = function(file) {
      if (input$speciesData == "Human") {
        file.copy("bigdata/cpm/HumanCountsMajorTissuesCor90.Norm.txt", file)
      }
      if (input$speciesData == "Macaque") {
        file.copy("bigdata/cpm/RhesusCountsMajorTissuesCor90.Norm.txt", file)
      }
      if (input$speciesData == "Mouse") {
        file.copy("bigdata/cpm/MouseCountsMajorTissuesCor90.Norm.txt", file)
      }
      if (input$speciesData == "Rat") {
        file.copy("bigdata/cpm/RatCountsMajorTissuesCor90.Norm.txt", file)
      }
      if (input$speciesData == "Rabbit") {
        file.copy("bigdata/cpm/RabbitCountsMajorTissuesCor90.Norm.txt", file)
      }
      if (input$speciesData == "Opossum") {
        file.copy("bigdata/cpm/OpossumCountsMajorTissuesCor90.Norm.txt", file)
      }
      if (input$speciesData == "Chicken") {
        file.copy("bigdata/cpm/ChickenCountsMajorTissuesCor90.Norm.txt", file)
      }
    }
  )
}

# Run the application 
shinyApp(ui = ui, server = server)

