#
# This is the Evo-Devo shiny web application. 
#

# Load libraries
library(shiny)
library(shinyjs)
#library(plyr)
#library(dplyr)
#library(data.table)
print(list.files())
# # Load functions
# source("helpers.R")
# 
# # Load file with codes
# human <- fread("data/listGenes/listHumanGenes.csv", header = FALSE)
# macaque <- fread("data/listGenes/listMacaqueGenes.csv", header = FALSE)
# mouse <- fread("data/listGenes/listMouseGenes.csv", header = FALSE)
# rat <- fread("data/listGenes/listRatGenes.csv", header = FALSE)
# rabbit <- fread("data/listGenes/listRabbitGenes.csv", header = FALSE)
# opossum <- fread("data/listGenes/listOpossumGenes.csv", header = FALSE)
# chicken <- fread("data/listGenes/listChickenGenes.csv", header = FALSE)
# 
# extraCols <- c("V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14")
# human[, extraCols] <- ""
# macaque[, extraCols] <- ""
# mouse[, extraCols] <- ""
# rat[, extraCols] <- ""
# rabbit[, extraCols] <- ""
# opossum[, extraCols] <- ""
# chicken[, extraCols] <- ""
# 
# listGenes <- list(human, macaque, mouse, rat, rabbit, opossum, chicken)
# 
# # new genes
# # Load file with codes
# human <- fread("data/newListGenes/human.final.txt", header = FALSE)
# macaque <- fread("data/newListGenes/rhesus.final.txt", header = FALSE)
# mouse <- fread("data/newListGenes/mouse.final.txt", header = FALSE)
# rat <- fread("data/newListGenes/rat.final.txt", header = FALSE)
# rabbit <- fread("data/newListGenes/rabbit.final.txt", header = FALSE)
# opossum <- fread("data/newListGenes/opossum.final.txt", header = FALSE)
# chicken <- fread("data/newListGenes/chicken.final.txt", header = FALSE)
# 
# newListGenes <- list(human, macaque, mouse, rat, rabbit, opossum, chicken)
# 
# # Load gene aliases file as data frame
# aliasDF <- read.csv("data/aliasDF.csv", header = FALSE, strip.white = TRUE, stringsAsFactors = FALSE)
# 
# ## Load data - for human
# humanLabels <- fread("data/humanLabels.csv", header = FALSE)
# # cpm
# # humanGenes <- fread("data/humanRounded_cpm.csv", header = FALSE)
# # rpkm
# humanGenes <- fread("data/humanRounded_rpkm.csv", header = FALSE)
# humanOrder <- fread("data/humanOrder.csv", header = FALSE)
# 
# ## Load data - for mouse
# mouseLabels <- fread("data/mouseLabels.csv", header = FALSE)
# # cpm
# # mouseGenes <- fread("data/mouseRounded_cpm.csv", header = FALSE)
# # rpkm
# mouseGenes <- fread("data/mouseRounded_rpkm.csv", header = FALSE)
# mouseOrder <- fread("data/mouseOrder.csv", header = FALSE)
# 
# ## Load data - for rat
# ratLabels <- fread("data/ratLabels.csv", header = FALSE)
# # cpm
# # ratGenes <- fread("data/ratRounded_cpm.csv", header = FALSE)
# # rpkm
# ratGenes <- fread("data/ratRounded_rpkm.csv", header = FALSE)
# ratOrder <- fread("data/ratOrder.csv", header = FALSE)
# 
# ## Load data - for opossum
# opossumLabels <- fread("data/opossumLabels.csv", header = FALSE)
# # cpm
# # opossumGenes <- fread("data/opossumRounded_cpm.csv", header = FALSE)
# # rpkm
# opossumGenes <- fread("data/opossumRounded_rpkm.csv", header = FALSE)
# opossumOrder <- fread("data/opossumOrder.csv", header = FALSE)
# 
# ## Load data - for rabbit
# rabbitLabels <- fread("data/rabbitLabels.csv", header = FALSE)
# # cpm
# # rabbitGenes <- fread("data/rabbitRounded_cpm.csv", header = FALSE)
# # rpkm
# rabbitGenes <- fread("data/rabbitRounded_rpkm.csv", header = FALSE)
# rabbitOrder <- fread("data/rabbitOrder.csv", header = FALSE)
# 
# ## Load data - for macaque
# macaqueLabels <- fread("data/macaqueLabels.csv", header = FALSE)
# # cpm
# # macaqueGenes <- fread("data/macaqueRounded_cpm.csv", header = FALSE)
# # rpkm
# macaqueGenes <- fread("data/macaqueRounded_rpkm.csv", header = FALSE)
# macaqueOrder <- fread("data/macaqueOrder.csv", header = FALSE)
# 
# ## Load data - for chicken
# chickenLabels <- fread("data/chickenLabels.csv", header = FALSE)
# # cpm
# # chickenGenes <- fread("data/chickenRounded_cpm.csv", header = FALSE)
# # rpkm
# chickenGenes <- fread("data/chickenRounded_rpkm.csv", header = FALSE)
# chickenOrder <- fread("data/chickenOrder.csv", header = FALSE)
# 
# # Define pallet
colTissue <- c(Brain = "#3399CC", Cerebellum = "#33CCFF", Heart = "#CC0000", Kidney = "#CC9900", Liver = "#339900", Ovary = "#CC3399", Testis = "#FF6600")

# Vector of species
species <- c("Human", "Macaque", "Mouse", "Rat", "Rabbit", "Opossum", "Chicken")

# Vector of species initials
speciesInitials <- c("ENSG000", "ENSMMUG", "ENSMUSG", "ENSRNOG", "ENSOCUG", "ENSMODG", "ENSGALG")
names(speciesInitials) <- species

# Chosen method for search - gene name or ensembl id
by <- "gene name"

# Define UI for application
ui <- fluidPage(
	
	# link css styles
	shiny::tags$head(shiny::tags$link(rel = "stylesheet", type = "text/css", href = "styles.css")),
	# link font
	shiny::tags$head(shiny::tags$link(rel = "stylesheet", type = "text/css", href = "https://fonts.googleapis.com/css?family=Quicksand")),
	# link javascript file
	shiny::tags$head(shiny::tags$script(src = "scripts.js")),
	# favicon
	#shiny::tags$head(tags$link(rel="shortcut icon", href="favicon.ico")),

	useShinyjs(),

	absolutePanel(id = "absNav", top = 0, left= 0, right = 0, fixed = TRUE,
		shiny::tags$div(
			id = "topBar",
			shiny::tags$img(id = "mainPicture", src = "Mammalian_evolutionTransparent.png"),
			fluidRow(id = "title",
			         column(12, titlePanel("Evo-devo mammalian organs")),
			         column(12, 
			                shiny::tags$a(id = "link_hklab", href = "https://www.zmbh.uni-heidelberg.de/kaessmann", target="_blank",
			                              "Kaessmann Lab"),
			                shiny::tags$a(id = "link_lnc", href = "https://apps.kaessmannlab.org/lncRNA_app/", target="_blank", 
			                              "lncRNA app")
			         )
			),
			htmlOutput("gene_title"),
			htmlOutput("aliases"),
      shiny::tags$div(class = "pbusy",
                      shiny::tags$img(id = "spinner", src="ajax-loader.gif", width = "80px"),
                      p("Loading...")
      ),
			shiny::tags$a(id = "scroll", href = "#top", shiny::tags$img(id = "arrow", src = "arrow.png"), "top")
		)),

	#info text
  	shiny::tags$div(
    	class = "info-text",
		HTML(
		'<p>This Shiny app allows the interactive exploration of individual gene expression profiles across organs, developmental stages and species.
		<br>
			If you are using this app, please consult and cite our publication:<br>
			<span>Cardoso-Moreira et al. Gene expression across mammalian organ development. Nature (2019)</span><br>
			<a href="https://www.nature.com/articles/s41586-019-1338-5" target="_blank" rel="noopener">https://www.nature.com/articles/s41586-019-1338-5</a></p>'
		)
	),
	
	sidebarLayout(
		sidebarPanel(width = 4, class = "sideChoices",
			shiny::tags$div(class = "side_panel",
				# Panel for the search options
				# for specifying a gene and selecting it
				fluidRow(
					fluidRow(
						column(12, h4("Choose a gene:"))
					),
					# "gene" - create inputs of type textInput() and selectInput() for gene name search
					fluidRow(
						column(6, textInput("gene", label = h5("Gene name"), value = NULL, placeholder = "type gene")),
						column(6, selectInput("speciesI", label = h5("Species"), choices = speciesInitials, 
												selected = "Human"))
					),
					fluidRow(
						column(12, h4("or:"))
					),
					# "id" - create input of type textInput() for ensembl id search
					textInput("id", label = h5("Ensembl ID"), value = NULL, placeholder = "type ensembl id"),
					fluidRow(column(12, shiny::tags$hr())),
					fluidRow(
						# "by" - create input that chooses which method to use
						#column(1, h5("By")),
						shinyjs::hidden(column(6, selectInput("by", label = NULL, choices = c("ensembl id", "gene name")), 
											selected = "ensembl id")),
						# "click" - Create submittButton
						column(12, actionButton("click", label = "Search"), style = "text-align: center")
					)
				),
				# Panel for filtering options
				fluidRow(
					fluidRow(
						column(12, h4("Show me:"))
					),
					# "species" - create input of type checkboxGroupInput() for selection of species,
					column(6, checkboxGroupInput("species", label = h5("Species"),
													choices = species, 
													selected = species)),
					# "tissues" - create input of type checkboxGroupInput() for selection of tissues
					column(6, checkboxGroupInput("tissues", label = h5("Organ"),
													choices = names(colTissue),
													selected = names(colTissue))),
					fluidRow(column(12, shiny::tags$hr())),
					column(12, actionButton("update", label = "Update"), style = "text-align: center")
				),
				# Panel for downloading options
				fluidRow(
					fluidRow(
						column(12, h4("Downloads:"))
					),
					# downloading plots
					column(8, downloadButton(outputId = "report", label = "Download plots"), style = "margin-top: 20px"),
					fluidRow(column(12, shiny::tags$hr())),
					# downloading expression tables
					column(12, p("Download expression tables (RPKM)")),
					column(9, selectInput(inputId = "speciesData", label = NULL, choices = species)),
					column(3, downloadButton(outputId = "downloadData", label = ""))
				)
			)
		),
		
		# Panel for the plots
		mainPanel(width = 8, class = "sideChoices",
			shiny::tags$div(id = "to_hide",
				fluidRow(
					shiny::tags$div(uiOutput("human"),
									plotOutput("humanPlot"),
									htmlOutput("hr_human")),
					shiny::tags$div(uiOutput("macaque"),
									plotOutput("macaquePlot"),
									htmlOutput("hr_macaque")),
					shiny::tags$div(uiOutput("mouse"),
									plotOutput("mousePlot"),
									htmlOutput("hr_mouse")),
					shiny::tags$div(uiOutput("rat"),
									plotOutput("ratPlot"),
									htmlOutput("hr_rat")),
					shiny::tags$div(uiOutput("rabbit"),
									plotOutput("rabbitPlot"),
									htmlOutput("hr_rabbit")),
					shiny::tags$div(uiOutput("opossum"),
									plotOutput("opossumPlot"),
									htmlOutput("hr_opossum")),
					shiny::tags$div(uiOutput("chicken"),
									plotOutput("chickenPlot"),
									htmlOutput("hr_chicken"))
				)
			)
		)
	),
	shiny::tags$footer(
		HTML(
		'<div><p>Contact: <a href="mailto:m.moreira@zmbh.uni-heidelberg.de">m.moreira@zmbh.uni-heidelberg.de</a> | <a href="http://privacy.kaessmannlab.org/evodevo.html">Privacy Notice</a></p></div>'
		)
	)
)

# Define server logic required to draw the plots
server <- function(input, output, session) {

	observeEvent(input$gene, {
		updateSelectInput(session, "by", selected = "gene name")
	})

	observeEvent(input$id, {
		updateSelectInput(session, "by", selected = "ensembl id")
	})
	
	# Get information (ids and gene names) for each species
	## create another input in ui (input$gene and input$id)
	geneTitle <- eventReactive(input$click, {
		# if using ensembl id, use this funtion
		if (input$by == "ensembl id") {
			geneT <- getEnsembl(input$id, newListGenes, listGenes)
		}
		# if searching gene names, use this function
		else if (input$by == "gene name") {
			geneT <- getEnsembl(input$gene, newListGenes, listGenes, input$speciesI)
		}
		shinyjs::reset(id = "id")
		return(geneT)
	})
	
	# Create gene title string
	geneT <- eventReactive(input$click, {
		# show title and aliases elements
		if (any(geneTitle() == 0)) {
			if (input$by == "ensembl id") {
				return("Sorry, we cannot find that gene name.")
			} else {
				geneA <- trueGeneName(input$gene, aliasDF)
				if (any(geneA != 0)) {
					if (any(length(geneA) > 1)) {
						if (any(geneA[1] == "0")) {
							geneA <- geneA[-1]
							geneA <- paste(paste(geneA[-length(geneA)], collapse = ', '), "or", geneA[length(geneA)])
						} else {
							geneA <- geneA[1]
						}
						if (toupper(geneA[1]) == toupper(input$gene)) {
							return("Sorry, we cannot find that gene name for the selected species.")
						} else {
							return(paste("Sorry, we cannot find that gene name. \nTry", geneA, "instead."))
						}
					}
				} else {
					return("Sorry, we cannot find that gene name.")
				}
			}
		}
		else if (input$by == "ensembl id") {
			if (any(rightSpecies(input$id, newListGenes) != 0)) {
				x <- as.character(rightSpecies(input$id, newListGenes)[1,2])
			}
			else if (any(rightSpecies(input$id, listGenes) != 0)) {
				x <- as.character(rightSpecies(input$id, listGenes)[1,2])
			}
			initialSpecies <- substr(input$id, 1, 7)
			return(x)
		}
		else if (input$by == "gene name") {
			if (any(rightSpecies(input$gene, newListGenes, input$speciesI) != 0)) {
				x <- as.character(rightSpecies(input$gene, newListGenes, input$speciesI)[1,2])
			}
			else if (any(rightSpecies(input$gene, listGenes, input$speciesI) != 0)) {
				x <- as.character(rightSpecies(input$gene, listGenes, input$speciesI)[1,2])
			}
			initialSpecies <- input$speciesI
			return(x)
		}
		else {
			return("error")
		}
	})

	# 
	initialSpecies <- eventReactive(input$click, {
		if (input$by == "ensembl id") {
			iniSpe <- substr(input$id, 1, 7)
		}
		else if (input$by == "gene name") {
			iniSpe <- input$speciesI
		}
		iniSpe <- names(speciesInitials[which(grepl(iniSpe, speciesInitials))])	
		return(iniSpe)
	})
	
	# Check for synonyms
	geneAliases <- eventReactive(input$click, {
		synonyms <- trueGeneName(geneT(), aliasDF)
		if (any(length(synonyms) > 1)) {
			return(paste("Synonyms: ", paste(synonyms[], collapse = ', ')))
			#return(paste("Synonyms: ", paste(synonyms[-1], collapse = ', ')))
		} else { return(0) }
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
			#toggle(id = "hr_human", condition = ("Human" %in% input$species) & (is.data.frame(geneE[["human"]])) & (any(input$species > 0)))
			toggle(id = "human", condition = "Human" %in% input$species)
			#toggle(id = "hr_human", condition = "Human" %in% input$species)
	
			# if MOUSE is selected
			toggle(id = "mousePlot", condition = ("Mouse" %in% input$species) & (is.data.frame(geneE[["mouse"]])) & (any(input$species > 0)))
			#toggle(id = "hr_mouse", condition = ("Mouse" %in% input$species) & (is.data.frame(geneE[["mouse"]])) & (any(input$species > 0)))
			toggle(id = "mouse", condition = "Mouse" %in% input$species)
			#toggle(id = "hr_mouse", condition = "Mouse" %in% input$species)
	
			# if RAT is selected
			toggle(id = "ratPlot", condition = ("Rat" %in% input$species) & (is.data.frame(geneE[["rat"]])) & (any(input$species > 0)))
			#toggle(id = "hr_rat", condition = ("Rat" %in% input$species) & (is.data.frame(geneE[["rat"]])) & (any(input$species > 0)))
			toggle(id = "rat", condition = "Rat" %in% input$species)
			#toggle(id = "hr_rat", condition = "Rat" %in% input$species)
	
			# if OPOSSUM is selected
			toggle(id = "opossumPlot", condition = ("Opossum" %in% input$species) & (is.data.frame(geneE[["opossum"]])) & (any(input$species > 0)))
			#toggle(id = "hr_opossum", condition = ("Opossum" %in% input$species) & (is.data.frame(geneE[["opossum"]])) & (any(input$species > 0)))
			toggle(id = "opossum", condition = "Opossum" %in% input$species)
			#toggle(id = "hr_opossum", condition = "Opossum" %in% input$species)
	
			# if RABBIT is selected
			toggle(id = "rabbitPlot", condition = ("Rabbit" %in% input$species) & (is.data.frame(geneE[["rabbit"]])) & (any(input$species > 0)))
			#toggle(id = "hr_rabbit", condition = ("Rabbit" %in% input$species) & (is.data.frame(geneE[["rabbit"]])) & (any(input$species > 0)))
			toggle(id = "rabbit", condition = "Rabbit" %in% input$species)
			#toggle(id = "hr_rabbit", condition = "Rabbit" %in% input$species)
	
			# if MACAQUE is selected
			toggle(id = "macaquePlot", condition = ("Macaque" %in% input$species) & (is.data.frame(geneE[["macaque"]])) & (any(input$species > 0)))
			#toggle(id = "hr_macaque", condition = ("Macaque" %in% input$species) & (is.data.frame(geneE[["macaque"]])) & (any(input$species > 0)))
			toggle(id = "macaque", condition = "Macaque" %in% input$species)
			#toggle(id = "hr_macaque", condition = "Macaque" %in% input$species)
	
			# if CHICKEN is selected
			toggle(id = "chickenPlot", condition = ("Chicken" %in% input$species) & (is.data.frame(geneE[["chicken"]])) & (any(input$species > 0)))
			#toggle(id = "hr_chicken", condition = ("Chicken" %in% input$species) & (is.data.frame(geneE[["chicken"]])) & (any(input$species > 0)))
			toggle(id = "chicken", condition = "Chicken" %in% input$species)
			#toggle(id = "hr_chicken", condition = "Chicken" %in% input$species)
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

	output$aliases <- renderUI({
		if (any(geneAliases() != 0)) {
			p(geneAliases())
		}
	})

	# render the title for each species' plot
	output$human <- renderUI({ 
		if (any(geneTitle() == 0)) { return() }
		else if (!is.data.frame(geneEns()[["human"]])) {
			if (geneTitle()["human","id"] != "") {
				link <- paste("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", geneTitle()["human","id"], sep = "")
				tagList(h4("Human - ", geneTitle()["human","id"], " (no data available)"),
						shiny::tags$a(href=link, target="_blank", "ensembl link")
				)
			} else {
				#h4("Human - no correspondence found with searched species")
				h4("Human - no data available")
			}
		}
		else { 
			link <- paste("https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=", geneTitle()["human","id"], sep = "")
			tagList(h4("Human - " , geneTitle()["human","id"]),
					shiny::tags$a(href=link, target="_blank", "ensembl link")
			)
		}
	})
	
	output$mouse <- renderUI({
		if (any(geneTitle() == 0)) { return()}
		else if (!is.data.frame(geneEns()[["mouse"]])) {
			if (geneTitle()["mouse","id"] != "") {
				link <- paste("https://www.ensembl.org/Mus_musculus/Gene/Summary?g=", geneTitle()["mouse","id"], sep = "")
				tagList(h4("Mouse - ", geneTitle()["mouse","id"], " (no data available)"),
						shiny::tags$a(href=link, target="_blank", "ensembl link")
				)
			} else {
				h4("Mouse - no data available")
			}
		}
		else { 
			link <- paste("https://www.ensembl.org/Mus_musculus/Gene/Summary?g=", geneTitle()["mouse","id"], sep = "")
			tagList(h4("Mouse - " , geneTitle()["mouse","id"]),
					shiny::tags$a(href=link, target="_blank", "ensembl link")
			)
		}
	})
	
	output$rat <- renderUI({ 
		if (any(geneTitle() == 0)) { return()}
		else if (!is.data.frame(geneEns()[["rat"]])) {
			if (geneTitle()["rat","id"] != "") {
				link <- paste("https://www.ensembl.org/Rattus_norvegicus/Gene/Summary?g=", geneTitle()["rat","id"], sep = "")
				tagList(h4("Rat - ", geneTitle()["rat","id"], " (no data available)"),
						shiny::tags$a(href=link, target="_blank", "ensembl link")
				)
			} else {
				h4("Rat - no data available")
			}
		}
		else { 
			link <- paste("https://www.ensembl.org/Rattus_norvegicus/Gene/Summary?g=", geneTitle()["rat","id"], sep = "")
			tagList(h4("Rat - " , geneTitle()["rat","id"]),
					shiny::tags$a(href=link, target="_blank", "ensembl link")
			)
		}
	})
	
	output$opossum <- renderUI({ 
		if (any(geneTitle() == 0)) { return()}
		else if (!is.data.frame(geneEns()[["opossum"]])) {
			if (geneTitle()["opossum","id"] != "") {
				link <- paste("https://www.ensembl.org/Monodelphis_domestica/Gene/Summary?g=", geneTitle()["opossum","id"], sep = "")
				tagList(h4("Opossum - ", geneTitle()["opossum","id"], " (no data available)"),
						shiny::tags$a(href=link, target="_blank", "ensembl link")
				)
			} else {
				h4("Opossum - no data available")
			}
		}
		else { 
			link <- paste("https://www.ensembl.org/Monodelphis_domestica/Gene/Summary?g=", geneTitle()["opossum","id"], sep = "")
			tagList(h4("Opossum - " , geneTitle()["opossum","id"]),
					shiny::tags$a(href=link, target="_blank", "ensembl link")
			)
		}
	})
	
	output$rabbit <- renderUI({ 
		if (any(geneTitle() == 0)) { return()}
		else if (!is.data.frame(geneEns()[["rabbit"]])) {
			if (geneTitle()["rabbit","id"] != "") {
				link <- paste("https://www.ensembl.org/Oryctolagus_cuniculus/Gene/Summary?g=", geneTitle()["rabbit","id"], sep = "")
				tagList(h4("Rabbit - ", geneTitle()["rabbit","id"], " (no data available)"),
						shiny::tags$a(href=link, target="_blank", "ensembl link")
				)
			} else {
				h4("Rabbit - no data available")
			}
		}
		else { 
			link <- paste("https://www.ensembl.org/Oryctolagus_cuniculus/Gene/Summary?g=", geneTitle()["rabbit","id"], sep = "")
			tagList(h4("Rabbit - " , geneTitle()["rabbit","id"]),
					shiny::tags$a(href=link, target="_blank", "ensembl link")
			)
		}
	})

	output$macaque <- renderUI({ 
		if (any(geneTitle() == 0)) { return()}
		else if (!is.data.frame(geneEns()[["macaque"]])) {
			if (geneTitle()["macaque","id"] != "") {
				link <- paste("https://www.ensembl.org/Macaca_mulatta/Gene/Summary?g=", geneTitle()["macaque","id"], sep = "")
				tagList(h4("Macaque - ", geneTitle()["macaque","id"], " (no data available)"),
						shiny::tags$a(href=link, target="_blank", "ensembl link")
				)
			} else {
				h4("Macaque - no data available")
			}
		}
		else { 
			link <- paste("https://www.ensembl.org/Macaca_mulatta/Gene/Summary?g=", geneTitle()["macaque","id"], sep = "")
			tagList(h4("Macaque - " , geneTitle()["macaque","id"]),
					shiny::tags$a(href=link, target="_blank", "ensembl link")
			)
		}
	})

	output$chicken <- renderUI({
		if (any(geneTitle() == 0)) { return()}
		else if (!is.data.frame(geneEns()[["chicken"]])) {
			if (geneTitle()["chicken","id"] != "") {
				link <- paste("https://www.ensembl.org/Gallus_gallus/Gene/Summary?g=", geneTitle()["chicken","id"], sep = "")
				tagList(h4("Chicken - ", geneTitle()["chicken","id"], " (no data available)"),
						shiny::tags$a(href=link, target="_blank", "ensembl link")
				)
			} else {
				h4("Chicken - no data available")
			}
		}
		else { 
			link <- paste("https://www.ensembl.org/Gallus_gallus/Gene/Summary?g=", geneTitle()["chicken","id"], sep = "")
			tagList(h4("Chicken - " , geneTitle()["chicken","id"]),
					shiny::tags$a(href=link, target="_blank", "ensembl link")
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

	# render hr for HUMAN
	output$hr_human <- renderUI({
	  if (input$update == 0) { return() }
	  isolate({
	    if (any(input$tissues > 0) & (any(input$species > 0))) {
	      shiny::tags$hr()
	    }
	  })
	})
	
	# render hr for MACAQUE
	output$hr_macaque <- renderUI({
	  if (input$update == 0) { return() }
	  isolate({
	    if (any(input$tissues > 0) & (any(input$species > 0))) {
	      shiny::tags$hr()
	    }
	  })
	})
	
	# render hr for MOUSE
	output$hr_mouse <- renderUI({
	  if (input$update == 0) { return() }
	  isolate({
	    if (any(input$tissues > 0) & (any(input$species > 0))) {
	      shiny::tags$hr()
	    }
	  })
	})
	
	# render hr for RAT
	output$hr_rat <- renderUI({
	  if (input$update == 0) { return() }
	  isolate({
	    if (any(input$tissues > 0) & (any(input$species > 0))) {
	      shiny::tags$hr()
	    }
	  })
	})

	# render hr for RABBIT
	output$hr_rabbit <- renderUI({
	  if (input$update == 0) { return() }
	  isolate({
	    if (any(input$tissues > 0) & (any(input$species > 0))) {
	      shiny::tags$hr()
	    }
	  })
	})
	
	# render hr for OPOSSUM
	output$hr_opossum <- renderUI({
	  if (input$update == 0) { return() }
	  isolate({
	    if (any(input$tissues > 0) & (any(input$species > 0))) {
	      shiny::tags$hr()
	    }
	  })
	})
	
	# render hr for CHICKEN
	output$hr_chicken <- renderUI({
	  if (input$update == 0) { return() }
	  isolate({
	    if (any(input$tissues > 0) & (any(input$species > 0))) {
	      shiny::tags$hr()
	    }
	  })
	})

	# downloadHandler contains 2 arguments as functions, namely filename, content
	output$report <- downloadHandler(
			# For PDF output - "report.pdf"
			filename = function() {
				paste(geneT(), ".pdf", sep = "")
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
			paste(input$speciesData, "_rpkm.txt", sep = "")
		},
		content = function(file) {
			# if (input$speciesData == "Human") {
			# 	file.copy("bigdata/cpm/HumanCountsMajorTissuesCor90.Norm.txt", file)
			# }
			# if (input$speciesData == "Macaque") {
			# 	file.copy("bigdata/cpm/RhesusCountsMajorTissuesCor90.Norm.txt", file)
			# }
			# if (input$speciesData == "Mouse") {
			# 	file.copy("bigdata/cpm/MouseCountsMajorTissuesCor90.Norm.txt", file)
			# }
			# if (input$speciesData == "Rat") {
			# 	file.copy("bigdata/cpm/RatCountsMajorTissuesCor90.Norm.txt", file)
			# }
			# if (input$speciesData == "Rabbit") {
			# 	file.copy("bigdata/cpm/RabbitCountsMajorTissuesCor90.Norm.txt", file)
			# }
			# if (input$speciesData == "Opossum") {
			# 	file.copy("bigdata/cpm/OpossumCountsMajorTissuesCor90.Norm.txt", file)
			# }
			# if (input$speciesData == "Chicken") {
			# 	file.copy("bigdata/cpm/ChickenCountsMajorTissuesCor90.Norm.txt", file)
			# }
			if (input$speciesData == "Human") {
				file.copy("bigdata/rpkm/HumanRpkmMajorTissuesCor90.Norm.txt", file)
			}
			if (input$speciesData == "Macaque") {
				file.copy("bigdata/rpkm/MacaqueRpkmMajorTissuesCor90.Norm.txt", file)
			}
			if (input$speciesData == "Mouse") {
				file.copy("bigdata/rpkm/MouseRpkmMajorTissuesCor90.Norm.txt", file)
			}
			if (input$speciesData == "Rat") {
				file.copy("bigdata/rpkm/RatRpkmMajorTissuesCor90.Norm.txt", file)
			}
			if (input$speciesData == "Rabbit") {
				file.copy("bigdata/rpkm/RabbitRpkmMajorTissuesCor90.Norm.txt", file)
			}
			if (input$speciesData == "Opossum") {
				file.copy("bigdata/rpkm/OpossumRpkmMajorTissuesCor90.Norm.txt", file)
			}
			if (input$speciesData == "Chicken") {
				file.copy("bigdata/rpkm/ChickenRpkmMajorTissuesCor90.Norm.txt", file)
			}
		}
	)
}

# Run the application 
shinyApp(ui = ui, server = server)

