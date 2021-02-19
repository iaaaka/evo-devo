# Functions used in the app


# Function for summary of data
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm),
                     min   = min     (xx[[col]], na.rm=na.rm),
                     max   = max     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  datac <- plyr::rename(datac, c("mean" = measurevar))
  
  return(datac)
}

# Function that gets ensembl from gene name
# Returns: a vector (1 - id; 2 - display_name in caps), if valid; or 0, if invalid
getEnsembl <- function(geneName = NULL, species) {
  
  info <- c()
  geneName <- toupper(geneName)
  geneId <- ""
  
  for(spe in species) {
    if (count(spe[V1 == geneName]) > 0) {
      geneId <- geneName
      geneName <- spe[V1 == geneName, V2][1]
    }
  }
  
  # for(spe in species) {
  #   if (count(spe[V2 == geneName]) > 1) {
  #     geneLine <- spe[V2 == geneName][1]
  #   } else if (count(spe[V1 == geneName]) > 0) {
  #     geneLine <- spe[V1 == geneName]
  #   } else {
  #     geneLine <- spe[V2 == geneName]
  #   }
  #   info <- rbind(info, as.character(geneLine))
  # }

  for(spe in species) {
    if (count(spe[V1 == geneId]) > 0) {
      geneLine <- spe[V1 == geneId]
    } else if (count(spe[V2 == geneName]) > 1) {
      geneLine <- spe[V2 == geneName][1]
    } else {
      geneLine <- spe[V2 == geneName]
    }
    info <- rbind(info, as.character(geneLine))
  }
  
  a <- "character(0)"
  
  if (info[1] != a | info[2] != a | info[3] != a | info[4] != a | info[5] != a | info[6] != a | info[7] != a) {
    rownames(info) <- c("human", "macaque", "mouse", "rat", "rabbit", "opossum", "chicken")
    colnames(info) <- c("id", "display_name")
    info <- apply(info, 2, toupper)
    info <- as.data.frame(info, stringsAsFactors=FALSE)
    info
  } 
  else {
    return(0)
  }
}

# Function to plot graphic

plotFun <- function(speciesDF, colId, id, colourPallet, speciesOrder, xlabels, a) {
  
  ggplot(speciesDF, aes(x=colId, y=id, colour=Tissue, group=Tissue)) +
    geom_line(size=1.0) + geom_linerange(aes(ymin=min, ymax=max), alpha=.8) +
    geom_point() + theme_classic() + theme(axis.text.x = element_text(angle=90, vjust=0.5, size=8)) +
    scale_x_discrete(name="", labels = xlabels, limits=speciesOrder) + 
    geom_hline(aes(yintercept=1), colour="grey", linetype="dashed") +
    geom_vline(xintercept=which(xlabels == a), colour="grey", linetype="dashed") +
    ylab("CPM") + scale_colour_manual(values=colourPallet) + 
    theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank(), legend.position="right", axis.title.x = element_blank())
}

# Function that gets information for a gene of a certain species and summarizes data
# Args: gene id, species data.table, species' labels data.table
geneInfo <- function(id, speciesDT, speciesLabel) {
  # subset gene from table
  gene <- speciesDT[V1 == "Names" | V1 == id]
  
  # if there is data for the gene, return that data.table
  if (length(gene$V1) == 2) {
    # transpose values of gene
    gene <- transpose(gene)
    
    # create a data.table for labels with a column to rank them
    label <- speciesLabel
    label[1,1] <- "Names"
    label$V5 <- c(1:length(label$V1))
    
    # merge the two data frames by the row V1 (that holds the names of the tissues)
    gene <- merge(gene, label, by = "V1")
    # order the data frame
    gene <- gene[order(gene$V5),]
    # make the first row the colnames and erase that first row
    gene <- setDF(gene)
    colnames(gene) <- t(gene[1,])
    gene <- gene[-1,]
    # make gene column a numeric one
    gene[, 2] <- sapply(gene[, 2], as.numeric)
    
     # then summarize the values for that gene
     gene <- summarySE(gene, measurevar=id, groupvars=c("Time","Tissue"), na.rm=FALSE)
     gene$colId <- do.call(paste, c(gene[c("Tissue", "Time")], sep = "."))
     gene
   }
   # if there is no data for the gene, return number 4
   else {
     return(4)
   }
}

# Selection of tissues

tissues <- function(speciesList, inputTissues) {
  newList <- list()
  counter <- 1
  for(species in speciesList) {
    if (is.data.frame(species)){
      newList[[counter]] <- species[species$Tissue %in% inputTissues,]
    }
    else {
      newList[[counter]] <- species[1]
    }
    counter <- counter+1
  }
  names(newList) <- names(speciesList)
  return(newList)
}




# # Function that gets ensembl from gene name
# # Returns: a vector (1 - id; 2 - display_name in caps), if valid; or 0, if invalid
# library(httr)
# library(jsonlite)
# library(xml2)
# 
# getEnsembl <- function(geneName = NULL) {
#   
#   info <- c()
#   species <- c("human", "macaque", "mouse", "rat", "rabbit", "opossum", "chicken")
#   ##  ext <- paste0("/xrefs/symbol/", spe, "/", geneName, "?")
#   server <- "http://rest.ensembl.org"
#   
#   for(spe in species) {
#     ext <- paste0("/lookup/symbol/", spe, "/", geneName, "?")
#     
#     r <- GET(paste(server, ext, sep = ""), content_type("application/json"))
#     
#     #stop_for_status(r)
#     message_for_status(r, paste(geneName, ' - ', spe))
#     info <- rbind(info, as.character(fromJSON(toJSON(content(r)))[c("id", "display_name")]))
#   }
#   
#   if (info[1] != "NULL" | info[2] != "NULL" | info[3] != "NULL" | info[4] != "NULL" | info[5] != "NULL" | info[6] != "NULL" | info[7] != "NULL") {
#     rownames(info) <- species
#     colnames(info) <- c("id", "display_name")
#     info <- apply(info, 2, toupper)
#     info
#   }
#   else {
#     return(0)
#   }
# }
