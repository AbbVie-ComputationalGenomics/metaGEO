
getGEO_ABV <- function (GEO = NULL, filename = NULL, destdir = tempdir(), GSElimits = NULL,
    					GSEMatrix = TRUE, AnnotGPL = FALSE, getGPL = TRUE) {
    con <- NULL
    sprintf("GEO arg is: %s", GEO )
    if (!is.null(GSElimits)) {
        if (length(GSElimits) != 2) {
            stop("GSElimits should be an integer vector of length 2, like (1,10) to include GSMs 1 through 10")
        }
    }
    if (is.null(GEO) & is.null(filename)) {
        stop("You must supply either a filename of a GEO file or a GEO accession")
    }
    if (is.null(filename)) {
        GEO <- toupper(GEO)
        geotype <- toupper(substr(GEO, 1, 3))
        if (GSEMatrix & geotype == "GSE") {
            return(getAndParseGSEMatrices_ABV(GEO, destdir, AnnotGPL = AnnotGPL, getGPL = getGPL))
        }
        filename <- getGEOfile(GEO, destdir = destdir, AnnotGPL = AnnotGPL)
    }
    ret <- parseGEO(filename, GSElimits)
    return(ret)
}


getAndParseGSEMatrices_ABV <- function (GEO, destdir, AnnotGPL, getGPL = TRUE) {
    GEO <- toupper(GEO)
    stub = gsub("\\d{1,3}$", "nnn", GEO, perl = TRUE)
    gdsurl <- "https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/"
    b = getDirListing_ABV(sprintf(gdsurl, stub, GEO))
    message(sprintf("Found %d file(s)", length(b)))
    ret <- list()
    for (i in 1:length(b)) {
        message(b[i])
        destfile = file.path(destdir, b[i])
        if (file.exists(destfile)) {
            message(sprintf("Using locally cached version: %s", 
                destfile))
        }
        else {
            download.file(sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%s/%s/matrix/%s", 
                stub, GEO, b[i]), destfile = destfile, mode = "wb", 
                method = getOption("download.file.method.GEOquery"))
        }
        ret[[b[i]]] <- GEOquery:::parseGSEMatrix(destfile, destdir = destdir, 
            AnnotGPL = AnnotGPL, getGPL = getGPL)$eset
    }
    return(ret)
}


getDirListing_ABV <- function (url){
    message(url)
    a <- RCurl::getURL(url)
    if (grepl("HTML", a)) {
        message("# Processing HTML result page (behind a proxy?) ... ", 
            appendLF = FALSE)
	sa <- gsub("HREF", "href", a, fixed = TRUE)
        sa <- strsplit(sa, "href", fixed = TRUE)[[1]]
	sa <- sa[grep('matrix.txt.gz', sa)]
	b <- as.matrix(as.character(sapply(sa, function(a) gsub('.*GSE', 'GSE', gsub('matrix.txt.gz.*', 'matrix.txt.gz', a)))))
        message("OK")
    }
    else {
        tmpcon <- textConnection(a, "r")
        b <- read.table(tmpcon)
        close(tmpcon)
    }
    b <- as.character(b[, ncol(b)])
    return(b)
}


#' Get a gene expression set from GEO
#'
#' @param dataset A GEO gene expresseion set accession number (e.g., GSE1500)
#' @return A bioconductor gene expression set
getGSE <- function(dataset) {
	gset <- getGEO_ABV(GEO=dataset, GSEMatrix=TRUE)
	if (length(gset) > 1) idx <- grep("GPL[0-9]", attr(gset, "names")) else idx <- 1
	gset <- gset[idx]
	print(paste0('Dataset: ', dataset, ' returned: ', length(gset)))
	return(gset)
}


#' Get a range of gene expression sets
#'
#' @param datasets A character vector giving GEO
#' @return A list of bioconductor gene expression sets
getGSEs <- function(datasets) {
  mclapply(datasets, function(i) try(getGSE(i)), mc.cores=4)
}


#' Get a GEO Platform, or GPL, for a given GEO Series accession number
#' 
#' @param 'expSet' - A bioconductor ExpressionSet
#' @return A bioconductor GEO Platform entity
getGPL <- function(expSet) {
  print(expSet)
  gpl <- getGEO(expSet@annotation)
  return(gpl)
}


#' Get all GPL objects for the given datasets
#' 
#' @param expressionSets A list of GEO Expression Sets
#' @return A list of unique bioconductor GPL objects
getGPLs <- function(expressionSets) {
  library(parallel)
  gplSets <- mclapply(expressionSets, function(gse) try(getGPL(gse)), mc.cores=4)
  return(gplSets)
}


#' #' For a list of gene expression sets, return a corresponding set of matched gene annotation information
#' #' @param data A list of gene expression sets with annotation attribute giving the appropriate gene annotation file to use
#' #' @return A list of gene annotation objects of class GPL from GEOquery
#' getGeneAnnos <- function(data) {
#' 	annotations <- lapply(data, function(i) {GEOquery::getGEO(Biobase::annotation(i))})
#' 	return(annotations)
#' }


#' Extract gene expression data for a given gene
#' @param gene A gene annotation symbol
#' @param data  A list of gene expression sets
#' @return A list with two elements.  Each element is a list of the same length as data and gene_maps.  gene_data gives the gene expression data for any probe found to match the given gene.   gene_info gives the info for the probes found to match.  The match is made based on a grep of the gene symbol to any column with "gene" in the label (case insensitive).  This is very liberal and will probably pull out may probes you didn't want.   It will be important to display the matches made and encourage manual review.
getSingleGeneDataGrepGPL <- function(gene, data) {
	D=length(data)
	indices <- lapply(1:D, function(i) {
    	gene_text <- do.call(paste, lapply(data[[i]]@featureData@data[,grep('[gG][eE][nN[eE]', colnames(data[[1]]@featureData@data))], as.character))
    	grep(gene, gene_text)})
	gene_data <- lapply(1:D, function(i) data[[i]][indices[[i]],])
	#the gene info returned for some reason removes the important columns for the 4th test dataset.   Need to track down why but I think it is still returning the correct probes.
	gene_info <- lapply(1:D, function(i) data[[i]]@featureData@data[indices[[i]],])
	return(list(gene_data=gene_data, gene_info=gene_info))
}


#' Return index marking location of probes corresponding to a specified gene
#' @param gene A gene annotation symbol
#' @param data  A list of gene expression sets
#' @param gene_maps a list of GPL objects giving matching gene annotation information for the expression sets given in data
#' @return A list with two elements.  Each element is a list of the same length as data and gene_maps.  gene_data gives the gene expression data for any probe found to match the given gene.   gene_info gives the info for the probes found to match.  The match is made based on a grep of the gene symbol to any column with "gene" in the label (case insensitive).  This is very liberal and will probably pull out may probes you didn't want.   It will be important to display the matches made and encourage manual review.
getSingleGeneIndexGrepGPL <- function(gene, data) {
	D=length(data)
	MaybeGene <- which(apply(data[[1]]@featureData@data, 2, function(i) sum(is.element(as.character(i), c('BRCA1', 'BRCA2', 'CTCF', 'FTO')))) > 3)[1] #Finds out which column is the gene column
	gene_text <- try(as.character(data[[1]]@featureData@data[,MaybeGene])) # Makes it into text
	try(grep(paste0('\\<',gene,'\\>'), gene_text)) # returns index
}


#' Return index marking location of probes corresponding to a specified gene
#' @param gene A gene annotation symbol
#' @param data  A list of gene expression sets
#' @param gene_maps a list of GPL objects giving matching gene annotation information for the expression sets given in data
#' @return A list with two elements.  Each element is a list of the same length as data and gene_maps.  gene_data gives the gene expression data for any probe found to match the given gene.   gene_info gives the info for the probes found to match.  The match is made based on a grep of the gene symbol to any column with "gene" in the label (case insensitive).  This is very liberal and will probably pull out may probes you didn't want.   It will be important to display the matches made and encourage manual review.
getSingleGeneIndex450K <- function(gene, data) {
        D=length(data)
        index <- which(sapply(data$GSE32146@featureData@data$UCSC_RefGene_Name, function(i) sum(strsplit(as.character(i), ';')[[1]]==gene) > 0)) #index of any match to gene
	return(index)
	}




#' Extract gene expression data for a given gene
#' @param gene A gene annotation symbol
#' @param data  A list of gene expression sets
#' @param BioCdb This links the annotation infomation given in data to the appropriate bioconductor package for annotation of that particular platform.
#' @eturn A list with three elements.  Each element is a list of the same length as data and gene_maps.  gene_data gives the gene expression data for any probe found to match the given gene.   gene_info gives the info for the probes found to match.  The match is made by the curated bioconductor package available for that array platform.
getSingleGeneData <- function(gene, data, BioCdb) {
	D=length(data)
	output <- list(gene_data=lapply(1:D, function(i) NULL), gene_info=lapply(1:D, function(i) NULL))
	for(i in 1:D) {
		package <-  paste0(BioCdb[BioCdb$gpl == Biobase::annotation(data[[i]]),'bioc_package'], '.db')
		if(!require(package, character.only=T)) {
			source("http://bioconductor.org/biocLite.R")
			biocLite(package, suppressUpdates=T, lib.loc='~/R')
			library(package, character.only=T)
  		}
		output$gene_info[[i]] <- try(select(get(package), 
		                                    keys=gene, 
		                                    columns=c("ENTREZID","GENENAME","SYMBOL","PROBEID"), 
		                                    keytype="ALIAS"), 
		                             silent=T)
		if(class(output$gene_info[[i]]) != 'try-error') output$gene_data[[i]] <- try((data[[i]][output$gene_info[[i]]$PROBEID,]), silent=T)
		if(class(output$gene_data[[i]]) == 'try-error') {
			warning(paste('Warning: BioconductorDB probeID is not a feature label in this dataset.  Searching for a match in all feature metadata.'))
			index <- is.element(colnames(data[[i]]@featureData@data),c('Gene Symbol'))
			anyColumnMatchIndex <-  try(unique(na.omit(as.integer(apply(data[[i]]@featureData@data[index,], 2, function(j) match(output$gene_info[[i]]$PROBEID, j))))), silent=T)
			if(length(anyColumnMatchIndex) > 0) output$gene_data[[i]] <- data[[i]][anyColumnMatchIndex,]
		}
	}
	return(output)
}


#' Extract gene expression data for a given gene
#' @param gene A gene annotation symbol
#' @param data  A list of gene expression sets
#' @param BioCdb This links the annotation infomation given in data to the appropriate bioconductor package for annotation of that particular platform.
#' @return A list with three elements.  Each element is a list of the same length as data and gene_maps.  gene_data gives the gene expression data for any probe found to match the given gene.   gene_info gives the info for the probes found to match.  The match is made by the curated bioconductor package available for that array platform.
getSingleGeneIndex <- function(gene, data, BioCdb) {
	D=length(data)
	output <- list(indices=lapply(1:D, function(i) NULL))
	names(output$indices) <- names(data)
	for(i in 1:D) {
		package <-  paste0(BioCdb[BioCdb$gpl == Biobase::annotation(data[[i]]),'bioc_package'], '.db')
		if(data[[i]]@annotation == 'GPL13534') { #Returns a gene index for 450K array
                                print(paste('warning:', names(data)[i], 'Using 450K gene annotations'))
                                output$indices[[i]] <- getSingleGeneIndexGrepGPL(gene, data[i])
                                print('Does this look right?')
                                print(head(knitr::kable(data[[i]]@featureData@data[output$indices[[i]],])))

		} else {
			if(package == '.db') {  #This section deals with cases where there is no package but does not yet return a gene index
				print(paste('warning:', names(data)[i], 'does not have a bioc database, gene to probe matches may be incorrect'))
				output$indices[[i]] <- getSingleGeneIndexGrepGPL(gene, data[i])
				print('Does this look right?')
				print(head(knitr::kable(data[[i]]@featureData@data[output$indices[[i]],])))
			} else {
				if(!require(package, character.only=T)) { #This installs the package if needed
					source("http://bioconductor.org/biocLite.R")
					biocLite(package, suppressUpdates=T, lib.loc='~/R', lib='~/R')
					library(package, character.only=T, lib.loc='~/R')
				}
				gene_info <- try(select(get(package), 
				                        keys=gene, 
				                        columns=c("ENTREZID","GENENAME","SYMBOL","PROBEID"), 
				                        keytype="ALIAS"), 
				                 silent=T)
				if(class(gene_info) != 'try-error') { #First we get the indeces of probes that match the gene
 					output$indices[[i]] <- which(is.element(rownames(data[[i]]), gene_info$PROBEID))
 					if(length(output$indices[[i]]) == 0) { # This deals with when probe ID is in a non standard column
 						warning(paste('Warning: BioconductorDB probeID is not a feature label in this dataset.  Searching for a match in all feature metadata.'))
 						output$indices[[i]] <-  unique(unlist(apply(data[[i]]@featureData@data, 2, function(j) which(is.element(as.character(j), gene_info$PROBEID)))))
						print('Does this look right?')
						print(head(knitr::kable(data[[i]]@featureData@data[output$indices[[i]],])))
  					}
				}
			}
		}
	}
	return(output)
}


#' Define groupings of subjects within a single gene expression set
#' @param datum A single gene expression set
#' @param columns a character vector giving column names of phenoData to search
#' @param regex a character vector where each term uniquely defines a grouping of subjects when applied to the phenoData columns in datum.  Regular expressions are allowed.
getSingleExprSetAssignments <- function(datum, columns, regex) {
	if(length(columns) > 1) {
		metatext <- do.call(paste, lapply(datum@phenoData@data[,columns], as.character))
	} else metatext <- as.character(datum@phenoData@data[,columns])
 	groupAssignments <- lapply(regex, function(r) {grepl(r, metatext)})
	return(groupAssignments)
}


#' Extract data underlying the groupings of subjects within a single gene expression set
#' @param datum A single gene expression set
#' @param columns a character vector giving column names of phenoData to search
#' @param regex a character vector where each term uniquely defines a grouping of subjects when applied to the phenoData columns in datum.  Regular expressions are allowed.
getSingleExprSetIndividualData <- function(datum, columns, regex) {
	if(length(columns) > 1) {
		metatext <- do.call(paste, lapply(datum@phenoData@data[,columns], as.character))
	} else metatext <- as.character(datum@phenoData@data[,columns])
	groupAssignments <- lapply(regex, function(r) {datum@phenoData@data[grepl(r, metatext),]})
	return(groupAssignments)
}


#' Define grouping of subjects for a list of expressionSets and a matching list of regular expressions which define groups within each expressionSetgrouping of subjects for a list of expressionSets and a matching list of regular expressions which define groups within each expressionSet
#' @param data A list of ExpressionSets
#' @param columnList A list of character vectors each of which gives column names of relevant data in phenoData to search for regEx match
#' @param regexList A list of character vectors each of which give regular expressions that define the groupings within the corresponding expressionSet
getGroupings <- function(data, columnList, regexList) {
	D=length(data)
	if (length(columnList) != D) stop('The data and columnList input are not the same length')
	if (length(regexList) != D) stop('The data and regexList input are not the same length')
	groupings <- lapply(1:D, function(i) getSingleExprSetAssignments(data[[i]], columnList[[i]], regexList[[i]]))
	return(groupings)
}

#' get raw data for datasets, and normalize with SCAN.UPC.  Replace ExpressionSet data with normalized and return
#' @param data A list of ExpressionSets
getNormGSEs <- function(datasets) {
  library(SCAN.UPC)
  library(parallel)
  normData <- mclapply(datasets, function(gse) {
    SCAN(gse)
  }, mc.cores=4)
  return(normData)
}



