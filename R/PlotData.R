#' Boxplot for a single gene, multiple studies, grouped by categories
#' @param gene_data A gene expression set object that has been filtered for just the gene of interest
#' @param groups A list of logical vectors which determine assignment of individuals in each study to different groups


bpOneGene <- function(gene_data, groups, ...) {
par(mar=c(10,5,5,5))
boxplot(c(lapply(1:length(gene_data), function(gene) gene_data[[gene]]@assayData$exprs[groups[[gene]][[1]]]), lapply(1:length(gene_data), function(gene) gene_data[[gene]]@assayData$exprs[groups[[gene]][[2]]])), las=3, ...)
}

#' Boxplot for a single gene, multiple studies, grouped by categeoris.  Boxplots now have seperate axes as most studies will not be directly comparable on the same scale.  If multiple probes match the gene, all are included in same box.
#' @param data A gene expression set
#' @param geneindex An integer index giving the rows of data that match a given gene
#' @param groups A list of logical vectors which determine assignment of individuals in each study to different groups
#' @param names A list of vectors corresponding to the structure of groups giving appropriate names for labeling groups and studies
bpOneGeneSeperate <- function(data, geneindex, groups, names=NULL, mainlist, ...) {
N <- ceiling(sqrt(length(data)))
layout(matrix(1:N^2, nrow=N, ncol=N, byrow=T))
for(i in 1:length(data)) {
  par(mar=c(10,5,5,5))
  # log2 transform
  ex <- Biobase::exprs(data[[i]][geneindex$indices[[i]],])
  if(dim(ex)[1] > 0) {
    graphics::boxplot(list(ex[,groups[[i]][[1]]], ex[,groups[[i]][[2]]]), las=3, names=names[[i]],main=mainlist[i],...)
  } else {
  plot.new(); plot.window(c(0,1), c(0,1)); text(0.5,0.5,'Meaningful Error Message')
  }
 }
}

#' Barplot for a single gene, multiple studies, grouped by categories. In these, probes are  plotted seperately if there are mutliple.
#' @param data A gene expression set
#' @param geneindex An integer index giving the rows of data that match a given gene
#' @param groups A list of logical vectors which determine assignment of individuals in each study to different groups
#' @param names A list of vectors corresponding to the structure of groups giving appropriate names for labeling groups and studies

barOneGene <- function(data, geneindex, groups, names=NULL, mainlist, colors,...) {
N <- ceiling(sqrt(length(data)))
layout(matrix(1:(2*N^2), ncol=N*2, nrow=N, byrow=T))
for(i in 1:length(data)) {
  #par(mar=c(3,5,7,0))
  # log2 transform
  ex <- Biobase::exprs(data[[i]][geneindex$indices[[i]],])
  labs <- rownames(ex)
  ylim<-c(min(as.numeric(ex), na.rm=T), max(as.numeric(ex), na.rm=T))
  ylim <- c(ylim[1] - abs(.05*diff(ylim)), ylim[2])
  if(length(labs) == 1) labs <- NULL
  if(dim(ex)[1] > 0 & sum(groups[[i]][[1]]) > 0) {
    graphics::barplot(t(ex[,groups[[i]][[1]]]), las=3, main=paste(mainlist[i], names[[i]][1], sep='\n'), col=colors[1],names.arg=labs,beside=T, border=NA,ylim=ylim,xpd=F,width=1,...)
    if(dim(ex)[1] > 1) lapply(1:(dim(ex)[1] - 1), function(k) abline(v=sum(groups[[i]][[1]])*k, col='black', lty=2, lwd=2))
  } else  {plot.new(); plot.window(c(0,1), c(0,1)); text(0.5,0.5,'Check data format or empty groups.\nPlot not rendered', cex=1.5)}
  if(dim(ex)[1] > 0 & sum(groups[[i]][[2]]) > 0) {
    graphics::barplot(t(ex[,groups[[i]][[2]]]), las=3, main=paste(mainlist[i], names[[i]][2], sep='\n'), col=colors[2],names.arg=labs, beside=T,border=NA,ylim=ylim,xpd=F,width=1,...)
    if(dim(ex)[1] > 1) lapply(1:(dim(ex)[1] - 1), function(k) abline(v=sum(groups[[i]][[2]])*k, col='black', lty=2, lwd=2))
  } else  {plot.new(); plot.window(c(0,1), c(0,1)); text(0.5,0.5,'Check data format or empty groups.\nPlot not rendered', cex=1.5)}
 }
}



#' Violin plot for a single gene, multiple studies, grouped by categories. In these, probes are  plotted seperately if there are mutliple.
#' @param data A gene expression set
#' @param geneindex An integer index giving the rows of data that match a given gene
#' @param groups A list of logical vectors which determine assignment of individuals in each study to different groups
#' @param names A list of vectors corresponding to the structure of groups giving appropriate names for labeling groups and studies

violinOneGene <- function(data, geneindex, groups, names=NULL, mainlist, colors,...) {
long_data <- do.call(rbind, lapply(1:length(data), function(i) {
  library(lattice)
  mydata<-data[[i]][geneindex$indices[[i]],]
  #groupStatus <- do.call(paste0, groups[[i]])
  groupStatus <- do.call(paste0, lapply(1:length(groups[[i]]), function(j) gsub('FALSE', '', gsub('TRUE', names(groups[[i]])[j], groups[[i]][[j]]))))
  groupStatus[groupStatus == ''] <- 'No\nAssignment'
  Biobase::pData(mydata)$GROUP <- as.factor(groupStatus)
  wide_data.df<-data.frame(
                        cbind(
                                Biobase::pData(mydata)[,c("GROUP","geo_accession")],
                                t(Biobase::exprs(mydata))
                                ), Dataset=names(data)[i]
                        )
  long_data <- reshape2::melt(wide_data.df, value.name = 'expr', variable.name="Probe.ID",
                 id.vars = c("GROUP","geo_accession", 'Dataset'))
  long_data <- data.frame(long_data, SetByID=paste(long_data$Dataset, long_data$Probe.ID, sep='_'))
  return(long_data)
}))
print(bwplot(expr ~ GROUP | SetByID,
data = long_data,
groups=GROUP,
xlab='GROUP',
ylab='Distribution of raw GEO data',
layout=c(2,ceiling(length(unlist(geneindex))/2)),
outer=FALSE,
as.table=TRUE,
horizontal=FALSE,
strip=myStripStyle,
col=colors,
scales=list(relation='free', x=list(rot=45)),
panel = panel.superpose,
panel.groups = panel.violin))
}

#' Box plot for a single gene, multiple studies, grouped by categories. In these, probes are  plotted seperately if there are mutliple.
#' @param data A gene expression set
#' @param geneindex An integer index giving the rows of data that match a given gene
#' @param groups A list of logical vectors which determine assignment of individuals in each study to different groups
#' @param names A list of vectors corresponding to the structure of groups giving appropriate names for labeling groups and studies

boxOneGene <- function(data, geneindex, groups, names=NULL, mainlist, colors,...) {
library(lattice)
long_data <- do.call(rbind, lapply(1:length(data), function(i) {
  mydata<-data[[i]][geneindex$indices[[i]],]
  #groupStatus <- do.call(paste0, groups[[i]])
  groupStatus <- do.call(paste0, lapply(1:length(groups[[i]]), function(j) gsub('FALSE', '', gsub('TRUE', names(groups[[i]])[j], groups[[i]][[j]]))))
  groupStatus[groupStatus == ''] <- 'No Assignment'
  Biobase::pData(mydata)$GROUP <- as.factor(groupStatus)
  print("boxOneGene - long_data - wide_data.df creation start")
  wide_data.df<-data.frame(
                        cbind(
                                Biobase::pData(mydata)[,c("GROUP","geo_accession")],
                                t(Biobase::exprs(mydata))
                                ), Dataset=names(data)[i]
                        )
  print("boxOneGene - long_data - wide_data.df creation done")
  print("boxOneGene - long_data - melting start")
  long_data <- reshape2::melt(wide_data.df, value.name = 'expr', variable.name="Probe.ID",
                 id.vars = c("GROUP","geo_accession", 'Dataset'))
  print("boxOneGene - long_data - melting done")
  print("boxOneGene - long_data - building long_data df start")
  long_data <- data.frame(long_data, SetByID=paste(long_data$Dataset, long_data$Probe.ID, sep='_'))
  print("boxOneGene - long_data - building long_data df done")
  return(long_data)
}))
print(paste0('length long_data:', length(long_data)))
print(paste0('ceiling: ceiling(length(long_data)/2): ',  ceiling(length(long_data)/2)))
print(bwplot(expr ~ GROUP | SetByID,
data = long_data,
layout=c(2,ceiling(length(unlist(geneindex))/2)),
xlab='GROUP',
ylab='Distribution of raw GEO data',
outer=FALSE,
notch=T,
as.table=TRUE,
horizontal=FALSE,
strip=myStripStyle,
scales=list(relation='free', x=list(rot=45)),
par.settings = list(box.rectangle = list(fill= colors))))
}

myStripStyle <- function(which.panel, factor.levels,...) {
  study <- as.integer(as.factor(sapply(factor.levels, function(i) strsplit(i, '_')[[1]][1])))
  bgColors=rainbow(length(unique(study)))[study]
  panel.rect(0, 0, 1, 1,
             col = bgColors[which.panel],
             border = 1)
  panel.text(x = 0.5, y = 0.5,
             font=2,
             lab = factor.levels[which.panel],
             col = 'white', cex=0.8)
  }

#' Box plot for a single gene, multiple studies, grouped by categories. In these, probes are  plotted seperately if there are mutliple.
#' @param data A gene expression set
#' @param geneindex An integer index giving the rows of data that match a given gene
#' @param groups A list of logical vectors which determine assignment of individuals in each study to different groups
#' @param names A list of vectors corresponding to the structure of groups giving appropriate names for labeling groups and studies

stripOneGene <- function(data, geneindex, groups, names=NULL, mainlist, colors,...) {
  library(lattice)
  long_data <- do.call(rbind, lapply(1:length(data), function(i) {
    mydata<-data[[i]][geneindex$indices[[i]],]
    #groupStatus <- do.call(paste0, groups[[i]])
    groupStatus <- do.call(paste0, lapply(1:length(groups[[i]]), function(j) gsub('FALSE', '', gsub('TRUE', names(groups[[i]])[j], groups[[i]][[j]]))))
    groupStatus[groupStatus == ''] <- 'No Assignment'
    Biobase::pData(mydata)$GROUP <- as.factor(groupStatus)
    wide_data.df<-data.frame(
      cbind(
        Biobase::pData(mydata)[,c("GROUP","geo_accession")],
        t(Biobase::exprs(mydata))
      ), Dataset=names(data)[i]
    )
    long_data <- reshape2::melt(wide_data.df, value.name = 'expr', variable.name="Probe.ID",
                                id.vars = c("GROUP","geo_accession", 'Dataset'))
    long_data <- data.frame(long_data, SetByID=paste(long_data$Dataset, long_data$Probe.ID, sep='_'))
    return(long_data)
  }))
  bgcolors <- rainbow(length(unique(long_data$SetByID)))
  print(stripplot(expr ~ GROUP | SetByID,
               data = long_data,
               xlab='GROUP',
               ylab='Distribution of raw GEO data',
               layout=c(2,ceiling(length(unlist(geneindex))/2)),
               groups=GROUP,
               outer=FALSE,
               jitter=T,
               as.table=TRUE,
               horizontal=FALSE,
               col=colors,
               strip=myStripStyle,
               scales=list(relation='free', x=list(rot=45))))
}


#' Make full array of raw data grouped by study and by probeID
#' @param data A gene expression set
#' @param geneindex An integer index giving the rows of data that match a given gene
#' @param groups A list of logical vectors which determine assignment of individuals in each study to different groups
#' @param names A list of vectors corresponding to the structure of groups giving appropriate names for labeling groups and studies
tableOneGene <- function(data, geneindex, groups) {
library(lattice)
long_data <- do.call(rbind, lapply(1:length(data), function(i) {
  mydata<-data[[i]][geneindex$indices[[i]],]
  #groupStatus <- do.call(paste0, groups[[i]])
  groupStatus <- do.call(paste0, lapply(1:length(groups[[i]]), function(j) gsub('FALSE', '', gsub('TRUE', names(groups[[i]])[j], groups[[i]][[j]]))))
  groupStatus[groupStatus == ''] <- 'No Assignment'
  Biobase::pData(mydata)$GROUP <- as.factor(groupStatus)
  wide_data.df<-data.frame(
                        cbind(
                                Biobase::pData(mydata)[,c("GROUP","geo_accession")],
                                t(Biobase::exprs(mydata))
                                ), Dataset=names(data)[i]
                        )
  long_data <- reshape2::melt(wide_data.df, value.name = 'expr', variable.name="Probe.ID",
                 id.vars = c("GROUP","geo_accession", 'Dataset'))
  return(long_data)
}))
return((long_data))
}


#' Make ragged array of raw data grouped by study and grouping but individual de-identified
#' @param data A gene expression set
#' @param geneindex An integer index giving the rows of data that match a given gene
#' @param groups A list of logical vectors which determine assignment of individuals in each study to different groups
#' @param names A list of vectors corresponding to the structure of groups giving appropriate names for labeling groups and studies
raggedOneGene <- function(data, geneindex, groups) {
ragged_data <-lapply(1:length(data), function(i) {
  mydata<-lapply(1:length(groups[[i]]), function(j) {
  ss <-(Biobase::exprs(data[[i]][geneindex$indices[[i]],groups[[i]][[j]]]))
  tt <- lapply(1:dim(ss)[1], function(i) {x<-ss[i,]; return(x)})
  names(tt) <- rownames(ss)
  return(tt)
  })
  names(mydata) <- names(groups[[i]])
  return((mydata))
  })
names(ragged_data) <- names(groups)
ragged_data <-  unlist(unlist(ragged_data, recursive=F), recursive=F)
maxlength<-max(sapply(ragged_data, length))
for(i in 1:length(ragged_data)) {length(ragged_data[[i]]) <- maxlength}
ragged_data <- data.frame(ragged_data)
rownames(ragged_data) <- NULL
return(data.frame(ragged_data))
}

#' Make array of individual assignments
#' @param data A gene expression set
#' @param groups A list of logical vectors which determine assignment of individuals in each study to different groups
#' @param names A list of vectors corresponding to the structure of groups giving appropriate names for labeling groups and studies
groupingDataset <- function(data, dataset, groups) {
groupStatus <- do.call(paste0, lapply(1:length(groups[[dataset]]), function(j) gsub('FALSE', '', gsub('TRUE', names(groups[[dataset]])[j], groups[[dataset]][[j]]))))
groupStatus[groupStatus == ''] <- 'No Assignment'
return(cbind(groupStatus, Biobase::pData(data[[dataset]])))
}

