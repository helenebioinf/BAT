#!/usr/bin/env Rscript

## Loading packages
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))

## Global options
## - 'VERSION' is used in order to potentially allow for different
## calculations of the performance measures and depend that on the
## version

VERSION <- "0.1"
VERSION.DATE <- "2015-08-06"

main <- function(){
  cat(paste0("[INFO]\t", format(Sys.time(), "%a %b %d, %X, %Y"), "\tBAT_overview v", VERSION, " started\n"))
  
  ## Easy command line parsing
  option_list <- list(
    make_option(c("-i", "--input.filename"), type="character", default="",
                help="name of input file (summary file)", metavar="file"),
    make_option(c("-o", "--output.prefix"), type="character", default="",
                help="path/prefix of output file", metavar="file"),
    make_option(c("-g", "--groups"), type="character", default="g1,g2",
                help="identifier of first and second group, seperated by ,
                columns need to begin with one of the identifiers"),
    make_option(c("-m", "--mis"), type="character", default="NA",
                help="identifer for missing values, default NA")
    )
  
  opt <- parse_args(OptionParser(option_list=option_list,
                               epilogue=sprintf("Version: %s (%s)\n",
                               VERSION, VERSION.DATE)))

  ## Few argument checks
  #print(str(opt))
  cat(paste0("[INFO]\t", format(Sys.time(), "%a %b %d, %X, %Y"), "\tChecking flags\n"))
  
  if (opt$input.filename != ""){
    if (file.access(opt$input.filename, 4) == -1){
      stop(sprintf("Specified file (%s) does not exist or is not readable",
                   opt$input.filename))
    }
  } else {
    stop("Required to specify -i (--input.filename)")
  }
  
  if (opt$output.prefix == ""){
    stop("Required to specify -o (--output.prefix)")
  }
  opt$output.filename <- paste(opt$output.prefix, ".pdf", sep="")
  
  ## split group identifiers
  opt$group1.identifier <- unlist(strsplit(opt$groups, ","))[1]
  opt$group2.identifier <- unlist(strsplit(opt$groups, ","))[2]
  
  ## Read input file
  cat(paste0("[INFO]\t", format(Sys.time(), "%a %b %d, %X, %Y"), "\tReading input ", opt$input.filename, "\n"))
  data <- read.table(file=opt$input.filename, header=TRUE, stringsAsFactors=FALSE, comment.char="", na.strings=opt$mis)
  
  ## Check number of columns, column names, and classes
  if (ncol(data) < 2){
    stop("Input file is invalid: incorrect number of columns, at least methylation rates of one sample per group is required")
  }
	
  ## Get group dependend columns
  opt$g1;
  opt$g2;
  opt$g1 <- grep(paste("^", opt$group1.identifier, sep=""),colnames(data))
  opt$g2 <- grep(paste("^", opt$group2.identifier, sep=""),colnames(data))
  if(length(opt$g1) == 0){
    stop("Required group identfier for group 1 not found in input.")
  }
  if(length(opt$g2) == 0){
    stop("Required group identfier for group 2 not found in input.")
  }
  
  ## Mean/diff methylation per group
  if(length(opt$g1) > 1){
    data$mean1 <- rowMeans(data[,opt$g1], na.rm=TRUE)
  } else {
    data$mean1 <- data[,opt$g1]
  }
  if(length(opt$g2) > 1){
    data$mean2 <- rowMeans(data[,opt$g2], na.rm=TRUE)
  } else {
    data$mean2 <- data[,opt$g2]
  }
  data$diff <- data$mean1-data$mean2
  
  #open pdf
  cat(paste0("[INFO]\t", format(Sys.time(), "%a %b %d, %X, %Y"), "\tPlot statistics to ", opt$output.filename, "\n"))
  pdf(file=opt$output.filename)
  myavgmethyl(opt, data)
  myhierclust(opt, data)
  mybinwiselevel(opt, data)
  mymean(opt, data)
  myscatter(opt, data)
  mydiff(opt, data)
  invisible(dev.off())
}

#hierarchical clustering
myhierclust <- function(opt, data){
	d <- dist(t(data[,c(opt$g1,opt$g2)]), method="euclidean")
	fit <- hclust(d, method="ward.D2")
	plot(fit, xlab='', sub='', main='Hierachical Clustering')
}

#binwise average methylation level distribution
mybinwiselevel <- function(opt, data){
	data[,as.character(opt$group1.identifier)] <- cut(data$mean1, breaks=seq(0,1,0.2), right=F, include.lowest=T)
	data[,as.character(opt$group2.identifier)] <- cut(data$mean2, breaks=seq(0,1,0.2), right=F, include.lowest=T)
	suppressWarnings(print(ggplot(na.omit(reshape(data[,c(opt$group1.identifier,opt$group2.identifier)], varying=c(opt$group1.identifier,opt$group2.identifier), v.names="bins", direction="long", timevar="group", times=c(opt$group1.identifier,opt$group2.identifier))), aes(x=group, fill=bins)) + geom_bar(aes(y=(..count..)), width=0.5) + scale_fill_manual(values=c("#FEE5D9","#FCAE91","#FB6A4A","#DE2D26","#A50F15")) + theme_bw(14) + theme(text=element_text(family='Helvetica'), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black")) + xlab('') + ylab('count') + ggtitle("Binned group-average methylation rate")))
}

#average methylation level per groups
myavgmethyl <- function(opt, data){
	tmp <- melt(colMeans(data[,c(colnames(data)[c(opt$g1,opt$g2)])], na.rm=T), id=1)
	tmp$group <- ifelse(grepl(paste("^", opt$group1.identifier, sep=""),rownames(tmp)), tmp$group <- opt$group1.identifier, tmp$group <- opt$group2.identifier)
	suppressWarnings(print(ggplot(tmp, aes(x=group, y=value)) + geom_boxplot(aes(fill=group), width=0.5) + theme_bw(14) + scale_fill_manual(values=c("#9E9E9E","#009e73")) + theme(legend.position='none', panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black"), legend.key=element_blank()) + xlab('') + ylab('average methylation rate') + ggtitle("Average methylation rate per sample")))
	tmp <- NULL
}


#smoothed scatterplot of methylation levels
myscatter <- function(opt, data){
	Lab.palette <- colorRampPalette(c("darkblue", "lightblue", "green", "yellow", "darkred"), space = "Lab")
	summary <- summary(lm(data$mean1~data$mean2))
	smoothScatter(data$mean1, data$mean2, colramp=Lab.palette, xlab=paste("average methylation",opt$group1.identifier,sep=" "), ylab=paste("average methylation",opt$group2.identifier,sep=" "), sub=paste("adjusted R-squared:",round(summary$adj.r.squared,4),sep=""), main="Scatter plot of group-average methylation rates")
}

#average methylation
mymean <- function(opt, data){
	suppressWarnings(print(ggplot(na.omit(reshape(data[,c("mean1","mean2")], varying=c("mean1","mean2"), v.names="rate", direction="long", timevar="group", times=c("mean1","mean2"))), aes(x=rate, color=group)) + geom_line(stat="density", size=1.3) + theme_bw(14) + scale_color_manual(values=c("#9E9E9E","#009e73"), guide=guide_legend(title=NULL), labels=c(opt$group1.identifier,opt$group2.identifier)) + theme(panel.grid=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black")) + xlab('average methylation') + scale_x_continuous(limits=c(0,1)) + ggtitle("Group-average methylation rate distribution")))
}

#difference of average methylation
mydiff <- function(opt, data){
	suppressWarnings(print(ggplot(na.omit(data), aes(x=diff)) + geom_histogram(aes(fill=..count..), binwidth=0.02) + scale_fill_gradient(low="#A50F15", high="#FEE5D9") + theme_bw(14) + theme(panel.grid=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black")) + xlab(paste("methylation difference", opt$group1.identifier, "-", opt$group2.identifier, sep=" ")) + scale_x_continuous(limits=c(-1,1)) + ggtitle("C-wise difference of\ngroup-average methylation rates")))
}

#difference of average methylation
mycirco <- function(opt, data){
	suppressWarnings(print(ggplot(na.omit(data), aes(x=diff)) + geom_histogram(aes(fill=..count..), binwidth=0.02) + scale_fill_gradient(low="#A50F15", high="#FEE5D9") + theme_bw(14) + theme(panel.grid=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.line.x=element_line(colour="black"), axis.line.y=element_line(colour="black")) + xlab(paste("methylation difference", opt$group1.identifier, "-", opt$group2.identifier, sep=" ")) + scale_x_continuous(limits=c(-1,1)) + ggtitle("C-wise difference of\ngroup-average methylation rates")))
}




main()
