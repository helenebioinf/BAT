#!/usr/bin/env Rscript

## Loading packages
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressWarnings(suppressPackageStartupMessages(library("gplots")))
suppressPackageStartupMessages(library("RColorBrewer"))

## Global options
## - 'VERSION' is used in order to potentially allow for different
## calculations of the performance measures and depend that on the
## version

VERSION <- "0.1"
VERSION.DATE <- "2015-05-11"

main <- function(){
  
  ## Easy command line parsing
  option_list <- list(
    make_option(c("-i", "--input.filename"), type="character", default="",
                help="name of input file", metavar="file"),
    make_option(c("-o", "--output.filename"), type="character", default="",
                help="prefix of output files", metavar="file"),
    make_option(c("-p", "--groups"), type="character", default="g1,g2",
                help="identifier of groups (comma-seperated)")
    )
  
  opt <- parse_args(OptionParser(option_list=option_list,
                               epilogue=sprintf("Version: %s (%s)\n",
                               VERSION, VERSION.DATE)))
  
  ## Few argument checks
  #print(str(opt))
  if (opt$input.filename != ""){
    if (file.access(opt$input.filename, 4) == -1){
      stop(sprintf("Specified file (%s) does not exist or is not readable",
                   opt$input.filename))
    }
  } else {
    stop("Required to specify -i (--input.filename)")
  }
  
  if (opt$output.filename == ""){
    stop("Required to specify -o (--output.filename)")
  }
  opt$pdf <- paste0(opt$output.filename, ".pdf")
  opt$txt <- paste0(opt$output.filename, ".txt") 
  
  ## Read input file
  cat("Read data .. ")
  data <- read.table(file=opt$input.filename, header=TRUE, comment.char="")
  cat("done.\n")
  
  ## Check number of columns, column names, and classes
  if (ncol(data) < 7){
    stop("Input file is invalid: incorrect number of columns, at least chr, start, end, unique identifier, grouping label, C.count and one methylation rate is required")
  }
	
  ## split group identifiers
  opt$group1.identifier <- unlist(strsplit(opt$groups, ","))[1]
  opt$group2.identifier <- unlist(strsplit(opt$groups, ","))[2]
  
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
  opt$group.length <- length(opt$g1)+length(opt$g2)

  
  ## Mean/diff methylation per group
  if(length(opt$g1) > 1){
	  data[,paste0('mean_', eval(parse(text="opt$group1.identifier")))] <- rowMeans(data[, opt$g1])
  }else{
  	  data[,paste0('mean_', eval(parse(text="opt$group1.identifier")))] <- data[, opt$g1]
  }
  if(length(opt$g2) > 1){
  	  data[,paste0('mean_', eval(parse(text="opt$group2.identifier")))] <- rowMeans(data[, opt$g2])
  }else{
  	  data[,paste0('mean_', eval(parse(text="opt$group2.identifier")))] <- data[, opt$g2]
  }
  
  #open pdf
  pdf(file=opt$pdf)
  cat("Length distributions .. ")
  myitem(opt, data)
  cat("done.\nBoxplots .. ")
  myboxplots(opt, data)
  cat("done.\nHeatmaps/Trees .. ")
  myheatmap(opt, data)
  cat("done.\nWrite tables .. ")
  myprint(opt, data)
  cat("done.\n")
  dev.off()
}

#boxplots of mean methylation level per item per sample
myboxplots <- function(opt, data){
	suppressMessages(suppressWarnings(print(ggplot(melt(data[-c(1:4,6)]), aes(x=variable, y=value)) + geom_boxplot(width=0.8, outlier.shape=NA, fill=c(rep("#9E9E9E",length(opt$g1)),rep("#009e73",length(opt$g2)),"#9E9E9E","#009e73")) + theme_cfg + theme(axis.text.x = element_text(angle=90, hjust=1)) + xlab('') + ylab('average methylation rate') + ggtitle('Average methylation rate of\nall annotation items per sample'))))
	if (length(unique(data$item)) > 1) {
		lapply(unique(data$item), function(x) {
			tmp <- subset(data, item==x)
			suppressMessages(suppressWarnings(print(ggplot(melt(tmp[-c(1:4,6)]), aes(x=variable, y=value)) + geom_boxplot(width=0.8, outlier.shape=NA, fill=c(rep("#9E9E9E",length(opt$g1)),rep("#009e73",length(opt$g2)),"#9E9E9E","#009e73")) + theme_cfg + theme(axis.text.x = element_text(angle=90, hjust=1)) + xlab('') + ylab('average methylation rate') + ggtitle(paste0('Average methylation rate\nper ', x, ' annotation item per sample')))))
		})
	}
}

#heatmap with methylation rates per samples for each item
myheatmap <- function(opt, data){
	if (length(unique(data$item)) > 1) {
		lapply(unique(data$item), function(x) {
			tmp <- data[data$item==x,-c(1:6)]
			row.names(tmp) <- data[data$item==x,4]
			if (nrow(tmp)>1 && ncol(tmp)>1){
				if (nrow(tmp)>50000){
					d <- dist(t(tmp[,1:opt$group.length]), method="euclidean")
					fit <- hclust(d, method="ward.D2")
					plot(fit, xlab='', sub='', main=paste('Hierachical Clustering - ', x))
				}
				else{
					print(heatmap.2(as.matrix(tmp[,1:opt$group.length]), Colv=T, Rowv=T, col=c(rev(brewer.pal(9,"Blues")),brewer.pal(9,"Reds")), scale='none', trace='none', key=TRUE, density.info='none', dendrogram='col', cexRow=0.6, cexCol=0.6, margins=c(8,8), main=x, ColSideColors=c(rep("#9E9E9E",length(opt$g1)),rep("#009e73",length(opt$g2)))))
				}
			}
		})
	}
	tmp <- data[,-c(1:6)]
	row.names(tmp) <- data[,4]
	if (nrow(tmp)>50000){
		d <- dist(t(tmp[,1:opt$group.length]), method="euclidean")
		fit <- hclust(d, method="ward.D2")
		plot(fit, xlab='', sub='', main='Hierachical Clustering all')
	}
	else{
		print(heatmap.2(as.matrix(tmp[,1:opt$group.length]), Colv=T, Rowv=T, col=c(rev(brewer.pal(9,"Blues")),brewer.pal(9,"Reds")), scale='none', trace='none', key=TRUE, density.info='none', dendrogram='col', labRow = NULL, cexCol=0.6, margins=c(8,8), main="all", ColSideColors=c(rep("#9E9E9E",length(opt$g1)),rep("#009e73",length(opt$g2)))))
	}
}

#describe items
myitem <- function(opt, data){
	if (length(unique(data$item)) < 9) {
		suppressWarnings(print(ggplot(data, aes(x=count)) + geom_line(aes(color=item), stat="density", size=1) + xlab("item length 0-99% quantile [Cs]") + theme_cfg + scale_color_brewer(palette="Dark2") + scale_x_continuous(limits=quantile(data$count, seq(0,1,0.01))[c(1,100)]) + theme(axis.text.x = element_text(angle=60, hjust=1)) + ggtitle('Length of covered\nannotation items')))
		suppressWarnings(print(ggplot(data, aes(x=(end-start))) + geom_line(aes(color=item), stat="density", size=1) + xlab("item length 0-99% quantile [nt]") + theme_cfg + scale_color_brewer(palette="Dark2") + scale_x_continuous(limits=quantile((data$end-data$start), seq(0,1,0.01))[c(1,100)]) + theme(axis.text.x = element_text(angle=60, hjust=1)) + ggtitle(paste('Length of covered\nannotation annotation items'))))
	} else {
		lapply(unique(data$item), function(x) {
			tmp <- subset(data, item==x)
			suppressWarnings(print(ggplot(tmp, aes(x=count)) + geom_line(color="#A50F15", stat="density", size=1) + xlab("item length 0-99% quantile [Cs]") + theme_cfg + scale_x_continuous(limits=quantile(tmp$count, seq(0,1,0.01))[c(1,100)]) + theme(axis.text.x = element_text(angle=60, hjust=1)) + ggtitle(paste('Length of ', x, ' annotation items'))))
			suppressWarnings(print(ggplot(tmp, aes(x=(end-start))) + geom_line(color="#A50F15", stat="density", size=1) + xlab("item length 0-99% quantile [nt]") + theme_cfg + scale_color_brewer(palette="Dark2") + scale_x_continuous(limits=quantile((tmp$end-tmp$start), seq(0,1,0.01))[c(1,100)]) + theme(axis.text.x = element_text(angle=60, hjust=1)) + ggtitle(paste('Length of ', x, ' annotation items'))))
			})
	}
	lapply(unique(data$item), function(x) {
		tmp <- subset(data, item==x)
		if (nrow(tmp)>1 && ncol(tmp)>1){
			suppressWarnings(print(ggplot(tmp, aes(x=(end-start), y=count)) + geom_point(color="#A50F15") + xlab("item length 0-99% quantile [nt]") + ylab("item length 0-99% quantile [Cs]") + theme_cfg + theme(axis.text.x = element_text(angle=60, hjust=1)) + ggtitle(paste('Length of ', x, ' annotation items')) + scale_x_continuous(limits=quantile((tmp$end-tmp$start), seq(0,1,0.01))[c(1,100)]) + scale_y_continuous(limits=quantile(tmp$count, seq(0,1,0.01))[c(1,100)])))
		}
	})
}

#write regions
myprint <- function(opt, data){
	write.table(file=opt$txt, format(data, scientific=F, trim=T), col.names=T, row.names=F, quote=F, sep="\t")
}

## make nice plots
theme_cfg <- theme_bw() + theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line = element_line(colour = "black"),
    legend.key = element_blank(),
    axis.line.x = element_line(colour="black"),
    axis.line.y = element_line(colour="black")
  )




main()
