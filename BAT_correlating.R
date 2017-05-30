#!/usr/bin/env Rscript
setEPS()

## Loading packages
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("reshape2"))
suppressPackageStartupMessages(library("plyr"))
suppressPackageStartupMessages(library("grid"))

## Global options
## - 'VERSION' is used in order to potentially allow for different
## calculations of the performance measures and depend that on the
## version

VERSION <- "0.1"
VERSION.DATE <- "2015-08-25"

main <- function(){
  
  ## Easy command line parsing
  option_list <- list(
    make_option(c("-e", "--expression.filename"), type="character", default="",
                help="path/file of expression input file", metavar="file"),
    make_option(c("-m", "--methylation.filename"), type="character", default="",
                help="path/file of methylation input file", metavar="file"),
    
    make_option(c("-s", "--statistics.filename"), type="character", default="",
                help="path/prefix of statistics file", metavar="file"),            
    make_option(c("-o", "--output.filename"), type="character", default="",
                help="path/prefix of output file", metavar="file"),
    make_option(c("-p", "--group1.identifier"), type="character", default="",
                help="identifier of first group and columns need to begin with the identifier"),
    make_option(c("-q", "--group2.identifier"), type="character", default="",
                help="identifier of second group and columns need to begin with the identifier"),
    make_option(c("-r", "--group3.identifier"), type="character", default="",
                help="identifier of third group and columns need to begin with the identifier")
    )
  
  opt <- parse_args(OptionParser(option_list=option_list,
                               epilogue=sprintf("Version: %s (%s)\n",
                               VERSION, VERSION.DATE)))
  
  opt$colors <- c('#7f7d7a','#6da6c8','#2d8c00')
  
  ## Few argument checks
  #print(str(opt))
  if (opt$expression.filename != ""){
    if (file.access(opt$expression.filename, 4) == -1){
      stop(sprintf("Specified file (%s) does not exist or is not readable",
                   opt$expression.filename))
    }
  } else {
    stop("Required to specify -e (--expression.filename)")
  }
  
  if (opt$methylation.filename != ""){
    if (file.access(opt$methylation.filename, 4) == -1){
      stop(sprintf("Specified file (%s) does not exist or is not readable",
                   opt$methylation.filename))
    }
  } else {
    stop("Required to specify -m (--methylation.filename)")
  }
  
  if (opt$output.filename == ""){
    stop("Required to specify -o (--output.filename)")
  }
  
  if (opt$statistics.filename == ""){
    stop("Required to specify -s (--statistics.filename)")
  }
  
  if (opt$group1.identifier == ""){
    stop("Required to specify -p (--group1.identifier)")
  }
  if (opt$group2.identifier == ""){
    stop("Required to specify -q (--group2.identifier)")
  }
  
  # prepare output statistics file
  text <- "plot.pdf\tENSG_ID\tmethylation_region\tadj.R2\trho\trho.p-val"
  write.table(file=opt$statistics.filename, text, row.names=F, quote=F, sep='\t')
  
  ## Read input file
  expr <- read.table(file=opt$expression.filename, col.names=c("sample","group","id","expression"))
  methyl <- read.table(file=opt$methylation.filename, col.names=c("sample","group","region","id","methylation"))
  
  data <- merge(expr, methyl, by=c("sample","group","id"))
  ddply(data, .(id, region), function(x){myOut(x, opt)})
}

myOut <- function(data, opt){
	out <- NULL
	out$plot_out	<- paste(opt$output.filename, "_", data$region[1], "_", data$id[1], ".pdf", sep="")
	out$stats_out	<- opt$statistics.filename;
   out$text_tmp	<- paste(opt$output.filename, "_", data$region[1], "_", data$id[1], ".text.tmp.eps", sep="")
	out$expr_tmp	<- paste(opt$output.filename, "_", data$region[1], "_", data$id[1], ".expr.tmp.eps", sep="")
	out$methyl_tmp	<- paste(opt$output.filename, "_", data$region[1], "_", data$id[1], ".methyl.tmp.eps", sep="")
	out$cor_tmp		<- paste(opt$output.filename, "_", data$region[1], "_", data$id[1], ".cor.tmp.eps", sep="")
  	 
  	## rows
	out$g1 <- grep(opt$group1.identifier, data$group)
	out$g2 <- grep(opt$group2.identifier, data$group)
	out$g3 <- NULL
	if (opt$group3.identifier != "") {
		out$g3 <- grep(opt$group3.identifier, data$group)
	}
	
	## Statistics
	text <- myStat(data[c(out$g1, out$g2),], opt, out)
	
	## Correlation sub-plots
	myExpr(data[c(out$g1, out$g2, out$g3),], opt, out)
	myMethyl(data[c(out$g1, out$g2, out$g3),], opt, out)
	myCorr(data[c(out$g1, out$g2, out$g3),], opt, out)
}

myStat <- function(data, opt, out){
	# statistics
	r2 <- round(summary(lm(expression~methylation, data=data))$adj.r.squared, 2)
	test <- suppressWarnings(cor.test(data$expression, data$methylation, method="spearman"))
	
	text <- paste(unique(data$id), "\n", unique(data$region), "\nadj.R2=", r2, "\nrho=", round(test$estimate, 2), "\nrho.p-val=", round(test$p.value, 4), sep="")
	stats <- paste(unique(out$plot_out), "\t", unique(data$id), "\t", unique(data$region), "\t", r2, "\t", round(test$estimate, 2), "\t", round(test$p.value, 4), sep="")
    
	write.table(file=out$stats_out, stats, col.names=F, row.names=F, quote=F, sep='\t', append=T)
	
	postscript(out$text_tmp, width=3, height=3)
		plot(-5:5, -5:5, type='n', xlab=(''), ylab=(''), xaxt='n', yaxt='n', axes=F)
		text(0,0,text)
   invisible(dev.off())
}

myExpr <- function(data, opt, out){
	postscript(out$expr_tmp, width=3.5)
	suppressWarnings(print(ggplot(data, aes(group, expression)) + stat_boxplot(geom ='errorbar', width=0.25) + geom_boxplot(aes(fill=group)) + theme_bw(18) + xlab('') + ylab('expression') + scale_fill_manual(values=opt$colors) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), legend.position='none', panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.line.x = element_line(colour="black"), axis.line.y = element_line(colour="black"), plot.margin=unit(c(0,0,0,0),'cm'))))
	invisible(dev.off())
}

myMethyl <- function(data, opt, out){
	postscript(out$methyl_tmp, height=3.5)
   suppressWarnings(print(ggplot(data, aes(group, methylation)) + stat_boxplot(geom ='errorbar', width=0.25) + geom_boxplot(aes(fill=group)) + theme_bw(18) + xlab('') + ylab('methylation') + scale_fill_manual(values=opt$colors) + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.position='none', panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.line.x = element_line(colour="black"), axis.line.y = element_line(colour="black"), plot.margin=unit(c(0,0.5,0,0),'cm')) + coord_flip(ylim=c(-0.05,1.05))))
	invisible(dev.off())
}

myCorr <- function(data, opt, out){
	postscript(out$cor_tmp)
   suppressWarnings(print(ggplot() + geom_point(data=data[c(out$g1,out$g2),], aes(x=methylation, y=expression, color=group), size=3) + theme_bw(18) + xlab('') + ylab('') + coord_cartesian(xlim=c(-0.05,1.05)) + scale_color_manual(values=opt$colors) + theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), legend.title=element_blank(), legend.key=element_blank(), legend.position=c(0.9,0.9), panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border=element_blank(), panel.background=element_blank(), axis.line.x = element_line(colour="black"), axis.line.y = element_line(colour="black"), plot.margin=unit(c(0,0.5,0,0),'cm')) + stat_smooth(data=data[c(out$g1,out$g2),], aes(x=methylation, y=expression, color=group), method='lm', se=FALSE, colour='grey', lty='dashed') + guides(fill=guide_legend(override.aes=list(shape=22, size=3))) + geom_point(data=data[out$g3,], aes(x=methylation, y=expression, color=group), size=3)))
	invisible(dev.off())
}

main()
