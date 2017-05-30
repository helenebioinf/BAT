#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

## Loading packages
library(ggplot2)

## make nice plots
theme_cfg <-
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    legend.key = element_blank(),
    axis.line.x = element_line(colour="black"),
    axis.line.y = element_line(colour="black")
  )

## Read input file
data <- read.table(file=args[1], header=F, col.names=c('chr','start','end','q_val','diff','C','group1','group2'))
data$length <- data$end-data$start

## Plot statistics
pdf(args[2])
#difference histogram
suppressWarnings(ggplot(data, aes(x=diff)) + geom_histogram(binwidth=0.02, fill='darkgreen', color='black') + xlab(paste0("mean methylation difference\n", args[3], " - ", args[4])) + theme_cfg + scale_x_continuous(limits=c(-1,1)))
#q_val vs difference
suppressWarnings(ggplot(data, aes(x=diff, y=q_val)) + geom_point(alpha=.5, color='darkgreen') + scale_y_log10() + xlab(paste0("mean methylation difference\n", args[3], " - ", args[4])) + ylab("q-value") + theme_cfg + scale_x_continuous(limits=c(-1,1)))
#mean1 vs mean2
suppressWarnings(ggplot(data, aes(x=group1, y=group2)) + geom_point(alpha=.5, color='darkgreen') + coord_fixed() + xlab(paste0("methylation rate in DMRs - ", args[3])) + ylab(paste0("methylation rate in DMRs - ", args[4])) + theme_cfg + scale_x_continuous(limits=c(0,1)) + scale_y_continuous(limits=c(0,1)))
#length distribution Cs
suppressWarnings(ggplot(data, aes(x=length)) + geom_line(stat="density", size=1, color='darkgreen') + xlab("DMR length [nt]") + theme_cfg)
#length distribution nt
suppressWarnings(ggplot(data, aes(x=C)) + geom_line(stat="density", size=1, color='darkgreen') + xlab("DMR length [Cs]") + theme_cfg)
#nt vs Cs
suppressWarnings(ggplot(data, aes(x=length, y=C)) + geom_point(alpha=.5, color='darkgreen') + xlab("DMR length [nt]") + ylab("DMR length [Cs]") + theme_cfg)
invisible(dev.off())
