# BAT - Bisulfite Analyss Toolkit

[WEBSITE](http://www.bioinf.uni-leipzig.de/Software/BAT/)

## Introduction

Cytosine DNA methylation is a biochemical process that has been shown
to play an important roll in gene expression and cell
differentiation. Recently, a number of whole-genome bisulfite
sequencing (WGBS) and targeted bisulfite sequencing (i.e., RRBS)
protocols have made it possible to precisely and accurately capture
this major epigenetic modification.

Here, a modular bisulfite analysis toolkit (BAT) is introduced. It
tackles the major tasks for analyzing bisulfite sequencing data:
mapping, extraction of the methylation information (referred to as
methylation calling), and differential methylation analysis as well as
downstream analyses like integration of the methylation data with
annotation and gene expression data. Each part of this analysis
workflow is modular and can easily be customized or extended by other
bisulfite- or NGS-related tools, but can also be used *as is* with the
additional benefit of many automatically generated graphics by the
modules of BAT.


## Modules
[**Mapping**](http://www.bioinf.uni-leipzig.de/Software/BAT/mapping)

The first module comprises read mapping including pre- and
postprocessings. This includes conversion to BAM format, mapping 
statistics and merging of multiple mapping runs.

[**Calling**](http://www.bioinf.uni-leipzig.de/Software/BAT/calling)

The second module covers the extraction of methylation information from 
the alignments, filtering for positions of interest, e.g. CG context, 
and conversion of methylation information for visualisation.

[**Analysis**](http://www.bioinf.uni-leipzig.de/Software/BAT/analysis)

In the third module basic analysis of two groups of a single sample or 
up to multiple samples are performed.

[**DMRs**](http://www.bioinf.uni-leipzig.de/Software/BAT/dmrs)

Finally, the calling of DMRs is coverd by the fourth module. Basic statistics 
of the DMRs are provided and given expression information of genes, correlating 
DMRs can be calculated.

## Example data

The example data comprise the raw reads of one sample and the already
called, but not filtered reads of that sample and further 7
samples. The samples blong to two groups, each of four samples. The
data are adopted from a recent lymphoma publication
([link](http://www.nature.com/ng/journal/v47/n11/full/ng.3413.html)). The
unmapped sample consists of two sequencing runs. These reads can be
mapped to a reduced genome and merged prior to methylation calling. In
addition to the raw and called methylation data, a reduced reference
genome, some gene annotations and gene expression data are
provided. This will enable you, to run the entire toolkit on a
small example region.

For more information please go to the [website](http://www.bioinf.uni-leipzig.de/Software/BAT).
