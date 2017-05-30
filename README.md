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
postprocessings. For mapping, [BAT_mapping](http://www.bioinf.uni-leipzig.de/Software/BAT/mapping.md#mapping)
employs
[segemehl](http://www.bioinf.uni-leipzig.de/Software/segemehl/), a
performant and highly sensitive short-read aligner with a specialized
bisulfite mode
([link](http://bioinformatics.oxfordjournals.org/content/28/13/1698.long)),
but due to the modularity of BAT, this step could be exchanged by
running a different bisulfite-sensitive aligner. In addition to a
bisulfite-specific quality filtering, aligned reads are converted to
an indexed and sorted
[BAM](https://samtools.github.io/hts-specs/SAMv1.pdf) file. Basic
mapping statistics such as the number of mapped pairs/reads, the
number of reads with a single (unique-mapped) or multiple alternative
alignments (multi-mapped), the distribution of the multiplicity of
read alignment, and the distribution of the edit distance of read
alignments are calculated by
[BAT_mapping_stat](http://www.bioinf.uni-leipzig.de/Software/BAT/mapping.md#statistics). If there are multiple
datasets per sample (e.g., due to multiplexing on different lanes),
all alignment files corresponding to one sample can be merged using
[BAT_merging](http://www.bioinf.uni-leipzig.de/Software/BAT/mapping.md#merging). It also enables the addition of
dataset-specific read group information during the merging process.

[**Calling**](http://www.bioinf.uni-leipzig.de/Software/BAT/calling)

Following mapping, the methylation information needs to be extracted
from the alignments, referred to as methylation calling. First,
[BAT_calling](http://www.bioinf.uni-leipzig.de/Software/BAT/calling.md#calling) takes the alignments and generates
a VCF file that contains information for each cytosine including the
sequence context, coverage, detailed number of covering nucleotides,
and the estimated methylation rate.  Second, cytosine positions can be
filtered by coverage, genomic context, and methylation rate, using
[BAT_filtering](http://www.bioinf.uni-leipzig.de/Software/BAT/calling.md#filtering). The output is again in VCF
format but it is also provided as
[bedGraph](http://genome.ucsc.edu/FAQ/FAQformat.html#format1.8) file
with the estimated methylation rate in the fourth column, ready for
loading in [IGV](https://www.broadinstitute.org/igv/) or uploading to
the [UCSC genome browser](https://genome.ucsc.edu). Furthermore, the
coverage and methylation rate distributions for all and filtered
positions are illustrated as barplots.

[**Analysis**](http://www.bioinf.uni-leipzig.de/Software/BAT/analysis)

The third module covers the basic analysis of two
groups of a single sample or up to multiple samples. At first, various
helpful summary,
[bedGraph](http://genome.ucsc.edu/FAQ/FAQformat.html#format1.8) and
[bigWig](https://genome.ucsc.edu/goldenpath/help/bigWig.html) files
for all samples are created with
[BAT_summarize](http://www.bioinf.uni-leipzig.de/Software/BAT/analysis.md#summarize). Furthermore, a
[Circos](http://circos.ca) plot containing a methylation rate heatmap
for each sample could be provided. Overview plots comprising
hierarchical clustering, genome-wide average methylation rate
boxplots, correlation plots of mean group methylation rates per
position and distribution of position-wise group differences are
plotted.  Specific regions of interest or annotations, e.g.,
transcription factor binding sites (TFBS), CpG islands, promoter
regions, [BAT_annotation](http://www.bioinf.uni-leipzig.de/Software/BAT/analysis.md#annotation) can be used to get
an insight into the methylation of the samples in those regions. Basic
statistics, like length of annotation items (in nucleotides and Cs),
are calculated. In addition, the distribution of the average
group-wise methylation rates per annotation item, clustering heatmaps
containing all samples, and boxplots of the single sample average
methylation rates per annotation item are shown.

[**DMRs**](http://www.bioinf.uni-leipzig.de/Software/BAT/dmrs)

Finally, the calling analysis of DMRs is coverd by
[BAT_DMRcalling](http://www.bioinf.uni-leipzig.de/Software/BAT/dmrs.md#calling) and
[BAT_correlating](http://www.bioinf.uni-leipzig.de/Software/BAT/dmrs.md#correlating). The DMR calling tool
[metilene](http://www.bioinf.uni-leipzig.de/Software/metilene/)
identifies DMRs between two groups from one or more samples very
quickly and accurately. Subsequently, the raw metilene output can be
filtered and converted to BED-like or bedGraph format. Basic DMR
statistics including length distributions (in nucleotides and Cs),
distribution of group methylation differences, and scatterplots of
methylation means of group 1 vs. group 2 as well as methylation
difference vs. q-value of DMRs is illustrated. Finally,
[BAT_correlating](http://www.bioinf.uni-leipzig.de/Software/BAT/dmrs.md#correlating) facilitates the identification
of correlating DMRs (cDMRs), i.e., DMRs where the methylation change
correlates with a change in the expression of the associated
genes. However, it is not restricted to DMRs as input but can also be
used for inspecting other annotation items such as promoter regions or
TFBS. In result, linear and non-linear correlation effects are tested
and the results are reported as text file and correlation plot for
easy visual inspection.


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


## Obtaining

To download the BAT scripts, please use git to download the most recent development tree. Currently, the tree is hosted on github, and can be obtained via:

git clone https://github.com/helenebioinf/BAT


## Docker

If you prefer to not install all [dependencies](http://www.bioinf.uni-leipzig.de/Software/BAT/install), you can use the [BAT docker image](http://www.bioinf.uni-leipzig.de/Software/BAT/install.md/#docker). 
Dependencies and scripts are installed - simply pull the image. To test it, download 
the input data including the folder structure [here](http://www.bioinf.uni-leipzig.de/Software/BAT/BAT_example_structure.tar.gz) 
(985 MB), run the docker image and the run-script. For a quick start, have a look 
[here](http://www.bioinf.uni-leipzig.de/Software/BAT/install.md/#docker)


## Detailed Information
For more information please go to the [website](http://www.bioinf.uni-leipzig.de/Software/BAT).

## Contact

Please report any issues or questions to helene [at] bioinf [dot] uni-leipzig.de
