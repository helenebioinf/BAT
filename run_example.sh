PATH=$PATH:$(pwd)

###########
# mapping #
###########
echo "Mapping Module"
# mapping
BAT_mapping -g genomes/hg19/hg19_chr19.fa -q raw/S5.1_R1.fastq.gz -p raw/S5.1_R2.fastq.gz -i genomes/hg19/hg19_chr19 -o mapped/S5.1 -t 4
BAT_mapping -g genomes/hg19/hg19_chr19.fa -q raw/S5.2_R1.fastq.gz -p raw/S5.2_R2.fastq.gz -i genomes/hg19/hg19_chr19 -o mapped/S5.2 -t 4

# mapping statistics
BAT_mapping_stat --bam mapped/S5.1.bam --excluded mapped/S5.1.excluded.bam --fastq raw/S5.1_R1.fastq.gz >mapped/S5.1.stat
BAT_mapping_stat --bam mapped/S5.2.bam --excluded mapped/S5.2.excluded.bam --fastq raw/S5.2_R1.fastq.gz >mapped/S5.2.stat

# merging
BAT_merging -o mapped/S5.bam --bam mapped/S5.1.bam,mapped/S5.2.bam

###########
# calling #
###########
echo "Calling Module"
# calling
BAT_calling -d genomes/hg19/hg19_chr19.fa -q mapped/S5.bam -o called/

# filter_vcf
for i in called/S[1-8].vcf.gz; do
	o=`echo $i | sed 's/.vcf.gz/_CG.vcf.gz/'`;
	BAT_filter_vcf --vcf $i --out $o --context CG --MDP_min 10 --MDP_max 100;
done


############
# analysis #
############
echo "Analysis Module"
# summarize
BAT_summarize --in1 called/S1_CG.bedgraph,called/S2_CG.bedgraph,called/S3_CG.bedgraph,called/S4_CG.bedgraph --in2 called/S5_CG.bedgraph,called/S6_CG.bedgraph,called/S7_CG.bedgraph,called/S8_CG.bedgraph --groups control,case --h1 S1,S2,S3,S4 --h2 S5,S6,S7,S8 --out data/example  --cs genomes/hg19/hg19_chr19.chrom.sizes

# overview
BAT_overview.R  -i data/example_summary_control_case.bedgraph -o data/example_overview --groups control,case

# annotation
BAT_annotation -b genomes/hg19/TFBS.bed -i data/example_summary_control_case.bedgraph --groups control,case -o annotation/example_TFBS


########
# DMRs #
########
echo "DMR Module"
# DMR calling
BAT_DMRcalling -q data/example_metilene_control_case.txt -o DMRs/metilene -a control -b case

# cDMRs
bedtools intersect -wa -wb -a DMRs/metilene_qval.0.05.bed -b genomes/hg19/gene.bed | cut -f1-3,9 >DMRs/DMR_gene.txt
ls data/*.bw | grep -v "mean" | grep -v "diff" | sed 's/\t/\n/' >DMRs/methylation_files.list
ls expression/* | sed 's/\t/\n/' >DMRs/expression_files.list
echo -e 'S1\tcontrol\nS2\tcontrol\nS3\tcontrol\nS4\tcontrol\nS5\tcase\nS6\tcase\nS7\tcase\nS8\tcase' >DMRs/sample_to_group.txt

BAT_correlating -b DMRs/DMR_gene.txt -e DMRs/expression_files.list -m DMRs/methylation_files.list -g DMRs/sample_to_group.txt -i control,case -o DMRs/correlation
