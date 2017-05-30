echo "cleaning up mapped/"
rm mapped/S5.1.bam mapped/S5.1.bam.bai mapped/S5.1.excluded.bam mapped/S5.1.excluded.bam.bai mapped/S5.1.mapping.log mapped/S5.1.pdf mapped/S5.1.stat mapped/S5.1.unmapped.gz
rm mapped/S5.2.bam mapped/S5.2.bam.bai mapped/S5.2.excluded.bam mapped/S5.2.excluded.bam.bai mapped/S5.2.mapping.log mapped/S5.2.pdf mapped/S5.2.stat  mapped/S5.2.unmapped.gz
rm mapped/S5.bam mapped/S5.sam.gz mapped/S5.bam.bai mapped/S5.sam.idx

echo "cleaning up called/"
rm called/S5.calling.log
rm called/S1_CG.bedgraph called/S2_CG.bedgraph called/S3_CG.bedgraph called/S4_CG.bedgraph called/S5_CG.bedgraph called/S6_CG.bedgraph called/S7_CG.bedgraph called/S8_CG.bedgraph
rm called/S1_CG.pdf called/S2_CG.pdf called/S3_CG.pdf called/S4_CG.pdf called/S5_CG.pdf called/S6_CG.pdf called/S7_CG.pdf called/S8_CG.pdf
rm called/S1_CG.vcf.gz called/S2_CG.vcf.gz called/S3_CG.vcf.gz called/S4_CG.vcf.gz called/S5_CG.vcf.gz called/S6_CG.vcf.gz called/S7_CG.vcf.gz called/S8_CG.vcf.gz

echo "cleaning up data/"
rm data/example_diff_control_case.bedgraph data/example_diff_control_case.bw data/example_mean_case.bedgraph data/example_mean_case.bw data/example_mean_control.bedgraph data/example_mean_control.bw data/example_metilene_control_case.txt data/example_summary_control_case.bedgraph
rm data/example_S1.bedgraph data/example_S1.bw data/example_S2.bedgraph data/example_S2.bw data/example_S3.bedgraph data/example_S3.bw data/example_S4.bedgraph data/example_S4.bw data/example_S5.bedgraph data/example_S5.bw data/example_S6.bedgraph data/example_S6.bw data/example_S7.bedgraph data/example_S7.bw data/example_S8.bedgraph data/example_S8.bw
rm data/example_overview.pdf 

echo "cleaning up annotation/"
rm annotation/example_TFBS.pdf annotation/example_TFBS.txt

echo "cleaning up DMRs/"
rm DMRs/metilene.log DMRs/metilene.out DMRs/metilene_qval.0.05.bed DMRs/metilene_qval.0.05.bedgraph DMRs/metilene_qval.0.05.out DMRs/metilene_qval.0.05.pdf
rm DMRs/DMR_gene.txt DMRs/expression_files.list DMRs/methylation_files.list DMRs/sample_to_group.txt
rm DMRs/correlation.txt DMRs/correlation_chr19:10624503-10624742_ENSG00000180739.pdf DMRs/correlation_chr19:11077326-11077585_ENSG00000127616.pdf DMRs/correlation_chr19:11467267-11467765_ENSG00000105520.pdf DMRs/correlation_chr19:10625238-10625388_ENSG00000180739.pdf DMRs/correlation_chr19:11081481-11081743_ENSG00000127616.pdf DMRs/correlation_chr19:11516916-11517167_ENSG00000205517.pdf DMRs/correlation_chr19:10625393-10625445_ENSG00000180739.pdf DMRs/correlation_chr19:11084338-11084493_ENSG00000127616.pdf DMRs/correlation_chr19:11532291-11532486_ENSG00000198003.pdf DMRs/correlation_chr19:10793898-10794095_ENSG00000129351.pdf DMRs/correlation_chr19:11084524-11084841_ENSG00000127616.pdf DMRs/correlation_chr19:11557459-11557835_ENSG00000130175.pdf DMRs/correlation_chr19:10858170-10858760_ENSG00000079805.pdf DMRs/correlation_chr19:11086226-11086876_ENSG00000127616.pdf DMRs/correlation_chr19:11574129-11574509_ENSG00000196361.pdf DMRs/correlation_chr19:10928642-10928813_ENSG00000079805.pdf DMRs/correlation_chr19:11087248-11088068_ENSG00000127616.pdf DMRs/correlation_chr19:11590469-11591012_ENSG00000196361.pdf DMRs/correlation_chr19:10947344-10947483_ENSG00000099203.pdf DMRs/correlation_chr19:11090548-11091241_ENSG00000127616.pdf DMRs/correlation_chr19:11590469-11591012_ENSG00000267477.pdf DMRs/correlation_chr19:10947344-10947483_ENSG00000214212.pdf DMRs/correlation_chr19:11094281-11094942_ENSG00000127616.pdf DMRs/correlation_chr19:11591538-11591910_ENSG00000196361.pdf DMRs/correlation_chr19:10997976-10998526_ENSG00000142453.pdf DMRs/correlation_chr19:11142536-11143257_ENSG00000127616.pdf DMRs/correlation_chr19:11591538-11591910_ENSG00000267477.pdf DMRs/correlation_chr19:11006467-11007609_ENSG00000142453.pdf DMRs/correlation_chr19:11144943-11145419_ENSG00000127616.pdf DMRs/correlation_chr19:11605792-11605902_ENSG00000161914.pdf DMRs/correlation_chr19:11040212-11040659_ENSG00000130733.pdf DMRs/correlation_chr19:11159207-11159488_ENSG00000127616.pdf DMRs/correlation_chr19:11605792-11605902_ENSG00000267477.pdf DMRs/correlation_chr19:11040212-11040659_ENSG00000142444.pdf DMRs/correlation_chr19:11160325-11161217_ENSG00000127616.pdf DMRs/correlation_chr19:11072934-11073838_ENSG00000127616.pdf DMRs/correlation_chr19:11301235-11301809_ENSG00000197256.pdf 
