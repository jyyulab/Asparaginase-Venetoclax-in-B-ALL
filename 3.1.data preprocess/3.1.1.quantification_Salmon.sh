#!urs/bin/bash
#Coded by Xin Huang(xhuang@stjude.org)

#Reference genome: hg38 
#salmon v0.13.1
#picard v2.9.4
#java v1.8.0_60

for sample in sample1 sample2 sample3 ....
do
	java -jar ../path_to_picard/picard.jar SamToFastq I=../../path_to_bams/sample1.bam F=../../path_to_fastq/sample1_R1.fastq.gz F2=../../path_to_fastq/sample1_R2.fastq.gz VALIDATION_STRINGENCY=SILENT;salmon quant -i ../path_to_reference/salmon_usage/Salmon_index_hg38/ -l A -1 ../../path_to_fastq/sample1_R1.fastq.gz -2 ../../path_to_fastq/sample1_R2.fastq.gz -o ../../salmon/sample1_salmon
done


