#!/bin/bash

#Juan Santos, 2021 August

#this defines the samples to be analyzed
sample_list="endo_pol_x_pol_rep1  endo_pol_x_pol_rep2  endo_pol_x_wt_rep1  endo_pol_x_wt_rep2  endo_wt_x_wt_rep1  endo_wt_x_wt_rep2  pollen_pol_x_pol_rep1  pollen_pol_x_pol_rep2  pollen_pol_x_wt_rep1  pollen_pol_x_wt_rep2  pollen_wt_x_wt_rep1  pollen_wt_x_wt_rep2";





####### Step 0. Concatenating.


for sample in $sample_list  ; do

	cd /media/diskd/jiali_sRNA_grafting/REPLICATES_TOTAL/${sample};

	#cat *fq > sample.fastq;
		
done
wait




####### Step 1. Quality Control of raw files

for sample in $sample_list  ; do

	cd /media/diskd/jiali_sRNA_grafting/REPLICATES_TOTAL/${sample};

		mkdir QC_raw_1;
		mkdir QC_raw_2;

		fastqc *_1.fq -o QC_raw_1&
		fastqc *_2.fq -o QC_raw_2&
		
done
wait



####### Step 2a. Quality trimming of raw files

for sample in   $sample_list  ; do

	cd /media/diskd/jiali_sRNA_grafting/REPLICATES_TOTAL/${sample};

#adapters can be found with minion search-adapter -i *.fastq -show 3

#This does a first trimming
# More info on this adapaters can be found in https://www.neb.com/~/media/Catalog/All-Products/FAC109E8FD1341339AEADA0A081814C7/Datacards%20or%20Manuals/manualE7300.pdf
#OR
#https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html


	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 5 -m 18 -M 30 -o 3trimmed_1.fastq *_1.fq > summary_3trimming_1.txt;
	
	cutadapt -b GATCGTCGGACTGTAGAACTCTGAACGTGTAGATCTCGGTGGTCGCCGTATCATT -O 5 -m 18 -M 30 -o 3trimmed_2.fastq *_2.fq > summary_3trimming_2.txt;

	cat 3trimmed_1.fastq | awk 'NR%4==2{print}' > trimmed_1.txt;
	cat 3trimmed_2.fastq | awk 'NR%4==2{print}' > trimmed_2.txt;

	#gzip ${sample}*.fastq;

done
wait





####### Step 3. Quality Control of trimmed seqs

for sample in  $sample_list; do

	cd /media/diskd/jiali_sRNA_grafting/REPLICATES_TOTAL/${sample};

	mkdir QC_trimmed_1;

	fastqc 3trimmed_1.fastq -o QC_trimmed_1&

	mkdir QC_trimmed_2;

	fastqc 3trimmed_2.fastq -o QC_trimmed_2&

done
wait


