#!/bin/bash

#Juan Santos, 2022 November

##################################################################################

#working directory where all the sample directories are allocated
working_dir="/media/diskb/nicolas_tmp/REPLICATES_TOTAL"

#this defines the samples to be analysed (dir names)
sample_list="mutant_rep1  mutant_rep2  wt_rep1  wt_rep2";
#sample_list="test_rep1";

##################################################################################


####### Step 1. Quality Control of raw files

for sample in $sample_list  ; do

	cd $working_dir/${sample};
	
	mkdir QC_raw_1;
	gunzip *_1.fq.gz;
	fastqc --outdir QC_raw_1  *_1.fq;

done
wait



####### Step 2. Adapter and size trimming of raw files

for sample in $sample_list  ; do

	cd $working_dir/${sample};

#if used adapters are unknown can be found with minion e.g. minion search-adapter -i *.fq -show 3
#hg clone https://galaxy-ntnu.bioinfo.no/toolshed_nels/repos/kjetil/package_minion_15_065

	cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -O 5 -m 18 -M 30 -o 3trimmed_1.fastq *_1.fq > summary_3trimming_1.txt; #performs adapter trimming and selects read population in the desired size range
	
	cat 3trimmed_1.fastq | awk 'NR%4==2{print}' > sample.trimmed.txt; #creates a raw txt file with size trimmed sRNAs for working downstream

done
wait



####### Step 3. Quality Control of trimmed seqs

for sample in $sample_list  ; do

	cd $working_dir/${sample};
	
	mkdir QC_trimmed_1;
	fastqc --outdir QC_trimmed_1 3trimmed_1.fastq ;

done
wait
