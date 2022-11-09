#!/bin/bash

#Juan Santos, 2022 November

##################################################################################
###########  Reference files to be used and working dir        ###################
##################################################################################

#bowtie1 indexes to be used
genome_non_sRNA_bw_idx="/home/jsantos/genome_ref/TAIR10/sRNA_filter/non_sRNA"
genome_bw_idx="/home/jsantos/genome_ref/TAIR10/bowtie1/TAIR10"


##################################################################################

#working directory
working_dir="/media/diskb/nicolas_tmp/REPLICATES_TOTAL"

#this defines the samples to be analysed (dir names)
sample_list="mutant_rep1  mutant_rep2  wt_rep1  wt_rep2";

##################################################################################


#Step 1. Mapping against structural RNA (non-sRNA sequences), discarding those, and keeping (filtering) the reads which are considered 'real' sRNAs

for sample in $sample_list  ; do

	cd $working_dir/${sample};

	bowtie -p 20 -v 1 -k 1 -a -r --chunkmbs 200 -x $genome_non_sRNA_bw_idx  sample.trimmed.txt -S sample.sam 2>&1 | tee bowtie_filter_report.txt;
	samtools sort -@ 20 sample.sam -o sample.sorted.bam;
	samtools index sample.sorted.bam;
	
	samtools view -b -f 4 sample.sorted.bam > sample.filtered.bam; #the unmapped are the ones to be kept
	samtools view -b -F 4 sample.sorted.bam > sample.non_sRNA.bam; #the mapped ones are structural RNAs - not interesting.

	samtools fastq sample.filtered.bam > sample.filtered.fastq; #
	cat sample.filtered.fastq | awk 'NR%4==2{print}' > sample.filtered.txt;
	rm sample.sam;

done
wait



#Step 2: Generating an alignment of sRNA mapping perfectly to the genome

for sample in $sample_list  ; do

	cd $working_dir/${sample};
	
	bowtie -p 20 -v 0 --best -r -x $genome_bw_idx  sample.filtered.txt -S sample.tmp.sam 2>&1 | tee bowtie_all_sRNA_report.txt;
	samtools view -b -F 4 sample.tmp.sam > sample.perfect.bam; #extract only mapped
	samtools sort -@ 20 sample.perfect.bam -o sample.perfect.sorted.bam;
	samtools index sample.perfect.sorted.bam;
	samtools fastq sample.perfect.sorted.bam > sample.perfect.fastq; #sRNA file
	
	cat sample.perfect.fastq | awk 'NR%4==2{print}' > sample.perfect.txt;
	rm sample.tmp* sample.perfect.bam;

done
wait



#Step 3: This makes a pivot table with three columns: unique sequences, abundance and lenght

for sample in $sample_list  ; do

	cd $working_dir/${sample};

	#makes a pivot table for all the different sRNA and sorts the output with the most frequent sRNAs on top

	sort sample.perfect.txt | uniq -c | sort -nr>tmp_sample.summary;

	awk '{ print $2 " " $1}' tmp_sample.summary>tmp_sample.swapped.summary;#this swaps columns

	cut -f 1 -d ' ' tmp_sample.swapped.summary | awk '{ print length }'>tmp_sample.length;#this makes a tmp file counting the length of the sRNA

	paste tmp_sample.swapped.summary tmp_sample.length >sample.perfect.pivot.table.txt;#this pastes the two previous files

	rm tmp_sample*;

done
wait
