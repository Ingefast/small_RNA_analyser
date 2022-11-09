#!/bin/bash

#Juan Santos, 2022 November


##################################################################################
###########  Reference files to be used and working dir        ###################
##################################################################################

#genomic assembly in fasta format
genome_assembly_fasta="/home/jsantos/genome_ref/TAIR10/bowtie1/TAIR10_chr_all.fasta"

#genomic chromosome sizes txt file
genome_chr_sizes="/home/jsantos/genome_ref/TAIR10/TAIR10.chrom.sizes"

#bowtie1 indexes to be used
genome_non_sRNA_bw_idx="/home/jsantos/genome_ref/TAIR10/sRNA_filter/non_sRNA"
genome_bw_idx="/home/jsantos/genome_ref/TAIR10/bowtie1/TAIR10"

#annotation files in bed format
genes_annot_bed="/home/jsantos/genome_ref/TAIR10/bed_landscape/TAIR10_genes_isoform1_annot.bed"
tes_annot_bed="/home/jsantos/genome_ref/TAIR10/bed_landscape/TAIR_transposable_elements_6col.bed"

##################################################################################

#working directory where all the sample directories are allocated
working_dir="/media/diskb/nicolas_tmp/REPLICATES_TOTAL"

#this defines the samples to be analysed (dir names)
sample_list="mutant_rep1  mutant_rep2  wt_rep1  wt_rep2";
#sample_list="test_rep1";

##################################################################################



### Step 1 . Selecting sRNA size_1830 range and mapping

for sample in $sample_list  ; do

	cd $working_dir/${sample};
	
	mkdir size_1830; cd size_1830;

	#this cuts without adapter in desired size
	cutadapt --no-trim -m 18 -M 30 -o sample.size_sel.fastq ../sample.perfect.fastq > summary_NOtrimming.txt;

	#Mapping with protocol for handling multi-mapped reads is f (fractional-seeded guide)
	ShortStack --align_only --bowtie_cores 20 --mismatches 0  --nohp --mmap f --outdir shortstack --genome $genome_assembly_fasta --readfile sample.size_sel.fastq;

done
wait



### Step 2 . Calculating coverage and normalizing it to RPM.

for sample in $sample_list  ; do

	cd $working_dir/${sample}/size_1830/shortstack;
	
	samtools sort -@ 20 sample.size_sel.bam -o sample.size_sel.sorted.bam;
	samtools index sample.size_sel.sorted.bam;

	#transforms bam to bed
	bamToBed -i sample.size_sel.sorted.bam > sample.size_sel.sorted.bed;

	#calculation normalization factor, reads per million
	librarySize=$(wc -l ../../sample.perfect.txt | awk '{total+=$1}END{print total}');
	b=1000000;
	factor=`echo "$b / $librarySize "|bc -l`;
	
	#Creates a normalised coverage file in bed format (4 columns bed)
	bedtools genomecov -i sample.size_sel.sorted.bed -g $genome_chr_sizes -bga -scale $factor | grep -v "mitochondria\|chloroplast" > sample.coverage.bed;

done
wait



### Step 3 . Calculating diverse final output; bedGraph, RPM values on genomic features and more

for sample in $sample_list  ; do

	cd $working_dir/${sample}/size_1830/shortstack;
	
	#Creates a bedGraph file with normalized values
	bedtools genomecov -ibam sample.size_sel.sorted.bam -bg -scale $factor > sample.bedGraph&

	#Creates a coverage file for differential expression of genes with read counts (not normalised) per gene
	bedtools multicov -D -bams sample.size_sel.sorted.bam -bed  $genes_annot_bed> gene.counts.txt&

	#Creates a coverage file for differential expression of TEs with read counts (not normalised) per TE
	bedtools multicov -D -bams sample.size_sel.sorted.bam -bed  $tes_annot_bed> te.counts.txt&
	
	#Calculates genomic RPM values for genes
	bedtools map -a $genes_annot_bed -b sample.coverage.bed -c 4 -o mean>genomic_RPM_values.genes.txt&

	#Calculates genomic RPM values for tes
	bedtools map -a $tes_annot_bed -b sample.coverage.bed -c 4 -o mean>genomic_RPM_values.tes.txt&
	
	#making list of reads id with sequences (two columns, id after sample.size_sel.fastq)
	samtools view sample.size_sel.sorted.bam| cut -f 1,10 >srna_ids.txt&

	#extracting reads mapping into genes (column 4 is read id after sample.size_sel.fastq)
	intersectBed -abam sample.size_sel.sorted.bam -b $genes_annot_bed -wa -wb -bed|cut -f 1,2,3,4,13,14,15,16> gene.reads.txt&

	#extracting reads mapping into TEs  (column 4 is read id after sample.size_sel.fastq)
	intersectBed -abam sample.size_sel.sorted.bam -b $tes_annot_bed -wa -wb -bed|cut -f 1,2,3,4,13,14,15,16> te.reads.txt&

done
wait




