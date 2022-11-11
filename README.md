![This is an image](/images/srna_title.png)


# INTRODUCTION

This is a general pipeline for analysis of small RNA (sRNA) sequencing data from Illumina. It has mainly been developed with focus on the *Arabidopsis* TAIR10 genome and some related species, but it should be fully functional with other organisms. The scripts do not require to pass command-line arguments; settings like input data and reference genomic files have to be specified by editing the script in a text editor. Therefore, some very basic knowledge of linux and R is required. The scripts are generously commented in hashes (#) with complementary suggestions and hints (worth to read as a complement to this instruction). It produces basic background output for customised downstream analysis.
 
# SUPPORTED PLATFORMS

Linux, Mac OS

Shell scripts (*.sh) of this software were developed and tested using GNU bash (v4.4.20) in a Ubuntu 18.04 linux system. R scripts were developed using the R console (v4.1.1) in macOS Monterey.

# DEPENDENCIES

The following tools need to be installed and ideally available in the PATH environment. The pipeline is fully functional and tested with the following versions of the packages listed below. Other versions are very likely functional as well, but a detailed compatibility review of older and newer versions has not been done here. 

fastqc (v0.11.9)

cutadapt (v4.1 with Python 3.8.5)

bowtie (v1.3.1)

ShortStack (v3.8.5)

bedtools (v2.26.0)

samtools (v1.3.1)

# SETTING UP THE WORKING DIRECTORY AND THE GENOMIC REFERENCE FILES

The raw data (Illumina single-end fastq files) should be allocated in sample folders under a parent directory (/**REPLICATES_TOTAL**) following the file structure below. Single replicates are the basic units of this pipeline. Merging of replicates within conditions is usually considered after a first evaluation of the results, the merging itself can be performed at different stages (from raw sequences to fully processed files) depending on data structure, experimental design and taste. It is not handled here. Intermediary and final output files will be generated in respective sample folders. In the following example four samples are used: two conditions (*mutant* and *wt*) with two replicates each (*rep1* and *rep2*). 

```
    ├── REPLICATES_TOTAL
        └── mutant_rep1
        │   └── mutant_Rep1_FKDL210177288-1a-3_H3TM7DSX2_L4_1.fq.gz
        └── mutant_rep2
        │   └── mutant_Rep2_FKDL210177288-1a-4_H3TM7DSX2_L4_1.fq.gz
        └── wt_rep1
        │   └── Col_Rep1_FKDL210177288-1a-1_H3TM7DSX2_L4_1.fq.gz
        └── wt_rep2
            └── Col_Rep2_FKDL210177288-1a-2_H3TM7DSX2_L4_1.fq.gz
```

The paths to the working and sample directories, and reference genomic files have to be specified by editing the header of the respective shell script (*.sh).

```
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
```
If removing structural RNAs is wanted before downstream processing, a fasta file with the selected structural RNAs (including e.g. pre-tRNA, snoRNA, snRNA, rRNA) has to be prepared (with e.g. bedtools getfasta) and bowtie-indexed accordingly. Otherwise the pipeline can easily be modified skipping this step.

Annotation files (gtf or gff3) have to be transformed to a 6-column bed format looking like in the one below (**TAIR10_genes_isoform1_annot.bed**).

```
Chr1	3631	5899	AT1G01010	.	+
Chr1	5928	8737	AT1G01020	.	-
Chr1	11649	13714	AT1G01030	.	-
Chr1	23146	31227	AT1G01040	.	+
Chr1	28500	28706	AT1G01046	.	+
```
Two annotation bed files, one for genes and one for transposable elements (TE), are to be used here. To prepare bed files out of gtf or gff3 files is not straightforward. The [gff2bed](https://bedops.readthedocs.io/en/latest/content/reference/file-management/conversion/gff2bed.html) tool from BEDOPS suit is an option. Another possibility, often more pragmatic, is to process it with a combination of linux regular expressions and/or manual editing in a text editor.

Chromosome sizes should also be specified in a reference file (**TAIR10.chrom.sizes**) in the followin way.
```
Chr1	30427671
Chr2	19698289
Chr3	23459830
Chr4	18585056
Chr5	26975502
```

Annotation files in bed format, a fasta file with structural RNAs, and a text chromosome size file for the TAIR10 **Arabidopsis** genome are provided under [example](/example/genomic_reference_files).

Ordinary assembly fasta files and bowtie indexed fasta files (bowtie-build) should be available for the relevant genome and are usually downloadable from general or organism-specific genome repositories like ([TAIR10](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release)).

# INSTALLATION

Shell scripts can be cloned and run directly on a linux server.

```
git clone https://github.com/Ingefast/small_RNA_analyser.git
cd small_RNA_analyser
```

# WORKFLOW

## 1. Quality Control and Adapter Trimming.

Read quality (fastqc) is assessed before and after adapter trimming (cutadapt). Additionaly a size range trimming is performed, e.g. selecting a typical population of reads between 18 and 30 nt long.

Usage:
```
nohup bash sRNA.QC_Trimmer.sh
```
A final file with trimmed sequences in raw text format is generated (**sample.trimmed.txt**).

## 2. Removal of structural RNAs
Trimmed sequences belonging to structural RNAs and 'authentic' sRNA reads are segregated by mapping (bowtie) the reads to a fasta file with structural RNAs. Unmapped reads are considered for downstream analysis and, after a second mapping, only those aligning with zero  mismatches to the genome are kept (**sample.perfect.txt**).

Usage:
```
nohup bash sRNA.Filter.sh
```
In addition a table containing unique sequence signatures, their abundance and length is generated (**sample.perfect.pivot.table.txt**).

```
TCCGCTGTAGCACACAGGC 173995	19
GCGGACTGCTCGAGCTGC 136325	18
TCCGCTGTAGCACTTCAGGC 112283	20
CAGCGGACTGCTCGAGCTGC 103107	20
TCCGCTGTAGCACTTCAGGCT 51141	21
```

## 3. Mapping of selected sRNA size ranges.

Once the desired population of sRNA reads have been filtered, further processing can be done focusing on customised sRNA sizes. Apart of  the whole range of sRNA reads define  here (18-30nt), other usual choices are the 21-22 nt and 24 nt categories. The latter can be taylored in the **sRNA.Mapper.range_1830nt.sh** script by replacing the text string 'size_1830nt' to e.g. and 'size_2122nt' and changing the cutadapt size trimming settings as:
```
cutadapt --no-trim -m 21 -M 22 -o sample.size_sel.fastq ../sample.perfect.fastq > summary_NOtrimming.txt;
```
At this step time reads are mapped using ShortStack which is a wrapper of bowtie designed to handle multimapping reads.

Usage:
```
nohup bash sRNA.Mapper.range_1830nt.sh
```

Important output files are the following:

1. A table of genomic features and their expression values in Read Per Million (RPM).

**genomic_RPM_values.genes.txt**, **genomic_RPM_values.tes.txt**
```
Chr1    23146   31227   AT1G01040       .       +       1.245824286
Chr1    28500   28706   AT1G01046       .       +       0
Chr1    31170   33153   AT1G01050       .       -       0.5813833333
Chr1    33666   37840   AT1G01060       .       -       0.7474928571
```

2. Similar tables but with raw count values.

**gene.counts.txt**, **te.counts.txt**

Combining the two kind of tables across sample above one can easily build a general table for the whole experiment like **total_table.size_24.genes.txt** in the [output](/example/output) folder. Minimal formatting of this table allows for differential expression analysis in e.g. [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html).

3. Table with sRNA read mapping coordinates, unique read identity number (fourth column), and feature affected.

**gene.reads.txt**, **te.reads.txt**
```
Chr1	16990	17013	776187	Chr1	16882	17009	AT1TE00020
Chr1	17004	17027	84575	Chr1	16882	17009	AT1TE00020
Chr1	17004	17027	84578	Chr1	17023	18924	AT1TE00025
Chr1	17856	17874	178759	Chr1	17023	18924	AT1TE00025
```

4. Table with a unique read identity number and their sequence linking to the tables above.

**srna_ids.txt**
```
520486  GATTCAAGACTTCTCGGTACT
84989   GTACTGCAAAGTTCTTCCGCC
155073  ACTGCAAAGTTCTTCCGCCTGAT
231101  ACTGCAAAGTTCTTCCGCCT
```
5. A browsable bedGraph file **sample.bedGraph** (see Figure 1A).

# DOWNSTREAM ANALYSES

The processed data provided so far opens a wide range of analytical possibilities such as studying the relationship between sRNA expression and occurrence of diverse epigenetic marks and/or expression and silencing of genes,transposable elements.

Two common basic analyses are presented here as examples.

## 1. Checking data structure and replicability in the data set. Multivariate analysis and correlogram.
With **sRNA.correlogram_plotter.r** it is possible to analyse the table **total_table.size_24.genes.txt** to plot the correlation between particular sRNA expression values across all the samples. This is a good way to check for deviant samples and assess replicability. A scatterplot in the lower diagonal panel is presented and Pearson correlation coefficients in the upper panel (Figure 2B). The script **sRNA.ordination_analyser.r** performs a multivariate analysis using the same  input table to evaluate the similarities between samples. By default the samples are ordinated with Nonmetric multidimensional scaling (NMDS) as implemented in the R library vegan (Figure 2C). Principal component analysis (PCA) and redundancy analysis (RDA) are also available as alternative ordination approaches.

![This is an image](/images/figure1.png)

*Figure 1*. (A) bedGraph files of a wildtype in *Capsella* for different sRNA sizes. (B) Correlogram of 24nt sRNA values over genes in three conditions with two replicates each. (C) NMDS diagram of the same dataset.

## 2. sRNA size distribution over over genes and transposable elements.
Understanding the role of sRNA of particular sizes on the expression of genes and TEs is central in any sRNA study. The script **sRNA.size_distribution_plotter.r** plots the abundance (RPM) of sRNA reads of different sizes over selected genomic features (genes and TEs). Inputs are the previously generated **gene.reads.txt** and **te.reads.txt** files. Additionally, a file with total number of mapped reads for each sample has to be created manually (**read_n_baseline.txt**) in order to establish a baseline for normalisation, 

```
sample_name	mapped
mutantA_rep1	646980
mutantA_rep2	991033
mutantB_rep1	491573
mutantB_rep2	546385
wt_rep1	390199
wt_rep2	999514
```
The number of mapped reads can be estimated with e.g.
```
wc -l sample.perfect.txt
```

![This is an image](/images/figure2.png)

*Figure 2*. sRNA size distribution over genes and transposable elements in *Arabidopsis*.


 
# REFERENCES

Modified versions of this pipeline have been used to process the sRNA datasets in the following papers:

1. Martinez G et al (2018). Paternal easiRNAs regulate parental genome dosage in *Arabidopsis*. **Nature Genetics** 50 (2) 193-198.

2. Wang Z et al (2020). Polymerase IV Plays a Crucial Role in Pollen Development in *Capsella*. **Plant Cell** 32 (4) 950-966.

# CONTACT
juan.santos at slu.se
