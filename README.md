![This is an image](/images/srna_title.png)
### *Juan Santos*, September 2022


# INTRODUCTION

This is a general pipeline for analysis of small RNA (sRNA) sequencing data from Illumina. It has mainly been developed with focus on the *Arabidopsis* TAIR10 genome and some related species, but it should be fully functional with other organisms. The scripts do not require to pass command-line arguments; settings like input data and reference genomic files have to be specified by editing the script in a text editor. Therefore, some very basic knowledge of linux and R is required. The scripts are generously commented in hashes (#) with complementary suggestions and hints (worth to read as a complement to this instruction).
 
# SUPPORTED PLATFORMS

Linux, Mac OS

Shell scripts (*.sh) of this software were developed and tested using GNU bash (v4.4.20) in a Ubuntu 18.04 linux system. R scripts were developed using the R console (v4.1.1) in macOS Monterey.

# DEPENDENCIES

The following tools need to be installed and ideally available in the PATH environment. The pipeline is fully functional and tested with the following versions of the packages listed below. Other modern versions are very likely functional as well, but a detailed compatibility review of older and newer versions has not been done here. 

fastqc (v0.11.9)

cutadapt (v4.1 with Python 3.8.5)

bowtie (v1.3.1)

ShortStack (v3.8.5)

bedtools (v2.26.0)

samtools (v1.3.1)

# SETTING UP THE WORKING DIRECTORY AND THE GENOMIC REFERENCE FILES

The raw data (fastq files) are allocated in sample folders under a parent directory (/**REPLICATES_TOTAL**) following the file structure below. Single replicates are the basic units of this pipeline. Merging of replicates within conditions is usually considered after a first evaluation of the results, the merging itself can be performed at different stages (from raw sequences to fully processed files) depending on data structure, experimental design and taste. Most intermediary and final output files will be generated in respective sample folders. As example four samples are used: two conditions (*mutant* and *wt*) with two replicates each (*rep1* and *rep2*). 

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

The working directory, sample directories and reference genomic files have to be specified by editing the header of the respective shell script (*.sh).

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
#sample_list="test_rep1";

```

Ordinary fasta files and bowtie indexed fasta files (bowtie-build) should be available for the relevant assembly ([TAIR10](https://www.arabidopsis.org/download/index-auto.jsp?dir=%2Fdownload_files%2FGenes%2FTAIR10_genome_release))

If removing structural RNAs is wanted prior to downstream processing, a fasta file with the selected structural RNAs (including e.g. pre-tRNA, snoRNA, snRNA, rRNA) has to be prepared (bedtools getfasta) and bowtie-indexed accordingly. Otherwise the pipeline can be modified accordingly skipping this step.

Annotation files (gtf or gff3) have to be transformed to a 6-column bed format looking like in the one below (**TAIR10_genes_isoform1_annot.bed**).

```
Chr1	3631	5899	AT1G01010	.	+
Chr1	5928	8737	AT1G01020	.	-
Chr1	11649	13714	AT1G01030	.	-
Chr1	23146	31227	AT1G01040	.	+
Chr1	28500	28706	AT1G01046	.	+
Chr1	31170	33153	AT1G01050	.	-
Chr1	33666	37840	AT1G01060	.	-
Chr1	38752	40944	AT1G01070	.	-
Chr1	44677	44787	AT1G01073	.	+
Chr1	45296	47019	AT1G01080	.	-
```
Two annotation bed files, one for genes and one for transposable elements, are to be used here. To prepare bed files out of gtf or gff3 files is not straightforward. The gff2bed tool from BEDOPS suit is an option. Another possibility, sometimes more pragmatic, is processing it with a combination of linux regular expressions and/or manual editing in a text editor.

Chromosome sizes should also be specified in a reference file (**TAIR10.chrom.sizes**).
```
Chr1	30427671
Chr2	19698289
Chr3	23459830
Chr4	18585056
Chr5	26975502
```
# INSTALLATION

Shell scripts can be cloned and run directly on a linux server.

# WORKFLOW

## 1. Quality Control and Adapter Trimming.

Read quality (fastqc) is assessed before and after adapter trimming (cutadapt). An additional size range trimming is performed, typically selecting a population of reads between 18 and 30 nt long.

```
nohup bash sRNA.Filter.sh
```
A final file with trimmed reads in raw text format is generated (**sample.trimmed.txt**).

## 2. Removal of structural RNAs
Trimmed sequences belonging to structural RNAs and 'real' sRNA reads are separated by mapping (bowtie) the reads to a fasta file of only structural RNAs. Unmapped reads are considered for downstream analysis and, after a second mapping, only those aligning without mismatches to the genome are kept (**sample.perfect.txt**).

```
nohup bash sRNA.Filter.sh
```
In addition a table containing unique sequence signatures, their abundance and length is generated (**sample.perfect.pivot.table.txt**).

```
TCCGCTGTAGCACACAGGC 173995	19
GCGGACTGCTCGAGCTGC 136325	18
TCCGCTGTAGCACTTCAGGC 112283	20
CAGCGGACTGCTCGAGCTGC 103107	20
TCCGCTGTAGCACACAGGCC 102690	20
TCCGCTGTAGCACTTCAGGCT 51141	21
TCCGCTGTAGCACTTCAGGCTA 28841	22
TCCTCTGTAGCACACAGGC 27854	19
TCGATAAACCTCTGCATCCAG 27829	21
TCGGACTGCTCGAGCTGC 26830	18
```

## 3. Mapping of selected sRNA size ranges.

Once the working population of sRNA reads have been filtered, further processing can be done selecting customised sRNA sizes. In addition to the user defined whole range of sRNA reads (here 18-30nt), other usual choices are the 21-22 nt and 24 nt categories.
This time reads are mapped using ShortStack which is a wrapper of bowtie designed to also handle multimapping reads after the several optional algorithms.

```
nohup bash sRNA.Mapper.range_1830nt.sh
```

**genomic_RPM_values.genes.txt**, **genomic_RPM_values.tes.txt**
```
Chr1    23146   31227   AT1G01040       .       +       1.245824286
Chr1    28500   28706   AT1G01046       .       +       0
Chr1    31170   33153   AT1G01050       .       -       0.5813833333
Chr1    33666   37840   AT1G01060       .       -       0.7474928571
```

**gene.counts.txt**, **te.counts.txt**


**gene.reads**, **te.reads**
```
Chr1	16990	17013	776187	Chr1	16882	17009	AT1TE00020
Chr1	17004	17027	84575	Chr1	16882	17009	AT1TE00020
Chr1	17004	17027	84575	Chr1	17023	18924	AT1TE00025
Chr1	17856	17874	178759	Chr1	17023	18924	AT1TE00025
```
**srna_ids.txt**
```
520486  GATTCAAGACTTCTCGGTACT
84989   GTACTGCAAAGTTCTTCCGCC
155073  ACTGCAAAGTTCTTCCGCCTGAT
231101  ACTGCAAAGTTCTTCCGCCT
```


## 4. Checking replicability and data structure in the data set. Multivariate analysis  and correlogram.

![This is an image](/images/figure1.png)

## 5. sRNA size distribution over over genes and transposable elements.


![This is an image](/images/figure2.png)

*Figure 2. sRNA size distribution over genes and transposable elements*


 
# REFERENCES

Modified versions of this pipeline have been used to process the sRNA datasets in the following papers:

1. Martinez G et al (2018). Paternal easiRNAs regulate parental genome dosage in *Arabidopsis*. **Nature Genetics** 50 (2) 193-198.

2. Wang Z et al (2020). Polymerase IV Plays a Crucial Role in Pollen Development in *Capsella*. **Plant Cell** 32 (4) 950-966.

# CONTACT
juan.sverige at slu.se
