### *Juan Santos*, September 2022


# INTRODUCTION

This is a general pipeline for analysis of small RNA sequencing data. It has mainly been developed with focus on the *Arabidopsis* TAIR10 genome, but it should be fully functional with other organisms. The scripts do not require command line arguments to be run; input and reference genomic files have to be specified by editing the script in a text editor. Therefore, some very basic knowledge of linux and R is required. The processes are generously commented in hashes in the script with complementary instructions and hints (worth to read as a complement to this instruction).
 
# SUPPORTED PLATFORMS

Linux, Mac OS

The bash scripts (*.sh) of this software were designed and tested using GNU bash (v4.4.20) in a Ubuntu 18.04 linux system. The R scripts were tested in using R console (v4.1.1) in macOS Monterey.

# PREREQUISITES


fastqc (v0.11.9)

cutadapt (v4.1 with Python 3.8.5)

bowtie (v1.3.1)

ShortStack (v3.8.5)

bedtools (v2.26.0)

samtools (v1.3.1)

# SETTING UP THE WORKING ENVIRONMENT AND THE GENOMIC REFERENCE FILES

The raw data (fastq files) are allocated in sample folders following the file structure below. Most intermediary and final output files will be generated in respective sample folders. In the example there are four samples: two conditions (*mutant* and *wt*) with two replicates each (*rep1* and *rep2*).

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

The working environment and samples have to be specified by editing the script i.e.

```
#working directory where all the sample directories are allocated
working_dir="/media/diskb/nicolas_tmp/REPLICATES_TOTAL";

#this defines the samples to be analysed (dir names)
sample_list="atrm1_rep1  atrm1_rep2  wt_rep1  wt_rep2";
```

Relevant reference genomic files have to be specified as well the beginning of the shell scripts if needed. Ordinary fasta files and bowtie indexed fasta files (bowtie-build) should be available. If removing structural RNAs is wanted prior to downstream processing, a fasta file with the selected structural RNAs (including e.g. pre-tRNA, snoRNA, snRNA, rRNA) has to be prepared (bedtools getfasta) and bowtie-indexed accordingly.

Annotation files (gtf or gff3) have to be transformed to a 6-column bed format. Two annotation bed files, one for genes and one for transposable elements, are used here.

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

To prepare bed files out of gtf or gff3 files can be tricky. The gff2bed tool from BEDOPS suit is an option. Another, sometimes more pragmatic, is processing with a combination of linux regular expressions and/or manual editing in a text editor.

# INSTALLATION

Shell scripts can be downloaded and run directly on a linux server.

# USAGE

```
nohup bash sRNA.Filter.REPLICATES_TOTAL.sh
```
# SIZE DISTRIBUTION

![This is an image](/images/plot1.png)

# REFERENCES

Parts of this pipeline have been used in e.g. the following papers:

1. Martinez G et al (2018). Paternal easiRNAs regulate parental genome dosage in *Arabidopsis*. **Nature Genetics** 50 (2) 193-198.

2. Wang Z et al (2020). Polymerase IV Plays a Crucial Role in Pollen Development in *Capsella*. **Plant Cell** 32 (4) 950-966.


# CONTACT
juan.sverige at outlook.com

#EXAPLE OF GITHUB
https://github.com/phoenixding/rmethrafo/blob/master/README.md
