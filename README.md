### Juan Santos, November 2022

# INTRODUCTION

# SUPPORTED PLATFORMS

Linux, Mac OS

The shell scripts (*.sh) of this software were designed and tested using GNU bash (v4.4.20) in a Ubuntu 18.04 linux system. The R scripts were tested in using R console (v4.1.1) in macOS Monterey.

# PREREQUISITES

fastqc (v0.11.9)
cutadapt (v4.1 with Python 3.8.5)
bowtie (v1.3.1)
ShortStack (v3.8.5)
bedtools (v2.26.0)
samtools (v1.3.1)


# SETTING UPP THE PROJECT

The raw data (fastq files) are allocated in sample folders. Most intermediary and final output files will be generated in respective sample folders. In the example there are four samples: two conditions (mutant and wt) with two replicates each (rep1 and rep2).


REPLICATES_TOTAL/
├── mutant_rep1
│   └── mutant_Rep1_FKDL210177288-1a-3_H3TM7DSX2_L4_1.fq.gz
├── mutant_rep2
│   └── mutant_Rep2_FKDL210177288-1a-4_H3TM7DSX2_L4_1.fq.gz
├── wt_rep1
│   └── Col_Rep1_FKDL210177288-1a-1_H3TM7DSX2_L4_1.fq.gz
└── wt_rep2
    └── Col_Rep2_FKDL210177288-1a-2_H3TM7DSX2_L4_1.fq.gz


# SETTING UP GENOMIC REFERENCES

# INSTALLATION

# USAGE

step 0: setting up the project directory
I usually set up the directory for an RNA-seq analysis like this, if I'm starting with raw reads:

<include picture of file tree>

#TEST EXAMPLE


#REFERENCES

Parts of this pipeline have been used in e.g. the following papers:

Martinez G et al (2018). Paternal easiRNAs regulate parental genome dosage in <i>Arabidopsis</i>. Nature Genetics 50 (2) 193-198

Wang Z et al. (2020). Polymerase IV Plays a Crucial Role in Pollen Development in Capsella. Plant Cell 32 (4) 950-966.


#CONTACT
juan.sverige at outlook.com

#EXAPLE OF GITHUB
https://github.com/phoenixding/rmethrafo/blob/master/README.md
