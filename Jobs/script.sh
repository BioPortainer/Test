#!/bin/bash
#
# initial example of a pipeline script

# Download Samples
wget https://cs.wellesley.edu/~btjaden/Rockhopper/download/Example_Condition1.fastq
wget https://cs.wellesley.edu/~btjaden/Rockhopper/download/Example_Condition2.fastq

# Download Reference
wget https://raw.githubusercontent.com/BioPortainer/Test/master/Jobs/Mycoplasma_genitalium.fa

# Install tools from BioConda
conda install bwa bowtie2 fastqc -y

# Quality check
fastqc Example_Condition1.fastq
fastqc Example_Condition2.fastq

# Create a BWA index in the genomic reference
bwa index Mycoplasma_genitalium.fa

# Align the reads in the input file against the genomic reference
bwa aln -I -t 4 Mycoplasma_genitalium.fa Example_Condition1.fastq > out1_bwa.sai
bwa aln -I -t 4 Mycoplasma_genitalium.fa Example_Condition2.fastq > out2_bwa.sai

# Convert the alignment into a .sam file
bwa samse Mycoplasma_genitalium.fa out1.sai Example_Condition1.fastq > out1_bwa.sam
bwa samse Mycoplasma_genitalium.fa out2.sai Example_Condition2.fastq > out2_bwa.sam

# Create a bowtie2 index in the genomic reference
bowtie2-build Mycoplasma_genitalium.fa Mycoplasma_genitalium

# Align the reads in the input file against the genomic reference
bowtie2 -x Mycoplasma_genitalium -U Example_Condition1.fastq -S out1_bwt.sam
bowtie2 -x Mycoplasma_genitalium -U Example_Condition2.fastq -S out2_bwt.sam
