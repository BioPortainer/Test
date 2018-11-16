#!/bin/bash
#
# initial example of a pipeline script

echo "Initial example of a pipeline script"
echo " "

echo "Download Samples"
echo " "

# Download Samples
wget https://raw.githubusercontent.com/BioPortainer/Test/master/JobRunner/Example_Condition1.fastq
wget https://raw.githubusercontent.com/BioPortainer/Test/master/JobRunner/Example_Condition2.fastq

echo " "
echo "Download Reference"
echo " "

# Download Reference
wget https://raw.githubusercontent.com/BioPortainer/Test/master/JobRunner/Mycoplasma_genitalium.fa

echo " "
echo "Install tools from BioConda"
echo " "

# Install tools from BioConda
conda install bwa bowtie2 -y

echo " "
echo "Create a BWA index in the genomic reference"
echo " "

# Create a BWA index in the genomic reference
bwa index Mycoplasma_genitalium.fa

echo " "
echo "Align the reads in the input file against the genomic reference"
echo " "

# Align the reads in the input file against the genomic reference
bwa aln -I -t 4 Mycoplasma_genitalium.fa Example_Condition1.fastq > out1_bwa.sai
bwa aln -I -t 4 Mycoplasma_genitalium.fa Example_Condition2.fastq > out2_bwa.sai

echo " "
echo "Convert the alignment into a .sam file"
echo " "

# Convert the alignment into a .sam file
bwa samse Mycoplasma_genitalium.fa out1.sai Example_Condition1.fastq > out1_bwa.sam
bwa samse Mycoplasma_genitalium.fa out2.sai Example_Condition2.fastq > out2_bwa.sam

echo " "
echo "Create a bowtie2 index in the genomic reference"
echo " "

# Create a bowtie2 index in the genomic reference
bowtie2-build Mycoplasma_genitalium.fa Mycoplasma_genitalium

echo " "
echo "Create a bowtie2 index in the genomic reference"
echo " "

# Align the reads in the input file against the genomic reference
bowtie2 -x Mycoplasma_genitalium -U Example_Condition1.fastq -S out1_bwt.sam
bowtie2 -x Mycoplasma_genitalium -U Example_Condition2.fastq -S out2_bwt.sam

echo " "
echo "View create files"
echo " "

ls -lh

echo " "
echo "FINISH!"
echo " "
