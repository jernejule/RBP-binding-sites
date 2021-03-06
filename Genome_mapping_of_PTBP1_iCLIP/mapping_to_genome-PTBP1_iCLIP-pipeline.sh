#!/bin/bash -l
#$ -l h_vmem=16G
#$ -l tmem=16G
#$ -l h_rt=32:0:0
#$ -j y
#$ -S /bin/bash

PYTHONPATH=/home/programs/Python-2.7.5/bin
PERL5LIB=/home/programs/ActivePerl-5.16.3.1603-x86_64-linux-glibc-2.3.5-296746/perl/bin
export PATH=$PYTHONPATH:$PATH
export PATH=$PERL5LIB:$PATH
export PATH=/home/programs/tophat-2.0.9.Linux_x86_64:$PATH
export PATH=/home/programs/bowtie2-2.1.0:$PATH
export PATH=/home/programs/samtools-0.1.19:$PATH
export PATH=/home/programs/fastx_toolkit0.0.13:$PATH
export PATH=/home/programs/bedtools-2.17.0/bin:$PATH

data=$1
path=`pwd -P`

gunzip ${path}${data}.gz

# clip the adapter and discard clipped sequences and discard the sequences that are shorter then 17 nt + 5 random barcode + 4 experimental barcode
fastx_clipper -Q 33 -a AGATCGGAAG -C -n -l 26 -i  ${path}${data} -o  ${path}${data}-incomplete.fq
fastx_clipper -Q 33 -a AGATCGGAAG -c -n -l 26 -i  ${path}${data} -o  ${path}${data}-complete.fq

# fastq to fasta
fastq_to_fasta -Q 33 -n -i ${path}${data}-incomplete.fq -o ${path}${data}-incomplete.fa
fastq_to_fasta -Q 33 -n -i ${path}${data}-complete.fq -o ${path}${data}-complete.fa
rm ${path}${data}-incomplete.fq
rm ${path}${data}-complete.fq

# swap random barcodes to headers of fasta file
python ${path}swap_barcode_to_header.py ${path}${data}-incomplete.fa ${path}${data}-incomplete-barcodes.fa
python ${path}swap_barcode_to_header.py ${path}${data}-complete.fa ${path}${data}-complete-barcodes.fa
rm ${path}${data}-incomplete.fa
rm ${path}${data}-complete.fa

# map to hg19
bowtie2-align -x ~/bowtie-indexes/hg19/hg19 -f ${path}${data}-incomplete-barcodes.fa -S ${path}${data}-incomplete.sam
bowtie2-align -x ~/bowtie-indexes/hg19/hg19 -f ${path}${data}-complete-barcodes.fa -S ${path}${data}-complete.sam
rm ${path}${data}-incomplete-barcodes.fa
rm ${path}${data}-complete-barcodes.fa

# SAM to BAM
samtools view -hSb ${path}${data}-incomplete.sam > ${path}${data}-incomplete.bam
samtools view -hSb ${path}${data}-complete.sam > ${path}${data}-complete.bam
rm ${path}${data}-incomplete.sam
rm ${path}${data}-complete.sam

# convert bam to bed
bedtools bamtobed -i ${path}${data}-incomplete.bam > ${path}${data}-incomplete.bed
bedtools bamtobed -i ${path}${data}-complete.bam > ${path}${data}-complete.bed

# remove duplicates
cat ${path}${data}-incomplete.bed | sort -k1,1 -k2,2n -k4,4 | uniq > ${path}${data}-incomplete-uniq.bed
cat ${path}${data}-complete.bed | sort -k1,1 -k2,2n -k4,4 | uniq > ${path}${data}-complete-uniq.bed
rm ${path}${data}-incomplete.bed
rm ${path}${data}-complete.bed

# compress the original data
gzip ${path}${data}

