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
fastx_clipper -Q 33 -a AGATCGGAAG -C -n -l 26 -i  ${path}${data} -o  ${path}${data}-clipped.fq

# fastq to fasta
fastq_to_fasta -Q 33 -n -i ${path}${data}-clipped.fq -o ${path}${data}-clipped.fa
rm ${path}${data}-clipped.fq

# swap random barcodes to headers of fasta file
python ${path}swap_barcode_to_header.py ${path}${data}-clipped.fa ${path}${data}-barcodes.fa
rm ${path}${data}-clipped.fa

# map to hg19
bowtie2-align -x ~/bowtie-indexes/hg19/hg19 -f ${path}${data}-noBarcodes.fa -S ${path}${data}.sam
rm ${path}${data}-noBarcodes.fa

# SAM to BAM
samtools view -hSb ${path}${data}.sam > ${path}${data}.bam
rm ${path}${data}.sam

# convert bam to bed
bedtools bamtobed -i ${path}${data}.bam > ${path}${data}.bed

# remove duplicates
cat ${path}${data}.bed | sort -k1,1 -k2,2n -k4,4 | uniq > ${path}${data}-uniq.bed

# compress the original data
gzip ${path}${data}

