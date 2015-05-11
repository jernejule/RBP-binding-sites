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

# remove adapter seq and keep all reads that are longer then 26 nts which 17 nts (17nts seq + 5nts random barcode + 4nts experimental barcode)
fastx_clipper -Q 33 -n -l 26 -a AGATCGGAAG -i ${path}${data} -o ${path}${data}-clipped.fq

# convert fastq to fasta
fastq_to_fasta -Q 33 -n -i ${path}${data}-clipped.fq -o ${path}${data}-clipped.fa
rm ${path}${data}-clipped.fq

# swap random barcodes to headers of fasta file
python ${path}swap_barcode_to_header.py ${path}${data}-clipped.fa ${path}${data}-barcodes.fa
rm ${path}${data}-clipped.fa

# map to transcripts
bowtie2-align -x /home/bowtie-indexes/Human-GRCh38.p2-CDS-transcripts/Human-GRCh38.p2-CDS-transcripts -f ${path}${data}-barcodes.fa -S ${path}${data}.sam
rm ${path}${data}-barcodes.fa

# filter out all reads with more then 2 mismatches
samtools view -Sh ${path}${data}.sam | grep -e "^@" -e "XM:i:[012][^0-9]" > ${path}${data}-2mis.sam
rm ${path}${data}.sam

# SAM to BAM
samtools view -hSb ${path}${data}-2mis.sam > ${path}${data}-2mis.bam
rm ${path}${data}-2mis.sam

# convert bam to bed
bedtools bamtobed -i ${path}${data}-2mis.bam > ${path}${data}-2mis.bed

# remove duplicates
cat ${path}${data}-2mis.bed | sort -k1,1 -k2,2n -k4,4 | uniq > ${path}${data}-2mis-uniq.bed

# compress the original data
gzip ${path}${data}



