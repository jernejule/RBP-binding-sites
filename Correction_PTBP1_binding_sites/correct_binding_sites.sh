#!/bin/bash -l
#$ -l h_vmem=1G
#$ -l tmem=1G
#$ -l h_rt=1:0:0
#$ -j y
#$ -S /bin/bash

# This script will correct current binding sites (clusters) by find the closest 
# read-end position, downstream of each cluster, considering only reads for 
# which the # nucleotide preceding read start overlaps with the cluster. 
# In the first part it will consider complete iCLIP reads (that had the 3’ 
# adapter trimmed) and in the second part it will consider incomplete  reads 
# (including those that did not have the 3’  adapter trimmed).
#
# Input: existing clusters, comple and incomplete reads
#
# BEDtools software package needs to be installed first: 
# http://bedtools.readthedocs.org/en/latest/

export PATH=/home/programs/bedtools-2.17.0/bin:$PATH

clusters=$1
complete_reads=$2
incomplete_reads=$3
path=`pwd -P`

# 1. remove all reads that don't overlap with cluster ends
bedtools intersect -s -a  ${path}${complete_reads} -b  ${path}${clusters}.BED -wa >  ${path}${complete_reads}-candidates.bed
rm ${path}${complete_reads}-End.BED

# 2. Find closest ends
python getEnd-BED.py ${path}${complete_reads}-candidates.bed ${path}${complete_reads}-candidates-End.bed
bedtools closest -s -D a -io -iu -t first -a ${path}${clusters}.BED -b ${path}${complete_reads}-candidates-End.bed > ${path}${clusters}-closest-iCLIP-end.BED
rm ${path}${complete_reads}-candidates.bed
rm ${path}${complete_reads}-candidates-End.bed

# 3. correct ends if they are and if they are less then 40 nts away which is the maximum read length 
python correctClusterEnds.py ${path}${clusters}-closest-iCLIP-end.BED ${path}${clusters}-closest-iCLIP-end-corrected.bed
rm ${path}${clusters}-closest-iCLIP-end.BED

# 4. Find candidates from untrimmed reads (no adapter) and then find the closest End to uncorrected clusters by using PTB-iiCLIP without adapter removal
bedtools intersect -s -a ${path}${clusters} -b ${path}${incomplete_reads}.BED -wa > ${path}${incomplete_reads}-candidates.bed
python getEnd-BED.py ${path}${incomplete_reads}-candidates.bed ${path}${incomplete_reads}-candidates-End.bed
bedtools closest -s -D a -io -iu -t first -a ${path}${clusters}-closest-iCLIP-end-corrected.bed-NOT.bed -b ${path}${incomplete_reads}-candidates-End.bed > ${path}${clusters}-closest-iCLIP-end-corrected.bed-NOT.bed-closest.bed
rm ${path}${incomplete_reads}-candidates.bed
rm ${path}${incomplete_reads}-candidates-End.bed
rm ${path}${clusters}-closest-iCLIP-end-corrected.bed-NOT.bed

# 5. correct ends of no_adapters iCLIP reads if they are and if they are less then 40 nts away which is the maximum read length
python correctClusterEnds.py  ${path}${clusters}-closest-iCLIP-end-corrected.bed-NOT.bed-closest.bed  ${path}${clusters}-closest-iCLIP-end-corrected.bed-NOT.bed-closest-corrected.bed

# 6. merge clusters together: corrected from both (adapter and no_adapter) and NOT corrected after the previous step
cat ${path}${clusters}-closest-iCLIP-end-corrected.bed ${path}${clusters}-closest-iCLIP-end-corrected.bed-NOT.bed-closest-corrected.bed ${path}${clusters}-closest-iCLIP-end-corrected.bed-NOT.bed-closest-corrected.bed-NOT.bed >  ${path}${clusters}-merged.bed
rm ${path}${clusters}-closest-iCLIP-end-corrected.bed
rm ${path}${clusters}-closest-iCLIP-end-corrected.bed-NOT.bed-closest-corrected.bed
rm ${path}${clusters}-closest-iCLIP-end-corrected.bed-NOT.bed-closest-corrected.bed-NOT.bed

# 7. sort
sort -k1,1 -k2,2n -k5,5 ${path}${clusters}-merged.bed > ${path}${clusters}-merged-sorted.bed
rm ${path}${clusters}-merged.bed

# 8. merge clusters that are less then 20 nt apart (By selecting count >1 we can get all merged clusters)
bedtools merge -i ${path}${clusters}-merged-sorted.bed -s -d 20 -c 5,5,6 -o distinct,count,distinct > ${path}${clusters}-corrected.bed
rm ${path}${clusters}-merged-sorted.bed


