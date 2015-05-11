# RBP-binding-sites

#Trimming of adapter sequences
Before mapping the cDNAs, we removed random barcodes and trimmed the 3´ Solexa adapter sequence. Adapter sequences were trimmed by FASTX-Toolkit 0.0.13 adapter removal software, using the following parameters: -Q 33 -a AGATCGGAAG -c -n -l 26. For reads that did not contain parts of the adapter sequence, -C parameter was used, and these were analyzed separately.

#Genome mapping of PTBP1 iCLIP
We used UCSC hg19/GRCh37 genome assembly and bowtie2 2.1 alignment software with default settings. More than 85% (8,585,142 out of 9,634,025) of all cDNAs mapped uniquely to a single genomic position. After mapping, cDNAs with the same random barcode that mapped to the same starting position on the genome were considered to result from PCR amplification and were collapsed to a single cDNA.

#Transcriptome mapping of EIF4A3 iCLIP
For mapping, we compiled a set of representative mRNA sequences from BioMart Ensembl Genes 79, where we used the longest mRNA sequence available for each gene. We mapped to these mRNAs with Bowtie2.1 alignment software, allowing 2 mismatches. More than 50% (11,935,475 of 23,040,243) of cDNAs mapped to a unique mRNA. After mapping, cDNAs with the same random barcode that mapped to the same starting position on an mRNA were considered to result from PCR amplification and were collapsed to a single cDNA.

#Classification of cDNA length
Only cDNAs that mapped to a unique genomic position were evaluated. These were separated into cDNAs that either did (complete) or did not contain parts of the 3´ Solexa primer (incomplete). For libraries produced with the 50 cycle sequencing, the complete cDNAs were further separated into those that were <30nt, 30-34nt or 34-39nt long after trimming. 

#Scripts description:
 - mapping_to_genome-PTBP1_iCLIP-pipeline.sh
 - mapping_to_transcriptome-eIFA3_iCLIP-pipeline.sh
Main scripts of the mapping pipelines.

 - swap_barcode_to_header.py
The script will remove random barcode and experimental barcode from fasta. Random barcode will be saved in the header of fasta file.
