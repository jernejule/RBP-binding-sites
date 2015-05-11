# RBP-binding-sites

Trimming of adapter sequences
Before mapping the cDNAs, we removed random barcodes and trimmed the 3´ Solexa adapter sequence. Adapter sequences were trimmed by FASTX-Toolkit 0.0.13 adapter removal software, using the following parameters: -Q 33 -a AGATCGGAAG -c -n -l 26. For reads that did not contain parts of the adapter sequence, -C parameter was used, and these were analyzed separately.

Genome mapping of PTBP1 iCLIP
We used UCSC hg19/GRCh37 genome assembly and bowtie2 2.1 alignment software with default settings. More than 85% (8,585,142 out of 9,634,025) of all cDNAs mapped uniquely to a single genomic position. After mapping, cDNAs with the same random barcode that mapped to the same starting position on the genome were considered to result from PCR amplification and were collapsed to a single cDNA.

Transcriptome mapping of EIF4A3 iCLIP
For mapping, we compiled a set of representative mRNA sequences from BioMart Ensembl Genes 79, where we used the longest mRNA sequence available for each gene. We mapped to these mRNAs with Bowtie2.1 alignment software, allowing 2 mismatches. More than 50% (11,935,475 of 23,040,243) of cDNAs mapped to a unique mRNA. After mapping, cDNAs with the same random barcode that mapped to the same starting position on an mRNA were considered to result from PCR amplification and were collapsed to a single cDNA.

Classification of cDNA length
Only cDNAs that mapped to a unique genomic position were evaluated. These were separated into cDNAs that either did (complete) or did not contain parts of the 3´ Solexa primer (incomplete). For libraries produced with the 50 cycle sequencing, the complete cDNAs were further separated into those that were <30nt, 30-34nt or 34-39nt long after trimming. 

Scripts description:
 - three_prime_fragment-script.sh
Main script of the pipeline.

 - SAMtoCollapsedSAMandBED.py 
The script will read .SAM file and write it to collapsed .SAM file ignoring "4" flag for strand and remove all 
reads with duplicated barcode for each position. Collapsed data will be also written to  a .BED file format.

 - remove_up_stream.py
The scrpt will remove everything up to the end of adapter sequence.

 - setBEDpositions.py
The script will set the chromosome and extend BED positions.

 - swap_barcodes.py
The script will read fasta file and remove random barcode and experimental barcode from fasta. Random barcodes 
will be saved to a new fasta file.

 - number_of_reads_per_nt.py
Script will add a number of reads from BED to each nucleotide position in the transcript. Results will be written into 2 files 
seperated by strand of the binding (same.bed and anti.bed).

 - three-prime-fragment-plots-binning-filtering.R
The script will import 3 prime fragments tables and added aditional columns with binned and normalised values. 
Each one of them will be ploted in 1 nucleotide, 10 nt and 30 nt resolution. Set all paths to "*_same.tab" tables 
from 3 prime fragments results.
