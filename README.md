# RBP-binding-sites

#Trimming of adapter sequences
Before mapping the cDNAs, we removed random barcodes and trimmed the 3´ Solexa adapter sequence. Adapter sequences were trimmed by FASTX-Toolkit 0.0.13 adapter removal software, using the following parameters: -Q 33 -a AGATCGGAAG -c -n -l 26. For reads that did not contain parts of the adapter sequence, -C parameter was used, and these were analyzed separately.

#Genome mapping of PTBP1 iCLIP
We used UCSC hg19/GRCh37 genome assembly and bowtie2 2.1 alignment software with default settings. More than 85% (8,585,142 out of 9,634,025) of all cDNAs mapped uniquely to a single genomic position. After mapping, cDNAs with the same random barcode that mapped to the same starting position on the genome were considered to result from PCR amplification and were collapsed to a single cDNA.

#Transcriptome mapping of EIF4A3 iCLIP
For mapping, we compiled a set of representative mRNA sequences from BioMart Ensembl Genes 79, where we used the longest mRNA sequence available for each gene. We mapped to these mRNAs with Bowtie2.1 alignment software, allowing 2 mismatches. More than 50% (11,935,475 of 23,040,243) of cDNAs mapped to a unique mRNA. After mapping, cDNAs with the same random barcode that mapped to the same starting position on an mRNA were considered to result from PCR amplification and were collapsed to a single cDNA.

#Classification of cDNA length
Only cDNAs that mapped to a unique genomic position were evaluated. These were separated into cDNAs that either did (complete) or did not contain parts of the 3´ Solexa primer (incomplete). For libraries produced with the 50 cycle sequencing, the complete cDNAs were further separated into those that were <30nt, 30-34nt or 34-39nt long after trimming. 

#Normalisation of data for drawing of density graphs
All normalisations were performed in R (version 3.1.0) together with “ggplot2” and “smoother” package for final graphical output.

- For analysis of EIFA3 iCLIP, each density graph (RNA map) shows a distribution of cDNA start and end positions relative to positions of exon-exon junctions in mRNAs. All exon-exon junctions within coding regions were taken into account, apart from those that junctioned to the first or last exon in the mRNA. The number of cDNAs starting or ending at each position on the graph was normalised by the number of all cDNAs mapped to representative mRNAs, the mRNA length and the number of examined exon-exon junction positions in the following way: 
RNAmap[n] = ((cDNAs[n] / sum(cDNAs)) * length(mRNAs) / count(exon_exon_junctions)
where [n] stands for a specific position on the density graph.
To draw the graph, we then used Gaussian method with a 3-nucleotide smoothing window. 

- For analysis of PTBP1 iCLIP and CLIP, each density graph (RNA map) shows a distribution of cDNA start and end positions relative to positions of its binding sites. We defined the binding sites in two ways: either by using the position of Y-tracts or extended crosslink clusters in introns. All introns in protein coding genes were taken into account. Due to highly variable abundance of intronic RNAs (and occasional presence of highly abundant non-coding transcripts, such as snoRNAs), we first divided counts at each binding site by dividing them by the count at the MaxCount. To define MaxCount, we examined the region of the binding site, as well as 120 nt 5´ and 3´ of the binding site, to find the nucleotide with the largest count of cDNA count starts or cDNA ends (according to whether starts or ends were plotted on the graph). This avoids from the highly abundant intronic snoRNAs or other abundant introns from dominating the results. The MaxCount-normalised counts of cDNAs starting or ending at each position on the graph was then normalised by the density of all MaxCount-normalised cDNA count/nt starting in the region 50-100 nt downstream of the binding site in the following way: 
RNAmap[n] = (MaxCount-normalised cDNAs[n]) / MaxCount-normalised cDNA density(outside the binding site)
where [n] stands for a specific position on the density graph.
To draw the graph, we then used Gaussian method with a 10-nucleotide smoothing window. 

# Main scripts
 - mapping_to_genome-PTBP1_iCLIP-pipeline.sh (Main script of the mapping pipeline to genome)
 - mapping_to_transcriptome-eIFA3_iCLIP-pipeline.sh (Main script of the mapping pipeline to transcriptome)
 - correct_binding_sites.sh (Main script for correction of binding sites by read-ends)
