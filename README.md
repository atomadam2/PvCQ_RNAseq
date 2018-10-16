# PvCQ_RNAseq
Code and additional data for P. vivax RNA-seq studies from patient samples performed in the Serre Lab

This repository contains all additional code and files used in the following manuscript:

In vivo P. vivax gene expression analyses reveal stage-specific chloroquine response and 2 differential regulation of male and female gametocytes

List of scripts with descriptions:

count_reads.pl - perl script used to count all reads that overlap annotated genes, taking into account read length and gene structure. Input: Bam file of reads aligned to genome, GFF file of genome annotation. Output: List of all annotated exons, with number of reads aligned

count_assemble.pl - perl script used to compile output files from count_reads.pl from many samples. This will create a matrix of counts for all genes for all samples that can be used directly for downstream applications (EdgeR, PCA, etc.) Input: All output files from count_reads.pl. Output: Count matrix of genes (rows) and samples (columns)

expression_simulation.R - R script used to simulate the gene expression differences in response to CQ between properly paired samples (same patient) vs randomly paired samples (different patients). Input: count matrix from count_assemble.pl with TPM conversion (all_patient_TPM.txt). Output: results of 100 random comparisons and final entry is the correct comparison.


List of files with descriptions:

all_patient_count.txt

all_patient_TPM.txt - All patient gene expression data converted to TPM (transcripts per million) with GeneID, Name of Gene, Number of Exons, and Gene Length (bp).

Pf_Pv_Poran_scRNA_TPM.txt - Single cell RNA-seq data from P. falciparum converted to TPM (transcripts per million) and gene names converted to P. vivax orthologs. Contains data for rings, troph30, troph 36, and schizont, 100 cells each. Data from Poran 2017, https://github.com/KafsackLab/scRNAseq-Malaria 

Pb_Pv_Reid_scRNA_TPM.txt - Randomly generated mock RNA-seq data using single-cell data from P. berghei converted to TPM (transcripts per million) and gene names converted to P. vivax orthologs. Each sample is the sum of 20 cells, with varying proportions of trophozoites and schizonts. Data from Reid 2018, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5871331/

Pf_Pv_Otto_RNA_TPM.txt - RNA-seq data from in vitro synchronized P. falciparum onverted to TPM (transcripts per million) and gene names converted to P. vivax orthologs. Data from Otto 2010, https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2958.2009.07026.x


