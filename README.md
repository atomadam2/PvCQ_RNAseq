# PvCQ_RNAseq
Code and additional data for P. vivax RNA-seq studies from patient samples performed in the Serre Lab.

This repository contains all additional code and files used in the following manuscript:

In vivo P. vivax gene expression analyses reveal stage-specific chloroquine response and 2 differential regulation of male and female gametocytes.

List of scripts with descriptions:

count_reads.pl - perl script used to count all reads that overlap annotated genes, taking into account read length and gene structure. Input: Bam file of reads aligned to genome, GFF file of genome annotation. Output: List of all annotated exons, with number of reads aligned.

count_assemble.pl - perl script used to compile output files from count_reads.pl from many samples. This will create a matrix of counts for all genes for all samples that can be used directly for downstream applications (EdgeR, PCA, etc.) Input: All output files from count_reads.pl. Output: Count matrix of genes (rows) and samples (columns).

expression_simulation.R - R script used to simulate the gene expression differences in response to CQ between properly paired samples (same patient) vs randomly paired samples (different patients). Input: count matrix from count_assemble.pl with TPM conversion (all_patient_TPM.txt). Output: results of 100 random comparisons and final entry is the correct comparison.

CQ_EdgeR.R - R script used to run EdgeR to find differentially expressed genes between samples before and after chloroquine. Input: count matrix with before chloroquine samples subsetted to have similar total read counts as after chloroquine (all_patient_count_subset.txt). Output: File with differentially expressed genes (Supplementary Table 4) and volcano plot (Figure 3C).

Pv48_Pv47_plot.R - R script used to make scatter plot comparing expression of Pv48 vs Pv47 between all samples. Input: TPM converted expression data from all untreated samples for just Pv48 and Pv47 (Pv48_Pv47.txt). Output: Scatterplot of Pv48 vs Pv47 expression (Figure 2A).

PCA_untreated.R - R script for PCA of all untreated samples. Input: count matrix from all patients (all_patient_count.txt). Output: File containing values for first 20 PCs and scatterplot of PC1 and PC2 from PCA for all untreated samples (Supplementary Figure 2).

CQ_microscopy_proportion.R - R script to make the boxplot and scatterplot of the effect of chloroquine on P. vivax as measured by microscopy and proportion of reads aligned to P. vivax. Input: values for % decrease as observed by microscopy and proportion of reads aligned (parasitemia_reads.txt). Output: Boxplot of % decrease (Figure 3A) and scatterplot of microscopy vs proportion of reads (Supplementary Figure 4).

PCA_chloroquine.R - R script for PCA of sample before and after chloroquine treatment. Input: count matrix with before chloroquine samples subsetted to have similar total read counts as after chloroquine (all_patient_count_subset.txt). Output: File containing values for first 20 PCs, scatterplot of PC1 and PC2 from PCA for all samples before and after chloroquine (Figure 3B), and scatterplot of PC1 and PC2 from PCA with arrows connecting samples before and after chloroquine (Supplementary Figure 6).

Gametocyte_Gene_Cluster.R - R script for 21 known gametocyte gene expression correlation across all samples. Input: Pearson's correlation coefficient matrix for all 21 genes (Gametocyte_Corr.txt). Output: pvclust determined p values for clustering, and Dendrograms and heatmap (Figure 2B).

List of files with descriptions:

all_patient_count.txt - All patient gene expression data in counts (reads aligned to genes) with GeneID, Name of Gene, Number of Exons, and Gene Length (bp).

all_patient_count_subset.txt - All patient gene expression data in counts (reads aligned to genes) after before chloroquine samples were subset to have the same total read numbers as after chloroquine, with GeneID, Name of Gene, Number of Exons, and Gene Length (bp).

all_patient_TPM.txt - All patient gene expression data converted to TPM (transcripts per million) with GeneID, Name of Gene, Number of Exons, and Gene Length (bp).

Pf_Pv_Poran_scRNA_TPM.txt - Single cell RNA-seq data from P. falciparum converted to TPM (transcripts per million) and gene names converted to P. vivax orthologs. Contains data for rings, troph30, troph 36, and schizont, 100 cells each. Data from Poran 2017, https://github.com/KafsackLab/scRNAseq-Malaria 

Pb_Pv_Reid_scRNA_TPM.txt - Randomly generated mock RNA-seq data using single-cell data from P. berghei converted to TPM (transcripts per million) and gene names converted to P. vivax orthologs. Each sample is the sum of 20 cells, with varying proportions of trophozoites and schizonts. Data from Reid 2018, https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5871331/

Pf_Pv_Otto_RNA_TPM.txt - RNA-seq data from in vitro synchronized P. falciparum onverted to TPM (transcripts per million) and gene names converted to P. vivax orthologs. Data from Otto 2010, https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-2958.2009.07026.x

parasitemia_reads.txt - Data for effect of chloroquine on P. vivax as measure by microscopy and by proportion of reads aligned.

Pv48_Pv47.txt - TPM converted expression values for 2 gametocyte genes (Pv48 and Pv47),

Gametocyte_Corr.txt - Pearson's correlation coefficients for 21 known gametocyte genes and their expression across samples.



