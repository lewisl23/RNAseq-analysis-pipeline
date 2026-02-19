# RNA seq analysis of GSE268197 

This is a RNA seq analysis coursework for MSc Bioinformatics that analyses the 
RNA sequencing datasets from the paper "Effect of GRP75 deficiency on gene expression in DN3 thymocytes".
The study sequences wildtype and Hspa9 cKO DN3 thymocytes using Illumina HiSeq 2000 with 3 samples
for each condition. This project performed differential gene expression and functional enrichment analysis 
on the datasets to understand changes in gene expression and biological insights related to the 
change.

## Methods

### Raw reads from Illumina HiSeq 2000
Taken from "Effect of GRP75 deficiency on gene expression in DN3 thymocytes" \
Accessed through GSE268197 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE268197

Wildtype
- GSM8287690
- GSM8287691
- GSM8287692

cKO (knockout)
- GSM8287693
- GSM8287694
- GSM8287695

### 1. Quality control and alignment
- Illumina reads undergoes quality control using **FastQC** and **MultiQC** to combine the QC reports together
- Alignment using **STAR** with the mouse reference genome to produce aligned BAM file
- Count the number of reads that overlaps with the genome annotation file using **featureCounts**

### 2. Exploratory analysis with PCA and Heatmap
- Data normalisation with rlog to ensures different samples are comparable
- 3D and 2D PCA analysis to identify whether the 2 conditions have distinct variations or 
unwanted variations within the 3 samples of each condition
- Heatmap is drawn with clustering to identify the distance of relationships between different samples

### 3. Differential gene expression analysis
- Data normalisation with Size factors to ensure comparatability
- Differential gene expression with **DEseq2** to output the MA plot of the datasets, which
allows visualisation of the overall diffferential gene expression of the control vs knockout
- Removal of insignificant differences using adjusted p-value of 0.5 (multiple testing problem) 
and selecting the top 10 genes with the highest absolute log2 fold change (heatmap)

### 4. Functional enrichment analysis using GeneOntology (GO)
- Using **mart** to access "mmusculus_gene_ensembl" and Gene Ontology database
- Search for biological process, molecular function, and cellular component