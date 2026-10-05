# SUMseq_analysis
Analysis of paired scRNA and scATAC data from hiPSC-derived monocytes (iMono) and primary monocytes (classical, intermediate, and non-classical) generated with an optimized SUM-sequencing protocol.

The scripts and analyses are part of my Undergraduate Capstone Project at Erasmus University College, I assisted with the data gathering during my research internship at Erasmus Medical Center.

# Background
In modern research hiPSCs are commonly used as models for monocyte subsets and other immune cells. There is also increasing attention on stem cell transplantation therapies for chronic disorders, however there is a lack of research comparing simultaneous scRNA and scATAC data between primary monocyte and these hiPSC models. This project aims to address this gap, including addressing any proliferation or pluripotency-associated open loci. This project uses SUM-sequencing to gather the simultaneous transcriptome and chromatin profiles (single-cell multi-omic sequencing, Yildiz et al. 2026, DOI: 10.1038/s41596-025-01310-0). The project aims to answer the following questions:

1. How does the iMono transcriptome compare to primary monocyte subsets?
2. Are iMono cells transcriptionally homogeneous, or do subpopulations exist?
3. Do iMonos retain accessibility at pluripotency-associated loci?

# Structure

snakemake_configuration: Contains custom snakemake configuration files, edited to fit the data and version of snakemake used in this project 
RNA_analysis: scRNA QC metrics/graphs, clustering, differential expression. Also contains shell scripts for running with slurm.
ATAC_analysis: scATAC QC metrics/graphs

# Pipeline
Data preprocessing was done using the 9.26.1 version of the SnakeMake pipeline for SUM-sequencing from the Zaugg Group (https://git.embl.org/grp-zaugg/SUMseq) which includes BCL conversion to fastq files, demultiplexing, and STARsolo alignment for per-sample count matrices. It also includes QC of the sequencing read depth. Running it required adaptions to the ErasmusMC system, as seen in the configuration files.

# Order of Analysis

RNA
1. [`01_qc_RNA.py`](RNA_analysis/01_qc_RNA.py) - QC, Filtering, Doublet Detecting, Normalization, Feature Selection
2. [`02_clustering_RNA.py`](RNA_analysis/02_clustering_RNA.py) - PCA, neighbor UMAP, Leiden Clustering, Marker Gene Expression
3. [`DESeq2_SUMseq.R`](RNA_analysis/DESeq2_SUMseq.R) - Differential expression analysis
4. [`ClusterProfiler_SUMseq.R`](RNA_analysis/ClusterProfiler_SUMseq.R) - Enrichment

ATAC
1. [`ArchR_SUMseq.R`](ATAC_analysis/ArchR_SUMseq.R) - QC, Peak Calling, Cell Clustering on Accessibility, Marker Scores, Gene Activity Scores, Motif Enrichment

scATAC + scRNA (Gene Regulatory Network Analysis) **Script has not been added yet**
1. SCENIC+ - Inferring Gene Regulatory Networks

# Data Limitations
<img width="672" height="147" alt="image" src="https://github.com/user-attachments/assets/a4391994-157b-461c-b190-cf11bf80c6b6" />

Due to extremely low median gene count/cell, lack of expression of marker genes such as CD14, CD16, CD300E, CD45, etc. the primary monocyte's RNA data has been excluded from this analysis.

# iMono Quality Control Results (Scanpy)
<img width="995" height="912" alt="image" src="https://github.com/user-attachments/assets/7767e5b6-6611-46b1-bb97-213128aa6ad3" />
<img width="997" height="957" alt="image" src="https://github.com/user-attachments/assets/81fd86db-8646-4e2e-927c-8f79b8cf15a8" />

Violin plots support the low depth concerns in the iMono data, UMIs/nucleus is clustering at below 500. However, the original SUM-seq paper's hiPSC-derived macrophages also has a low depth of around 300 UMI/nucleus. As expected the nuclei cluster at near 0% mitochondrial reads. Due to low depth, it is likely a lot of this analysis will be exploratory, due to increased risk of statistical noise and artificats impacting the analysis.

<img width="992" height="352" alt="image" src="https://github.com/user-attachments/assets/6cbcc28f-3533-45e9-8678-31e03f84b187" />

This scrublet simulation histogram shows that scrublet detection will not work on this iMono sample. The left graph shows that the sample doublet scoring generally falls below scrublet's 0.28's threshold, and on the right scrublet's simulated doublets also fall below the 0.28 threshold. Therefore, even if doublets are present in the iMono sample, they won't get detected and therefore won't get removed. Scrublet will not be run on this sample.


**Note: As this analysis is ongoing, the scripts and organization may change at any time**
