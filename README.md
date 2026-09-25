# SUMseq_analysis
Analysis of paired scRNA and scATAC data from hiPSC-derived monocytes (iMono) and primary monocytes (classical, intermediate, and non-classical) generated with an optimized SUM-sequencing protocol.

The scripts and analyses are part of my Underaduate Capstone Project at Erasmus University College, I assisted with the data gathering during my research internship at Erasmus Medical Center.

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
1. SCENIC+ - Infering Gene Regulatory Networks

**Note: As this analysis is ongoing, the scripts and organization may change at any time**
