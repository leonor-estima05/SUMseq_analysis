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
<img width="650" height="591" alt="image" src="https://github.com/user-attachments/assets/7767e5b6-6611-46b1-bb97-213128aa6ad3" />
<img width="650" height="624" alt="image" src="https://github.com/user-attachments/assets/81fd86db-8646-4e2e-927c-8f79b8cf15a8" />

Violin plots support the low depth concerns in the iMono data, UMIs/nucleus is clustering at below 500. However, the original SUM-seq paper's hiPSC-derived macrophages also has a low depth of around 300 UMI/nucleus. As expected the nuclei cluster at near 0% mitochondrial reads. Due to low depth, it is likely a lot of this analysis will be exploratory, due to increased risk of statistical noise and artificats impacting the analysis.

<img width="992" height="352" alt="image" src="https://github.com/user-attachments/assets/6cbcc28f-3533-45e9-8678-31e03f84b187" />

This scrublet simulation histogram shows that scrublet detection will not work on this iMono sample. The left graph shows that the sample doublet scoring generally falls below scrublet's 0.28's threshold, and on the right scrublet's simulated doublets also fall below the 0.28 threshold. Therefore, even if doublets are present in the iMono sample, they won't get detected and therefore won't get removed. Scrublet will not be run on this sample.

This sample was filtered using min_cells = 3 for gene filtering, min_genes = 150 for cell filtering, N_PCS = 15 for principal component inclusion (although the figure indicated PC1 accounted for nearly all data variation) and mitochondrial cut off of over 5% as these are nuclei-only reads.

# iMono Clustering and Gene Expression Results (Scanpy)
Leiden Clustering was the chosen clustering method, and at first 3 resolutions were tested (0.2, 0.5, 1.0) to then choose the appropriate one.
<img width="1292" height="377" alt="image" src="https://github.com/user-attachments/assets/6780f0b0-89e9-4b71-9ed3-01f7b02c0028" />

Resolution 1.0 was overdone at 11 clusters, it is likely that with a low-information sample such as this most clusters are arbitrary. Resolution 0.5 was chosen for further clustering because further analysis revealed each cluster was expressing sufficiently different marker genes. To determine which clusters to merge, violin plots were created based on gene categories (monocyte-associated, macrophage-associated, pluripotent-associated, and proliferation-associated genes) to determine which clusters most expressed each gene's category. 

<img width="1357" height="342" alt="image" src="https://github.com/user-attachments/assets/15a228f0-78b6-44e6-a5d6-95cf6c782545" />

Clusters 0, 1, 2, 5 were grouped into one "myeloid-like" cluster as they were expressing genes highly associated with both monocytes and macrophages (so myeloid cells in general) such as CD97 and SLC11A1. Cluster 3 was kept and expressed the most pluripotency-associated genes such as XACT and ROR1. Cluster 4 was kept as a low-information cluster as it was not expressing various category-associated genes. The final clustering is below.
<img width="1050" height="680" alt="image" src="https://github.com/user-attachments/assets/601aef2c-5f20-4863-a475-c5115b58526d" />

Based on this clustering and a previous tracer analysis on the iMono sample which showed only a subsect of the sample was expressing pluripotency and proliferation-associated genes it is likely there is a less- or undifferentiated group of myeloid cells in the sample.


**Note: As this analysis is ongoing, the scripts and organization may change at any time**
