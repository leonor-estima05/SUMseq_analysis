# SUMseq_pilot
Analysis of scRNA and scATAC data of iPSC-derived monocytes and primary monocytes.

The snakemake pipeline by Zaugg Group (https://git.embl.org/grp-zaugg/SUMseq) was used to first demultiplex the data. Configuration files (config_leo.yaml, config.v8+.yaml) and the bash script were edited and used to run snakemake. Python and R scripts are all post-demultiplex analyses intended for either scRNA or scATAC data.

Data was collected following the SUM-Sequencing protocol by Yildiz et al., 2026 (DOI: 10.1038/s41596-025-01310-0) after optimization. 
