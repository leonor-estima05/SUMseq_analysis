# SUMseq_pilot
Analysis of scRNA and scATAC data of iPSC-derived monocytes and primary monocytes.

The snakemake pipeline by Zaugg Group was used to first demultiplex the data. Configuration files (config_leo.yaml, config.v8+.yaml) and the bash script were edited and used to run snakemake. Python and R scripts are all post-demultiplex analyses intended for either scRNA or scATAC data.

Data was collected following the SUM-Sequencing protocol by Yildiz et al., 2026 (DOI: 10.1038/s41596-025-01310-0) after optimization. Following data collection the following files were used to create a full computational analysis of the scRNA and scATAC data, with the ultimate goal of comparing the transcriptome and regulatory genome of iPSC-derived monocytes and primary monocytes at the single-cell resolution.
