#!/usr/bin/env bash
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=20

source ~/.bashrc
source "/home/l.estima/miniforge3/etc/profile.d/conda.sh"

conda activate snakemake

cd /home/projects/icell/sumseq_pilot/RNA/SUMseq/src/workflow

DRYRUN=""

snakemake \
  --profile profiles/slurm \
  --configfile example/input/config.yaml \
  --rerun-incomplete \
  $DRYRUN

conda deactivate
