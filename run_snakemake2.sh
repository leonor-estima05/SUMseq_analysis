#!/usr/bin/env bash
#SBATCH --mem=100GB
#SBATCH --cpus-per-task=20

source ~/.bashrc
source /home/l.estima/miniforge3/etc...

conda activate snakemake

cd /home/projects/icell/sumseq_pilot/RNA/SUMseq/src/workflow

DRYRUN="--dry-run"

snakemake \
  --profile profiles/slurm \
  --configfile example/input/config_leo.yaml \
  --rerun-incomplete \
  $DRYRUN

conda deactivate
