#!/usr/bin/env bash
#SBATCH --mem=100GB
#SBATCH --cpus-per-task=20
#SBATCH --output=/home/projects/icell/sumseq_pilot/ATAC/logs/%x_%j.out 
#SBATCH --error=/home/projects/icell/sumseq_pilot/ATAC/logs/%x_%j.err

set -euo pipefail

source ~/.bashrc
source "/home/l.estima/miniforge3/etc/profile.d/conda.sh"

conda activate snakemake

cd /home/projects/icell/sumseq_pilot/RNA/SUMseq/src/workflow

DRYRUN="--dry-run"

snakemake \
  --profile profiles/slurm \
  --configfile example/input/config.yaml \
  --rerun-incomplete \
  --printshellcmds \
  --show-failed-logs \
  $DRYRUN

conda deactivate
