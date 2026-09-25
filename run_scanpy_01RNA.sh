#!/usr/bin/env bash
#SBATCH --mem=50GB
#SBATCH --cpus-per-task=10
#SBATCH --output=/home/projects/icell/le_RNA_analysis/scanpy/logs/scanpy_%j.out
#SBATCH --error=/home/projects/icell/le_RNA_analysis/scanpy/logs/scanpy_%j.err

set -eo pipefail

source ~/.bashrc
source "/home/l.estima/miniforge3/etc/profile.d/conda.sh"

cd "/home/projects/icell/le_RNA_analysis/scanpy"

conda activate scanpy_env

export MPLBACKEND=Agg

python -u 01_qc_RNA.py

conda deactivate




