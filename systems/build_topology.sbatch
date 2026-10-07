#!/usr/bin/env bash
#SBATCH --job-name=oprm1_topo
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/topology_%j.log
# ParmEd conversion of an 85k-atom system takes ~20 min; the ORCA login node kills it.
# Usage: sbatch systems/build_topology.sbatch [--write inactive]
set -eo pipefail
cd /scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pilot
python -u systems/build_topology.py --root systems --starts active inactive "$@"
