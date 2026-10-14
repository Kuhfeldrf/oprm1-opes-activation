#!/usr/bin/env bash
#SBATCH --job-name=build_free
#SBATCH --partition=normal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/build_free_%j.log
set -eo pipefail
cd /scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pilot
python -u systems/build_topology.py --root systems --starts free_dam free_nlx --lig DAM:damgo:N1 NLX:naloxone:N04
