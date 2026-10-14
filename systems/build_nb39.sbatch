#!/usr/bin/env bash
#SBATCH --job-name=build_nb39
#SBATCH --partition=normal
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --output=logs/build_nb39_%j.log
set -eo pipefail
cd /scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pilot
python -u systems/build_topology.py --root systems --starts nb39_dam nb39_nlx --nb --lig DAM:damgo:N1 NLX:naloxone:N04
