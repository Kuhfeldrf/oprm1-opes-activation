#!/usr/bin/env bash
#SBATCH --job-name=g3_qc
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --output=logs/g3_qc_%j.log
set -eo pipefail
cd /scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pilot
mkdir -p results/g3_qc
for s in free_dam free_nlx nb39_dam nb39_nlx; do
  echo "===== membrane QC $s"
  python analysis/membrane_qc.py systems/$s/prod.tpr systems/$s/prod.xtc --stride 20 \
      --json results/g3_qc/membrane_qc_$s.json 2>/dev/null | tail -25
done
for s in nb39_dam nb39_nlx; do
  echo "===== Nb39 contacts $s"
  python analysis/nb_contacts.py systems/$s
done
echo DONE
