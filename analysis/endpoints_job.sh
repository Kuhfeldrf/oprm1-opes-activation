#!/usr/bin/env bash
#SBATCH --job-name=oprm1_endpoints
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/endpoints_%j.log
# §6.4 decision point: rebuild the CV2 reference from the equilibrated inactive basin,
# re-measure both unbiased legs against it, and run the bilayer QC over both legs.
set -eo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
cd $ROOT; source env/plumed_env.sh
mkdir -p results
$PY cv/rebuild_ref_equilibrated.py systems/inactive cv/ref_inactive_npxxy.pdb
for s in active inactive; do
  D=$ROOT/systems/$s; W=$ROOT/results/remeasure_$s; mkdir -p $W; cd $W
  $PY $ROOT/systems/render_plumed.py $D $ROOT/cv/plumed/plumed_measure.dat measure.dat STRIDE=1
  rm -f COLVAR_measure*
  plumed driver --plumed measure.dat --ixtc $D/prod.xtc --timestep 100 > driver.log 2>&1
  cp COLVAR_measure $ROOT/results/COLVAR_unbiased_$s.dat
  cd $ROOT
  $PY analysis/membrane_qc.py $D/prod.tpr $D/prod.xtc --stride 20 --json results/membrane_qc_$s.json \
      > results/membrane_qc_$s.txt 2>/dev/null
done
$PY analysis/endpoints.py --active results/COLVAR_unbiased_active.dat \
    --inactive results/COLVAR_unbiased_inactive.dat --skip 10 --out results/endpoints.json
echo DONE
