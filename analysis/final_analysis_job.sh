#!/usr/bin/env bash
#SBATCH --job-name=oprm1_final
#SBATCH --partition=normal
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=logs/final_%j.log
# §10 analysis of a finished OPES run: FES (all data and second half), convergence and
# walker error, time series and path overlay, PLUMED REWEIGHT_BIAS cross-check, and the
# bilayer QC on every biased walker trajectory.
# Usage: sbatch analysis/final_analysis_job.sh runs/prod2 results/fes_prod2
set -eo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
cd $ROOT; source env/plumed_env.sh
RUN=${1:?run dir}; OUT=${2:?out dir}; SKIP_FES=${SKIP_FES:-0}
read IB AB < <($PY -c "import json; b=json.load(open('results/endpoints.json'))['opes_proposal']['basin_boxes_for_crossing_count']; print(','.join(map(str,b['inactive'])), ','.join(map(str,b['active'])))")
F="--grid 0.48,1.56,0.0,0.52 --bins 109 --inactive-box $IB --active-box $AB"
if [ "$SKIP_FES" != 1 ]; then
$PY analysis/fes.py --active $RUN/active/w*/COLVAR.* --inactive $RUN/inactive/w*/COLVAR.* --out $OUT $F
$PY analysis/fes.py --active $RUN/active/w*/COLVAR.* --inactive $RUN/inactive/w*/COLVAR.* --out $OUT/second_half --burn 0.5 $F
fi
(cd analysis && $PY convergence.py $ROOT/$RUN --out $ROOT/$OUT --boxes $ROOT/results/endpoints.json --stride 50)

# PLUMED cross-check of the reweighting (runbook §7.6 route)
X=$ROOT/$OUT/plumed_reweight; mkdir -p $X
for s in active inactive; do
  C=$X/COLVAR_all_$s   # absolute: driver runs inside $X
  head -1 $RUN/$s/w0/COLVAR.0 > $C
  for f in $RUN/$s/w*/COLVAR.*; do grep -v '^#' $f >> $C; done
  sed -e "s#@COLVAR@#$C#g" -e "s#@GMIN@#0.48,0.0#" -e "s#@GMAX@#1.56,0.52#" -e "s#@BW@#0.02,0.01#" \
      -e "s#@FES_OUT@#$X/fes_plumed_$s.dat#" cv/plumed/plumed_reweight.dat > $X/reweight_$s.dat
  (cd $X && plumed driver --noatoms --plumed reweight_$s.dat --kt 2.577483 > driver_$s.log 2>&1) \
      || echo "WARNING: plumed reweight cross-check failed for $s (see $X/driver_$s.log)"
  rm -f $C
done

$PY analysis/compare_reweight.py $OUT || true

# bilayer / fold QC on each biased walker
for s in active inactive; do for k in 0 1 2 3; do
  $PY analysis/membrane_qc.py $RUN/$s/w$k/opes.tpr $RUN/$s/w$k/opes.xtc --stride 50 \
      > $OUT/membrane_qc_${s}_w$k.txt 2>/dev/null
done; done
echo DONE
