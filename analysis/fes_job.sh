#!/usr/bin/env bash
#SBATCH --job-name=oprm1_fes
#SBATCH --partition=short
#SBATCH --time=03:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --output=logs/fes_%j.log
# §10 analysis of one OPES run: reweighted 2D FES per start, start difference, quarters,
# crossings and figures. Grid spans the wall range (CV1 0.53-1.50 nm, CV2 <= 0.47 nm) plus
# margin; basin boxes are the measured unbiased endpoints (results/endpoints.json).
# Usage: sbatch analysis/fes_job.sh runs/prod2 results/fes_prod2 [burn_fraction]
set -eo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
cd $ROOT; source env/plumed_env.sh
RUN=${1:?run dir}; OUT=${2:?out dir}; BURN=${3:-0}
read IB AB < <($PY -c "import json; b=json.load(open('results/endpoints.json'))['opes_proposal']['basin_boxes_for_crossing_count']; print(','.join(map(str,b['inactive'])), ','.join(map(str,b['active'])))")
$PY analysis/fes.py --active $RUN/active/w*/COLVAR.* --inactive $RUN/inactive/w*/COLVAR.* \
    --out $OUT --burn $BURN --grid 0.48,1.56,0.0,0.52 --bins 109 \
    --inactive-box $IB --active-box $AB
echo DONE
