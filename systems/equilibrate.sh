#!/usr/bin/env bash
#SBATCH --job-name=oprm1_equil
#SBATCH --partition=long
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=32G
#SBATCH --output=logs/equil_%x_%j.log
# Minimise, heat, staged NPT restraint release, then UNBIASED production with the
# measurement-only PLUMED input (CV time series, no bias) - runbook §8.
# Usage: sbatch -J equil_<start> systems/equilibrate.sh <start> [prod_ns=100]
set -eo pipefail
START=$1; PROD_NS=${2:-100}
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source $ROOT/env/plumed_env.sh
D=$ROOT/systems/$START; M=$ROOT/systems/mdp
cd "$D"
NT=${SLURM_CPUS_PER_TASK:-16}
SEED=$(( 20261007 + $(printf '%d' "'${START:0:1}") ))
mk() { cat $M/_common.mdp "$M/$1.mdp" | sed -e "s/@SEED@/$SEED/" -e "s/@DT@/$2/" \
         -e "s/@NSTEPS@/$3/" -e "s/@FC@/$4/" > "$5"; }
# grompp that tolerates ONLY the known HMR/4 fs warning on Lipid21 C=C bonds (period
# 20 fs = 5 x dt, GROMACS's conservative threshold); any other warning stays fatal.
gpp() {
  if ! gmx grompp "$@" > grompp_try.log 2>&1; then
    n=$(grep -c "^WARNING" grompp_try.log || true)
    k=$(grep -A2 "^WARNING" grompp_try.log | grep -c "estimated oscillational period" || true)
    if [ "$n" -gt 0 ] && [ "$n" = "$k" ]; then
      echo "grompp: accepting $n oscillational-period warning(s) (Lipid21 C=C at 4 fs)"
      gmx grompp "$@" -maxwarn "$n"
    else
      cat grompp_try.log; return 1
    fi
  fi
}
run() { gmx mdrun -deffnm "$1" -ntmpi 1 -ntomp $NT -nb gpu -pme gpu -bonded gpu "${@:2}"; }

if [ ! -f min.gro ]; then
  mk min x x x min.mdp
  gmx grompp -f min.mdp -c system.gro -r system.gro -p system.top -n index.ndx -o min.tpr -maxwarn 1
  gmx mdrun -deffnm min -ntmpi 1 -ntomp $NT
fi
if [ ! -f nvt.gro ]; then
  mk nvt x x x nvt.mdp
  gmx grompp -f nvt.mdp -c min.gro -r min.gro -p system.top -n index.ndx -o nvt.tpr
  run nvt
fi
prev=nvt
# FC kJ/mol/nm^2 : dt : length(ps)
for stage in 1000:0.002:1000 500:0.002:1000 200:0.002:1000 100:0.004:2000 50:0.004:2000 10:0.004:2000; do
  IFS=: read FC DT PS <<< "$stage"
  name=npt_$FC
  if [ ! -f $name.gro ]; then
    NST=$(python3 -c "print(int(round($PS/$DT)))")
    mk npt "$DT" "$NST" "$FC.0" $name.mdp
    gpp -f $name.mdp -c $prev.gro -r $prev.gro -t $prev.cpt -p system.top -n index.ndx -o $name.tpr
    run $name
  fi
  prev=$name
done

# unbiased production with CV logging
NST=$(( PROD_NS * 250000 ))
if [ ! -f prod.tpr ]; then
  mk prod x "$NST" x prod.mdp
  gpp -f prod.mdp -c $prev.gro -t $prev.cpt -p system.top -n index.ndx -o prod.tpr
fi
$PY $ROOT/systems/render_plumed.py "$D" $ROOT/cv/plumed/plumed_measure.dat plumed_measure.dat
if [ -f prod.cpt ]; then CPI="-cpi prod.cpt"; else CPI=""; fi
run prod -plumed plumed_measure.dat $CPI -maxh 46
echo "=== $(date -Is) done"
grep -E "Performance" prod.log | tail -1
