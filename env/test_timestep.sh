#!/usr/bin/env bash
#SBATCH --job-name=dttest
#SBATCH --partition=short
#SBATCH --time=01:30:00
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:1
#SBATCH --mem=24G
#SBATCH --output=dttest_%x_%j.log
# usage: sbatch -J <name> dttest.sh <dt> <lincs_order> <lincs_iter> <ps>
set -eo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source $ROOT/env/plumed_env.sh
DT=$1; LO=$2; LI=$3; PS=$4; NAME=$SLURM_JOB_NAME
S=$ROOT/systems/active
mkdir -p $ROOT/scratch_dttest/$NAME && cd $ROOT/scratch_dttest/$NAME
NST=$(python3 -c "print(int(round($PS/$DT)))")
cat $ROOT/systems/mdp/_common.mdp $ROOT/systems/mdp/prod.mdp | sed -e "s/@NSTEPS@/$NST/" \
  -e "s/^dt .*/dt = $DT/" -e "s/^lincs-order .*/lincs-order = $LO/" -e "s/^lincs-iter .*/lincs-iter = $LI/" \
  -e "s/^nstxout-compressed .*/nstxout-compressed = 0/" > t.mdp
gmx grompp -f t.mdp -c $S/npt_10.gro -t $S/npt_10.cpt -p $S/system.top -n $S/index.ndx -o t.tpr -maxwarn 1 -quiet
GMX_MAXCONSTRWARN=-1 gmx mdrun -deffnm t -ntmpi 1 -ntomp 16 -nb gpu -pme gpu -bonded gpu || true
echo "RESULT $NAME dt=$DT order=$LO iter=$LI: lincs_warnings=$(grep -c 'LINCS WARNING' t.log) maxdev=$(grep -o 'max [0-9.e+-]*' t.log | awk '{print $2}' | sort -g | tail -1) $(grep Performance t.log)"
printf "T-Solute\nT-MEMB\nT-SOLV\n\n" | gmx energy -f t.edr 2>/dev/null | grep -E "^T-"
