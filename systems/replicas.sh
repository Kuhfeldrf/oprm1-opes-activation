#!/usr/bin/env bash
#SBATCH --job-name=oprm1_reps
#SBATCH --partition=long
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:4
#SBATCH --mem=64G
#SBATCH --output=logs/reps_%x_%j.log
# v0.3 (PLAN_v03.md): N independent UNBIASED replicas of one equilibrated system, one GPU
# each (gmx_mpi -multidir; no coupling between replicas). Each replica starts from the
# end of the restrained equilibration (npt_10) with fresh velocities (its own seed) and
# logs the measurement-only PLUMED observables (COLVAR_measure.<k>).
#
# Usage: sbatch -J reps_<system> systems/replicas.sh <system> <run_tag> <ns_per_replica>
#               [TEMPLATE=plumed_measure.dat] [SEED_BASE=20261300] [MAXH=70]
# MAXH must sit below the job --time. The 3-day default holds a free_* system (500 ns
# at a measured 282 ns/day is 43 h) but not an Nb39 one (165 ns/day, 73 h): submit
# those with --time=5-00:00:00 MAXH=118, or they stop at 70 h and need a restart.
# Restart: resubmit the same command; replicas continue from their checkpoints.
set -eo pipefail
SYSN=$1; TAG=$2; NS=$3; shift 3
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source $ROOT/env/plumed_env_patched.sh
NR=${SLURM_NTASKS:-4}
SYS=$ROOT/systems/$SYSN
RUN=$ROOT/runs/$TAG/$SYSN
mkdir -p "$RUN"; cd "$RUN"
declare -A KV=( [TEMPLATE]=plumed_measure.dat [SEED_BASE]=20261300 [MAXH]=70 )
for kv in "$@"; do KV[${kv%%=*}]=${kv#*=}; done
[ -f $SYS/npt_10.gro ] || { echo "FATAL: $SYSN not equilibrated (no npt_10.gro)"; exit 1; }

DIRS=()
for ((k=0; k<NR; k++)); do
  R=r$k; DIRS+=($R); mkdir -p $R
  if [ ! -f $R/md.tpr ]; then
    sed -e "s/^gen-vel .*/gen-vel                 = yes\ngen-temp                = 310\ngen-seed                = $((${KV[SEED_BASE]} + k))/" \
        -e "s/^continuation .*/continuation            = no/" \
        -e "s/@NSTEPS@/$(( NS * 250000 ))/" \
        <(cat $ROOT/systems/mdp/_common.mdp $ROOT/systems/mdp/prod.mdp) > $R/md.mdp
    $GMX_GROMPP grompp -f $R/md.mdp -c $SYS/npt_10.gro -p $SYS/system.top -n $SYS/index.ndx \
        -o $R/md.tpr -po $R/mdout.mdp -quiet -maxwarn 1   # Lipid21 C=C at 4 fs, see equilibrate.sh
    echo "replica $k of $SYSN from npt_10.gro, gen-seed $((${KV[SEED_BASE]} + k))" > $R/origin.txt
  fi
  $PY $ROOT/systems/render_plumed.py $SYS $ROOT/cv/plumed/${KV[TEMPLATE]} $R/plumed.dat
done
CPI=""; [ -f r0/md.cpt ] && CPI="-cpi md.cpt"
mpirun -np $NR --bind-to none gmx_mpi mdrun -multidir ${DIRS[@]} -deffnm md $CPI \
    -plumed plumed.dat -ntomp ${SLURM_CPUS_PER_TASK:-16} -nb gpu -pme gpu -bonded gpu \
    -maxh ${KV[MAXH]}
grep -h Performance r*/md.log | tail -$NR
