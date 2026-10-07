#!/usr/bin/env bash
#SBATCH --job-name=oprm1_opes
#SBATCH --partition=long
#SBATCH --time=3-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=16
#SBATCH --gres=gpu:4
#SBATCH --mem=64G
#SBATCH --output=logs/opes_%x_%j.log
# Multi-walker OPES_METAD over (d_tm36, rmsd_npxxy): N walkers share ONE bias
# (WALKERS_MPI), one GPU each, via gmx_mpi -multidir. Runbook §7.5, §9.
#
# Usage: sbatch -J opes_<start> systems/opes.sh <start> <run_tag> <ns_per_walker> \
#               SIGMA=0.0x,0.0y UW=a,b LW=c,d [BARRIER=50]
#   walker k starts from the unbiased production of <start> at time t_k (evenly spaced
#   over the second half), with fresh velocities (seed differs per walker).
# Restart: resubmit the same command; walkers continue from their checkpoints and the
# shared bias from STATE (RESTART=YES).
set -eo pipefail
START=$1; TAG=$2; NS=$3; shift 3
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source $ROOT/env/plumed_env_patched.sh
NW=${SLURM_NTASKS:-4}
SYS=$ROOT/systems/$START
RUN=$ROOT/runs/$TAG/$START
mkdir -p "$RUN"; cd "$RUN"
declare -A KV=( [BARRIER]=50 )
for kv in "$@"; do KV[${kv%%=*}]=${kv#*=}; done
for k in SIGMA UW LW; do [ -n "${KV[$k]}" ] || { echo "FATAL: $k= required"; exit 1; }; done

FIRST_RUN=0
[ -f w0/opes.tpr ] || FIRST_RUN=1
DIRS=()
for ((k=0; k<NW; k++)); do
  W=w$k; DIRS+=($W); mkdir -p $W
  if [ $FIRST_RUN = 1 ]; then
    # starting frame: second half of the unbiased leg, evenly spaced (ps)
    T0=$(python3 -c "print(int(50000 + 50000*($k+1)/$NW))")
    echo 0 | gmx_mpi trjconv -f $SYS/prod.xtc -s $SYS/prod.tpr -dump $T0 -o $W/start.gro -quiet
    sed -e "s/^gen-vel .*/gen-vel                 = yes\ngen-temp                = 310\ngen-seed                = $((20261100 + k))/" \
        -e "s/^continuation .*/continuation            = no/" \
        -e "s/@NSTEPS@/$(( NS * 250000 ))/" \
        <(cat $ROOT/systems/mdp/_common.mdp $ROOT/systems/mdp/prod.mdp) > $W/opes.mdp
    # trajectory every 100 ps is kept; COLVAR carries the CV time series at 2 ps
    $GMX_GROMPP grompp -f $W/opes.mdp -c $W/start.gro -p $SYS/system.top -n $SYS/index.ndx \
        -o $W/opes.tpr -po $W/mdout.mdp -quiet -maxwarn 1   # Lipid21 C=C at 4 fs, see equilibrate.sh
    echo "walker $k from t=$T0 ps of $START unbiased leg" > $W/origin.txt
  fi
  # with WALKERS_MPI only walker 0 writes KERNELS/STATE; on restart every walker reads it.
  # COLVAR is written as COLVAR.<k> (PLUMED appends the replica index in multi-sim mode).
  RLINE=""; RFILE=""
  if [ -f w0/STATE ]; then RLINE="RESTART"; RFILE="STATE_RFILE=$RUN/w0/STATE"; fi
  $PY $ROOT/systems/render_plumed.py $SYS $ROOT/cv/plumed/plumed_opes.dat $W/plumed.dat \
      INDEX=$SYS/index.ndx RESTART_LINE="$RLINE" STATE_RFILE="$RFILE" WALKERS=WALKERS_MPI \
      SIGMA=${KV[SIGMA]} UW=${KV[UW]} LW=${KV[LW]} BARRIER=${KV[BARRIER]}
done
cp $ROOT/cv/ref_inactive_npxxy.pdb ref_used.pdb   # provenance of the CV2 reference
CPI=""; [ -f w0/opes.cpt ] && CPI="-cpi opes.cpt"
mpirun -np $NW --bind-to none gmx_mpi mdrun -multidir ${DIRS[@]} -deffnm opes $CPI \
    -plumed plumed.dat -ntomp ${SLURM_CPUS_PER_TASK:-16} -nb gpu -pme gpu -bonded gpu \
    -maxh 70
grep -h Performance w*/opes.log | tail -$NW
