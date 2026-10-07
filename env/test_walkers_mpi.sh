#!/usr/bin/env bash
#SBATCH --job-name=opes_mwtest
#SBATCH --partition=short
#SBATCH --time=00:40:00
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:2
#SBATCH --mem=32G
#SBATCH --output=mwtest_%j.log
# machinery test: 2 walkers, WALKERS_MPI, 5000 steps (20 ps at 4 fs), from active nvt.gro
set -eo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source $ROOT/env/plumed_env_patched.sh
SYS=$ROOT/systems/active
cd $ROOT/scratch_mwtest
for k in 0 1; do
  mkdir -p w$k
  sed -e "s/^gen-vel .*/gen-vel                 = yes\ngen-temp                = 310\ngen-seed                = $((7 + k))/" \
      -e "s/^continuation .*/continuation            = no/" -e "s/@NSTEPS@/5000/" \
      -e "s/^nstxout-compressed .*/nstxout-compressed = 0/" \
      <(cat $ROOT/systems/mdp/_common.mdp $ROOT/systems/mdp/prod.mdp) > w$k/opes.mdp
  $GMX_GROMPP grompp -f w$k/opes.mdp -c $SYS/nvt.gro -p $SYS/system.top -n $SYS/index.ndx -o w$k/opes.tpr -po w$k/mdout.mdp -quiet -maxwarn 1
  $PY $ROOT/systems/render_plumed.py $SYS $ROOT/cv/plumed/plumed_opes.dat w$k/plumed.dat \
      INDEX=$SYS/index.ndx RESTART_LINE="" STATE_RFILE="" WALKERS=WALKERS_MPI \
      SIGMA=0.03,0.02 UW=1.60,0.60 LW=0.55,0.0 BARRIER=50
  sed -i "s/PACE=500/PACE=50/; s/STRIDE=500/STRIDE=50/" w$k/plumed.dat   # many kernels in 20 ps
done
mpirun -np 2 --bind-to none gmx_mpi mdrun -multidir w0 w1 -deffnm opes -plumed plumed.dat \
    -ntomp 8 -nb gpu -pme gpu -bonded gpu
grep -h -E "Performance" w*/opes.log
grep -h -i -E "GPU|PLUMED" w0/opes.log | head -12
