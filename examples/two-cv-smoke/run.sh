#!/usr/bin/env bash
# One command: grompp + 1000-step mdrun with PLUMED + automatic pass/fail check.
# Needs GROMACS >= 2025 (native -plumed) and a PLUMED kernel: export PLUMED_KERNEL=/path/libplumedKernel.so
# CPU only; ~10-60 s.
set -euo pipefail
cd "$(dirname "$0")"
GMX=${GMX:-gmx}
: "${PLUMED_KERNEL:?set PLUMED_KERNEL to libplumedKernel.so (PLUMED >= 2.9)}"
rm -f COLVAR_smoke smoke.tpr smoke.log smoke.edr smoke.gro smoke.cpt smoke.trr mdout.mdp \#*
$GMX grompp -f smoke.mdp -c solute.gro -p solute.top -o smoke.tpr -maxwarn 2 >grompp.out 2>&1
$GMX mdrun -s smoke.tpr -deffnm smoke -plumed plumed_smoke.dat -nb cpu -ntmpi 1 -ntomp "${NT:-4}" >mdrun.out 2>&1
python3 check.py COLVAR_smoke expected_output.json
