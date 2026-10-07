# Source before any GROMACS+PLUMED run on ORCA.
_PILOT=/scratch/kuhfeldr-Kuhfeld_temp
source $_PILOT/spack-install/linux-zen4/gromacs-2025.3-45hxd2evh2mxcjqopqa56cr3atcbwlrv/bin/GMXRC
export PLUMED_PREFIX=$_PILOT/oprm1-opes-activation/opt/plumed-2.10.1
export PATH=$PLUMED_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$PLUMED_PREFIX/lib:${LD_LIBRARY_PATH:-}
export PLUMED_KERNEL=$PLUMED_PREFIX/lib/libplumedKernel.so
export PY=$_PILOT/miniforge3/envs/oprm1-opes/bin/python
