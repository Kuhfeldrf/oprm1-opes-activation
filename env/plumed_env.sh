# Source before any GROMACS+PLUMED run on ORCA.
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
source $ROOT/spack-install/linux-zen4/gromacs-2025.3-45hxd2evh2mxcjqopqa56cr3atcbwlrv/bin/GMXRC
export PLUMED_PREFIX=$ROOT/oprm1-opes-activation/opt/plumed-2.10.1
export PATH=$PLUMED_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$PLUMED_PREFIX/lib:${LD_LIBRARY_PATH:-}
export PLUMED_KERNEL=$PLUMED_PREFIX/lib/libplumedKernel.so
export PY=$ROOT/miniforge3/envs/oprm1-opes/bin/python
