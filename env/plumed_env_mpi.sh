# Source before multi-walker (WALKERS_MPI) GROMACS+PLUMED runs on ORCA: gmx_mpi + MPI PLUMED.
ROOT=/scratch/kuhfeldr-Kuhfeld_temp
source $ROOT/spack-install/linux-zen4/gromacs-2025.3-47lswoc3a2ef5s5iq33gfyklfut6lt3u/bin/GMXRC
export PATH=/software/spack/v1.1.1/opt/spack/linux-zen4/openmpi-5.0.8-ruqvw4kfg3ftlc6nmn4gm6w7x2rwrtfl/bin:$PATH
export PLUMED_PREFIX=$ROOT/oprm1-opes-activation/opt/plumed-2.10.1-mpi
export PATH=$PLUMED_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$PLUMED_PREFIX/lib:${LD_LIBRARY_PATH:-}
export PLUMED_KERNEL=$PLUMED_PREFIX/lib/libplumedKernel.so
export PY=$ROOT/miniforge3/envs/oprm1-opes/bin/python
