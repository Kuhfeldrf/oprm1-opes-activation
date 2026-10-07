# Source before multi-walker OPES: PLUMED-PATCHED GROMACS 2025.0 (MPI, CUDA) + MPI PLUMED kernel.
# Use this, not plumed_env_mpi.sh, whenever WALKERS_MPI is in play (see NOTES.md, Day 1).
_PILOT=/scratch/kuhfeldr-Kuhfeld_temp
_SP=/software/spack/v1.1.1/opt/spack/linux-zen4
export PATH=$_SP/openmpi-5.0.8-ruqvw4kfg3ftlc6nmn4gm6w7x2rwrtfl/bin:$PATH
source $_PILOT/oprm1-opes-activation/opt/gromacs-2025.0-plumed/bin/GMXRC
export PLUMED_PREFIX=$_PILOT/oprm1-opes-activation/opt/plumed-2.10.1-mpi
export PATH=$PLUMED_PREFIX/bin:$PATH
export LD_LIBRARY_PATH=$PLUMED_PREFIX/lib:${LD_LIBRARY_PATH:-}
export PLUMED_KERNEL=$PLUMED_PREFIX/lib/libplumedKernel.so
export PY=$_PILOT/miniforge3/envs/oprm1-opes/bin/python
# grompp with 2025.3: 2025.0's grompp mis-parses ParmEd's multi-line ff19SB [ cmaptypes ]
# ("Unknown atomtype found at position 2 in cmap type"); 2025.0 mdrun reads 2025.3 tprs.
export GMX_GROMPP=$_PILOT/spack-install/linux-zen4/gromacs-2025.3-45hxd2evh2mxcjqopqa56cr3atcbwlrv/bin/gmx
