#!/usr/bin/env bash
#SBATCH --job-name=gmx_patched_build
#SBATCH --partition=normal
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/gromacs_patched_build_%j.log
# GROMACS 2025.0 patched with PLUMED 2.10.1 (runtime mode), MPI + CUDA.
#
# WHY: GROMACS 2025's NATIVE -plumed interface does not pass the multi-simulation
# communicator to PLUMED, so OPES WALKERS_MPI under -multidir silently runs N
# independent biases (verified 2026-10-07: walkers reported different opes.nker/zed).
# The classic patch passes it. 2025.0 is the newest version PLUMED 2.10.1 can patch.
set -eo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
PREFIX=$ROOT/opt/gromacs-2025.0-plumed
SP=/software/spack/v1.1.1/opt/spack/linux-zen4
GCC=/software/spack/gcc/opt/spack/linux-zen4/gcc-13.4.0-jzaajeossonz4fmu7l376bc7mzftzxam
export PATH=$GCC/bin:$SP/openmpi-5.0.8-ruqvw4kfg3ftlc6nmn4gm6w7x2rwrtfl/bin:$SP/cuda-12.9.0-vawg4hckgzbpqpwebuij3cxcl3dejzje/bin:$PATH
source /etc/profile.d/z00_lmod.sh || true
module load cmake/3.31.9-gcc-13.4.0 || true
source $ROOT/env/plumed_env_mpi.sh       # plumed (MPI build) for `plumed patch`
cd $ROOT/vendor
[ -f gromacs-2025.0.tar.gz ] || curl -sSfLO https://ftp.gromacs.org/gromacs/gromacs-2025.0.tar.gz
rm -rf gromacs-2025.0 && tar xzf gromacs-2025.0.tar.gz && cd gromacs-2025.0
plumed patch -p -e gromacs-2025.0 --runtime
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$PREFIX -DGMX_MPI=ON -DGMX_GPU=CUDA \
  -DCMAKE_CUDA_ARCHITECTURES="80;89" -DGMX_BUILD_OWN_FFTW=ON -DGMX_SIMD=AVX_512 \
  -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ -DGMX_OPENMP=ON \
  -DCMAKE_BUILD_TYPE=Release > $ROOT/logs/gromacs_patched_cmake.log 2>&1
make -j ${SLURM_CPUS_PER_TASK:-32} > $ROOT/logs/gromacs_patched_make.log 2>&1
make install > $ROOT/logs/gromacs_patched_install.log 2>&1
$PREFIX/bin/gmx_mpi --version | grep -E "GROMACS version|MPI library|GPU support|PLUMED"
$PREFIX/bin/gmx_mpi mdrun -h 2>&1 | grep -i -A1 plumed | head -4
echo BUILD_OK
