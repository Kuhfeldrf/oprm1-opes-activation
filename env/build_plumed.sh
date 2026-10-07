#!/usr/bin/env bash
#SBATCH --job-name=plumed_build
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --output=logs/plumed_build_%j.log
# Build PLUMED 2.10.1 from source with all modules (OPES is default-off in the
# conda-forge build). Serial (no MPI) kernel, matching the ~mpi GROMACS 2025.3
# spack build; GROMACS loads it at runtime via PLUMED_KERNEL (-plumed flag).
set -euo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
PREFIX=$ROOT/opt/plumed-2.10.1
source /etc/profile.d/z00_lmod.sh
cd "$ROOT/vendor"
rm -rf plumed-2.10.1 && tar xzf plumed-src-2.10.1.tgz && cd plumed-2.10.1
gcc --version | head -1
./configure --prefix="$PREFIX" --enable-modules=all --disable-mpi \
    --disable-external-blas --disable-external-lapack CXX=g++ CC=gcc \
    CXXFLAGS="-O3 -march=znver4" > "$ROOT/logs/plumed_configure.log" 2>&1
make -j "${SLURM_CPUS_PER_TASK:-16}" > "$ROOT/logs/plumed_make.log" 2>&1
make install > "$ROOT/logs/plumed_install.log" 2>&1
"$PREFIX/bin/plumed" info --long-version
"$PREFIX/bin/plumed" config module opes
ls -l "$PREFIX/lib/libplumedKernel.so"
echo BUILD_OK
