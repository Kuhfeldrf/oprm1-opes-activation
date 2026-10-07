#!/usr/bin/env bash
#SBATCH --job-name=plumed_mpi_build
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --output=logs/plumed_mpi_build_%j.log
# PLUMED 2.10.1, all modules (OPES on), built against the spack OpenMPI 5.0.8 that
# gmx_mpi (spack hash 47lswoc) links, so WALKERS_MPI works for multi-walker OPES.
set -eo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
PREFIX=$ROOT/opt/plumed-2.10.1-mpi
MPI=/software/spack/v1.1.1/opt/spack/linux-zen4/openmpi-5.0.8-ruqvw4kfg3ftlc6nmn4gm6w7x2rwrtfl
export PATH=$MPI/bin:$PATH
cd "$ROOT/vendor"
rm -rf plumed-2.10.1-mpi && mkdir plumed-2.10.1-mpi
tar xzf plumed-src-2.10.1.tgz -C plumed-2.10.1-mpi --strip-components=1
cd plumed-2.10.1-mpi
./configure --prefix="$PREFIX" --enable-modules=all --enable-mpi \
    --disable-external-blas --disable-external-lapack CXX=mpicxx CC=mpicc \
    CXXFLAGS="-O3 -march=znver4" > "$ROOT/logs/plumed_mpi_configure.log" 2>&1
grep -i "mpi" config.log | grep -i -E "found|enable|__PLUMED_HAS_MPI" | head -5 || true
make -j "${SLURM_CPUS_PER_TASK:-16}" > "$ROOT/logs/plumed_mpi_make.log" 2>&1
make install > "$ROOT/logs/plumed_mpi_install.log" 2>&1
export LD_LIBRARY_PATH=$PREFIX/lib:$LD_LIBRARY_PATH
"$PREFIX/bin/plumed" info --long-version
"$PREFIX/bin/plumed" config module opes
"$PREFIX/bin/plumed" config has mpi
echo BUILD_OK
