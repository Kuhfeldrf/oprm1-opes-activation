#!/usr/bin/env bash
#SBATCH --job-name=gmx_mpi_build
#SBATCH --partition=normal
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --output=logs/gromacs_mpi_build_%j.log
# GROMACS 2025.3 +mpi +cuda for multi-walker OPES (PLUMED WALKERS_MPI needs real
# MPI; the existing spack build is ~mpi / thread-MPI only). Same spack install
# tree as the parent repo so dependencies are reused.
set -eo pipefail
PILOT=/scratch/kuhfeldr-Kuhfeld_temp
ROOT=$PILOT/oprm1-opes-activation
SPACK_ENV=$ROOT/opt/spack-env-mpi
source /etc/profile.d/z00_lmod.sh
module load spack/v1.1.1
[ -d "$SPACK_ENV" ] || spack env create -d "$SPACK_ENV"
spack -e "$SPACK_ENV" config add "config:install_tree:root:$PILOT/spack-install"
spack -e "$SPACK_ENV" config add "config:build_jobs:32"
spack -e "$SPACK_ENV" add "gromacs@2025.3 +cuda cuda_arch=80,89 +mpi +openmp ^openmpi@5.0.8"
spack -e "$SPACK_ENV" concretize -f 2>&1 | tail -20
spack -e "$SPACK_ENV" install --fail-fast -j32
spack -e "$SPACK_ENV" find -lvd gromacs | head -5
P=$(spack -e "$SPACK_ENV" location -i gromacs)
echo "PREFIX $P"
ls "$P/bin"
echo BUILD_OK
