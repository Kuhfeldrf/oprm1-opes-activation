# Software versions (recorded 2026-10-07 on ORCA)

| Component | Version / build | Where |
|---|---|---|
| GROMACS | 2025.3-spack, mixed precision, CUDA 12.9 (sm_80, sm_89), thread-MPI, OpenMP, AVX_512, cuFFT | `/scratch/kuhfeldr-Kuhfeld_temp/spack-install/linux-zen4/gromacs-2025.3-45hxd2e…` (spack hash `45hxd2e`, built by the parent repo) |
| GROMACS (MPI) | 2025.3 +mpi +cuda, OpenMPI 5.0.8 — **building** (`env/build_gromacs_mpi.sh`), needed only for multi-walker OPES (`WALKERS_MPI`) | `opt/spack-env-mpi` |
| PLUMED | **2.10.1, built from source** with `--enable-modules=all` (OPES on), no MPI, internal BLAS/LAPACK, gcc 11.5 `-O3 -march=znver4` (`env/build_plumed.sh`) | `opt/plumed-2.10.1` |
| PLUMED ↔ GROMACS | GROMACS 2025 native interface (`gmx mdrun -plumed`), kernel loaded at runtime via `PLUMED_KERNEL`; **no patching** | `env/plumed_env.sh` |
| CUDA | 12.9.0 (default module) | |
| AmberTools | 26.0 (tleap, antechamber, sqm, packmol-memgen) | conda env `mor-pilot` |
| pdb2pqr / propka | 3.6.1 / 3.5.1 | conda env `mor-pilot` |
| ParmEd | 4.3.1 | conda env `mor-pilot` |
| gemmi / MDAnalysis / numpy / OpenMM | 0.7.5 / 2.10.0 / 2.5.3 / 8.6.1 | conda env `oprm1-opes` |

## Why PLUMED was built from source

conda-forge `plumed 2.10.1` installs and loads, but `plumed config module opes` reports
**`opes off (default-off)`** — OPES_METAD is not compiled in. Using it would have failed
only at production time. The source build reports `opes on`.

## §5 / §7.4 smoke test (2026-10-07, job 197759, 1× L40S)

Run on the parent repo's equilibrated DAMGO/μOR membrane system (102,002 atoms) as an
end-to-end coupling check before our own systems exist. Pass on all four conditions:

| Condition | Result |
|---|---|
| `COLVAR_smoke` written | yes |
| both CVs change | d_tm36 1.2616 → 1.2211 nm (target 1.20); rmsd_npxxy 0.0001 → 0.0664 nm (target 0.20) |
| both bias columns non-zero | r1.bias 0.95 → 0.11; r2.bias 9.99 → 4.47 kJ/mol |
| GPU offload | PP + PME on GPU; 195 ns/day (2 fs, 1000 steps, short-run estimate) |

Notes: `-update gpu` is unavailable because OPC is a 4-site (virtual-site) water model; update
and constraints run on the CPU. PLUMED forces this anyway.
