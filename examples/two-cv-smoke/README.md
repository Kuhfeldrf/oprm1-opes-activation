# two-cv-smoke — verify the CV machinery in seconds, without an allocation

**What it shows.** Both activation collective variables resolve to the intended atoms,
respond to motion, and can be biased, through the same code path as production:
GROMACS ≥ 2025 native `-plumed` interface plus a PLUMED kernel. Both CVs use the
production definitions. CV1 is the Cα(Arg167)–Cα(Thr281) distance. CV2 is the NPxxYA
(334–339) Cα RMSD against the inactive reference, aligned on the TM1–TM5 scaffold.
The production `system_ref.pdb` (canonical human numbering) and dual-weight
`ref_inactive_npxxy.pdb` are used unchanged.

**What it is not.** The system is the receptor plus DAMGO (4,634 atoms) cut out of the
production membrane system, with byte-identical force-field parameters (ff19SB +
GAFF2), **in vacuum**, run for 1 ps. The membrane, water and ions are gone, so the
dynamics are not physical. That does not matter for what is being tested: whether the
CVs are wired correctly. The vacuum system carries the solute's +13 charge, and
grompp notes this; it is harmless here.

## Run

```bash
export PLUMED_KERNEL=/path/to/libplumedKernel.so   # PLUMED >= 2.9
GMX=gmx NT=4 ./run.sh                              # ~5 s on 4 CPU cores
```

## Pass condition (runbook §7.4)

`check.py` compares the output against `expected_output.json` and exits 0 only if all of
the following hold:

1. `COLVAR_smoke` is written for every stride (≥ 100 rows).
2. Both CVs start at the known values for these coordinates: CV1 = 1.1990 nm,
   CV2 = 0.3844 nm, ± 0.002 nm. A selector resolving to the wrong atom fails here,
   which is the silent-failure mode §7.1 warns about.
3. **Both** CVs change by more than 0.02 nm, and CV1 moves toward its 1.00 nm restraint.
4. **Both** bias columns are non-zero. If the CVs print but the biases are zero, PLUMED
   is computing without being coupled to the integrator.

Production adds a fifth condition, GPU offload in the GROMACS log. The CPU-only example
cannot test it; see `env/versions.md` for the production smoke test.

Reference output (ORCA, 4 threads, 5.6 s):

```
#! FIELDS time d_tm36 rmsd_npxxy r1.bias r2.bias
 0.000000 1.198951 0.384432 39.581435 34.015322
 1.000000 1.065319 0.305006 4.266593 11.026313
PASS  (all 8 checks)
```

Rebuild from the production system: `python examples/make_example.py systems/active examples/two-cv-smoke`.
