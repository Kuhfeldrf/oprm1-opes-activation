# Plausibility log

Runbook §2: after each stage, (1) does what I found match what the runbook assumed,
(2) is the next step still right, (3) what would failure look like here. Where the
runbook and reality disagree, reality wins and is recorded here.

---

## Day 1 — 2026-10-07

### Scope decision (§3): Option A

Agonist-bound (DAMGO) μOR, sampled from an active start **and** an inactive start.
The apo and antagonist systems are only built if the agonist system converges early.
This follows the runbook's recommendation, and nothing found today argues against it.

### Stage 2 — PLUMED/GROMACS availability (§5)

**Found:** ORCA has **no** GROMACS or PLUMED module (`module avail`, login shell). Taken
literally, that is outcome 3. However, the parent repository
(`/scratch/kuhfeldr-Kuhfeld_temp`, MOR_peptide_binding) already has a spack-built
**GROMACS 2025.3 + CUDA**. GROMACS ≥2025 ships a native PLUMED interface (`mdrun -plumed`),
which loads a PLUMED kernel at runtime, so the situation is effectively outcome 2 with no
patching step.

**Deviation 1 — conda-forge PLUMED lacks OPES.** `plumed config module opes` →
`opes off (default-off)`. Using it would have failed at the first `OPES_METAD` line. PLUMED
2.10.1 was rebuilt from source with `--enable-modules=all` (`env/build_plumed.sh`), and
that build reports `opes on`.

**Deviation 2 — no real MPI in the GROMACS build.** It uses thread-MPI, so `WALKERS_MPI`
multi-walker OPES is impossible with it. A `+mpi` GROMACS build is running
(`env/build_gromacs_mpi.sh`); it also needs an MPI-enabled PLUMED kernel. If it fails, the
fallback is single-walker runs, with the reduced sampling stated in the verdict.

**Smoke test (§7.4):** passed on all four conditions. It ran on the parent repo's
equilibrated DAMGO system, not ours, because ours do not exist yet; see `env/versions.md`.
It will be repeated on our own system and committed as `examples/two-cv-smoke/`.

**What failure would have looked like:** bias columns flat at zero (kernel not loaded),
or an unresolved OPES action. Neither happened. The kernel-load banner appears in the log.

**Silent-failure modes hit and fixed during the smoke test:**
- PLUMED parses `key=value` tokens in PDB `REMARK` lines as arguments. A remark
  containing "occupancy=align" aborted initialisation with `Error converting align`. The
  reference builder now writes no `=` in remarks.
- OPC (4-site water) has virtual sites, so `-update gpu` is not allowed. This is not a
  misconfiguration; update runs on the CPU, as PLUMED requires anyway.
- **§7.1 numbering trap confirmed in the wild:** the parent system's topology numbers the
  receptor from 1, so canonical R167 is residue 101 (constant offset −66). `@CA-167` would
  have silently selected the wrong atom. Our build must restore canonical numbering in
  `system_ref.pdb`, and the build asserts it.

### Stage 1 — Structure selection (§4) and numbering (§6.5)

`cv/assert_numbering.py` (exit 0; full output in `structures/assert_numbering.out`). All 8
entries map to P35372 with zero offset. D166/R167/Y168/T281/N334–A339 are correct in every
receptor chain, and every author chain ID matches §6.5.

**Runbook vs reality:**

| Entry | Runbook said | Found |
|---|---|---|
| 9MQH | 72–354 modelled, gap 225–226, M266L/K271R | ✔ as stated |
| 9MQJ | completeness not verified | 67–354, gap 225–226, M266L/K271R. Title: "Locally-refined … with Nb6M, NabFab, and isoquinuclidine compound #020_E1". GPCRdb annotates #020_E1 as **Antagonist**, consistent with §4.2. |
| 9MQI | completeness not verified | 67–354, gap 225–226, M266L/K271R; partners Nb (C) + Fab (H, L) present |
| 9PXU | T281/NPxxY not verified | **receptor 68–354 fully resolved, no gaps**; T281 and NPxxY present. **Also carries M266L, K271R** (not listed in §4.2). BRIL is in auth chain R at negative numbers (−108…−8). **A Na⁺ ion is modelled** (sodium site). Partners: Nb (C), K, Fab (H, L). |
| 10TM | R,F chains; F158W | ✔. Two receptor copies (R, F), each with DAMGO (S, H). **No G protein or nanobody in the model.** Resolved 66–347. CLR modelled. |
| 9PPQ | state not confirmed | Receptor 67–352, no gaps, no SEQADV. **GPCRdb: state Active, ligand 0505 Agonist** (queried 2026-10-07), and on deposited coordinates it clusters with the active structures (below). §15 gap closed. |
| 8F7Q | numbering inferred | DBREF read directly: P35372 2–388, canonical. **Carries F158W (SEQADV "conflict")**, which the runbook missed — revert if 8F7Q is ever used. |
| 8Y72 | R chain, cross-check | ✔. No receptor SEQADV. |

**DAMGO is not "protein force field directly" (§8).** It contains D-Ala (DAL),
N-methyl-Phe (MEA) and a Gly-ol C-terminus (ETA). The parent repo parameterised it with
GAFF2/AM1-BCC (charge +1), and that is reused here (see force-field decision below).

### §6.4 cross-check on deposited coordinates — DECISION POINT PASSED

`cv/measure_deposited.py`. CV2 is computed with the dual-weight scheme of §7.2: superpose
TM1–TM5 scaffold Cα (134 atoms), then take the RMSD of NPxxYA Cα with no refit. Units Å.

| PDB | chain | state | CV1 R167–T281 Cα | CV2 vs 9MQH | CV2 vs 9PXU | R167 CZ–T281 OG1 |
|---|---|---|---|---|---|---|
| 9MQH | A | inactive | 6.75 | 0.00 | 0.63 | 5.49 |
| 9MQJ | A | inactive | 6.68 | 0.90 | 1.15 | 4.38 |
| 9MQI | A | inactive | 6.65 | 0.83 | 1.08 | 4.36 |
| 9PXU | R | inactive | 6.58 | 0.63 | 0.00 | 4.61 |
| 10TM | R/F | active | 12.24 | 4.09 | 3.78 | 12.23 |
| 8F7Q | R | active | 13.08 | 3.92 | 3.63 | 12.68 |
| 8F7Q | M | active | 13.18 | 4.08 | 3.74 | 12.62 |
| 8Y72 | R | active | 13.03 | 4.41 | 4.10 | 12.68 |
| 9PPQ | R | ? | 13.00 | 4.34 | 4.06 | 12.58 |

- Endpoints separate by **~5.5–6.6 Å on CV1**, well above the 2 Å failure threshold
  (§4.5), and by ~3 Å on CV2.
- The three independent active structures (10TM, 8F7Q, 8Y72) cluster within 1 Å on CV1
  and 0.5 Å on CV2. The coordinates report activation state cleanly, so §6.6 is not
  triggered.
- Scaffold superposition RMSD is 0.9–1.4 Å across all entries. The TM1–TM5 alignment set
  is stable across states, as §7.2 requires.
- 10TM, without a G protein, has a slightly smaller TM6 opening (12.2 Å) than the
  Gi-bound structures (13.0–13.2 Å). It is still fully on the active side.
- 9PPQ is geometrically active-like, and GPCRdb independently annotates it Active. It is
  not used as a start.
- GPCRdb annotates naloxone in 9PXU as "Agonist (partial)", an odd label for a
  prototypical MOR antagonist. Irrelevant here because naloxone is deleted, but if 9PXU
  is ever used as an "antagonist-bound" control (Option B), cite the pharmacology, not
  the GPCRdb tag.
- In every inactive structure, R167–T281 sidechain distance is 4.4–5.5 Å, consistent
  with the R3.50–T6.34 hydrogen bond being formed or nearly formed at 3.9 Å resolution.

### Start structure choice for Option A

- **Active start: 10TM chain R + DAMGO chain S.** Revert F158W.
- **Inactive start: 9PXU chain R**, deviating from the runbook's 9MQH for this purpose.
  The reasons: 3.4 Å vs 3.9 Å; no ECL2 gap, so no loop modelling; the sodium site is
  resolved, which §8 says to keep; and naloxone sits in the same orthosteric pocket that
  DAMGO will occupy. Strip BRIL, Nb, Fab and naloxone, and revert L266M and R271K.
  **DAMGO is placed by superposing 10TM onto 9PXU on the TM1–TM5 scaffold**, then
  checked for clashes. It is the same ligand in both starts, which start-independence
  requires. 9MQH stays the apo-inactive receptor if Option B is ever run.
- **Common modelled range: 69–346 for both** (10TM ends at 347; 9PXU starts at 68), with
  residue 68 → ACE and 347 → NME caps. This gives the two starts identical sequence and
  composition.
- No apo-active structure exists (§4.5); the active side is necessarily agonist-bound.
  This asymmetry is stated, not worked around.

### Force field decision (§8)

**Amber ff19SB + Lipid21 + OPC + GAFF2 (DAMGO)**, the runbook's documented alternative, is
used instead of CHARMM36m. Reasons:
1. A scripted, validated build exists in the parent repo: packmol-memgen → tleap →
   ParmEd → GROMACS, with per-species charge and atom checks. DAMGO in it held its
   orthosteric pose for 50 ns.
2. A CHARMM36m build would go through CHARMM-GUI, a web service that cannot be scripted
   from the cluster.

Every system and replica uses this one force field.

### Stage 4 — System preparation (§8), Day 1

**Pack:** packmol-memgen, POPC:CHL1 7:3 (30 mol% cholesterol), 150 mM NaCl, xy fixed at
95 Å for both starts, with the parent repo's calibrated settings (`--pbc`,
`apl_offset 1.15`, `nloop 200/40`). As with the parent's validated DAMGO pack, packmol
stopped at its iteration ceiling short of tolerance. Minimisation and restrained NPT are
what resolve this; whether they succeed is checked in the equilibration diagnostics.

**Silent traps hit and handled in `systems/build_topology.py`:**
1. **Charge.** packmol-memgen's charge estimate (+12) **ignored DAMGO's +1**. The
   true solute charge from tleap is **+13**, and a naive build is +1 non-neutral
   (caught by the neutrality assertion). The builder now measures the solute charge
   with tleap and sets Cl⁻ from it.
2. **DAMGO mol2 charges summed to +1.0020**, not +1, inherited from the parent; the
   parent's runs carried a 0.002 e net charge. Renormalised to exactly +1
   (−2.7×10⁻⁵ e/atom, `structures/params/renormalize_charge.py`; original kept).
3. **ff19SB's NME methyl carbon is named `C`**, not `CH3` (ACE's is `CH3`), and
   tleap aborts on it. Fixed in `clean.py` and normalised in the builder.
4. **Renumbering.** packmol-memgen renumbers the receptor from 1 (ACE68 → 1), and so
   does tleap. `system_ref.pdb` is written with **canonical numbers restored** (+67) and
   asserted against D166/R167/Y168/T281/N334–A339. PLUMED MOLINFO reports chain R =
   residues 68–347.
5. **Topology order** is receptor (1–4561), then DAMGO (4562–4634), then lipids,
   waters and ions. tleap's `combine` would otherwise append the ligand after the
   solvent.

**Identical composition (both starts):** 194 POPC, 84 CHL, 11,903 OPC waters, 20 Na⁺,
33 Cl⁻, 1 DAMGO; 84,511 atoms; box 95 × 95 × 86 Å. Adjustments to reach it:
- active: −12 waters, +1 Cl⁻ from water;
- inactive: −1 POPC from the fuller (upper) leaflet, farthest from the receptor, plus
  Na⁺/Cl⁻ from bulk waters;
- in both, molecules were removed from or converted in bulk, farthest from the receptor.

The system is ~85k atoms instead of the runbook's ~150k target. It is smaller because
the box is sized from the receptor (≥20 Å water beyond protein), and smaller is faster.

**HMR:** hydrogen mass 3.024 amu on solute and lipids (not water), which allows a 4 fs
timestep in production. This is a deviation from a plain 2 fs protocol, chosen for
throughput. Heating and the first three NPT stages still run at 2 fs.

**PLUMED selector validation (§7.1),** `plumed driver --igro systems/active/system.gro`:
d_tm36 = 1.2236 nm, matching the deposited 10TM value of 12.24 Å, so `@CA-167`/`@CA-281`
resolve to the right atoms. d_hbond = 1.2235 nm (deposited 12.23 Å). d_salt =
0.437 nm, na_site = 0 (no Na⁺ at D2.50 in active, as expected), lig_cont = 349.

### Multi-walker OPES machinery test — a silent failure caught (Day 1)

Test: 2 walkers, `WALKERS_MPI`, `gmx_mpi -multidir` (spack GROMACS 2025.3 +mpi with MPI
PLUMED kernel, **native** `-plumed` interface), 20 ps from the active NVT frame.

- It ran: GPU-aware MPI, one A30 per walker, **~360 ns/day per walker at 4 fs with PLUMED
  attached**. PLUMED accounts for 27% of wall time.
- **But the bias was not shared.** Walkers sharing one OPES bias must report identical
  `opes.nker`, `opes.zed` and `opes.rct`. They reported nker 32 vs 30 and zed 0.196 vs
  0.261. PLUMED's log says only "WALKERS_MPI: if multiple replicas are present…", with no
  walker count. GROMACS 2025's native PLUMED interface does not hand PLUMED the
  multi-simulation communicator, so each walker silently ran its own independent OPES.
  No error or warning is raised. A production "4-walker" run would have been 4 unshared
  biases, converging ~4× slower than planned, with nothing in the output saying so.
- **Response:** build GROMACS **2025.0 patched with PLUMED 2.10.1** (classic
  `plumed patch --runtime`, which wires the multisim communicator),
  `env/build_gromacs_plumed_patched.sh`. 2025.0 is the newest version PLUMED 2.10.1
  can patch. The unbiased endpoint legs keep the native-interface 2025.3 build, which
  is correct for single-simulation runs.
- **Acceptance test for the patched build:** the same 2-walker test must show
  identical `opes.nker`/`opes.zed`/`opes.rct` in both walkers' COLVAR at every line.

**Patched build: acceptance test passed.** GROMACS 2025.0 + PLUMED 2.10.1 patch
(`env/plumed_env_patched.sh`). PLUMED now logs "using multiple walkers, number of
walkers: 2". `opes.rct`, `opes.zed` and `opes.nker` are **identical in both walkers on
101/101 lines**, while d_tm36 differs on 100/101 lines (independent configurations, one
shared bias). Throughput is ~400 ns/day per walker (A30, 4 fs, PLUMED attached).

Two more wrinkles:
- **2025.0's grompp mis-parses ParmEd's multi-line ff19SB `[ cmaptypes ]`** ("Unknown
  atomtype found at position 2 in cmap type"); 2025.3's does not. Workaround: grompp
  with 2025.3 (`$GMX_GROMPP`), mdrun with patched 2025.0. 2025.0 reads 2025.3 tprs
  (verified: `gmx_mpi dump` → natoms 84511).
- In multi-sim mode PLUMED writes `COLVAR.<k>`, and only walker 0 writes
  KERNELS/STATE. Restarts point every walker at `w0/STATE`.

### Bilayer embedding check (user request, Day 1)

The receptor is embedded exactly as in `MOR_demo_analysis` and the parent ORCA runs:
packmol-memgen, ff19SB + Lipid21 + OPC, POPC:CHL1 7:3, apl_offset 1.15, nloop 200/40,
`--pbc`, OPM orientation, 150 mM NaCl. `analysis/membrane_qc.py` checks per frame:
P–P thickness, area per lipid (receptor footprint subtracted), receptor centring and
tilt, lipid atoms inside the TM bundle, water in the hydrophobic core, and TM helicity.

| active frame | P–P (Å) | APL (Å²) | lipids in bundle | core water | TM helix frac |
|---|---|---|---|---|---|
| as built | 38.1 | 58.7 | 0 | 0 | 0.88 |
| npt_1000 (1 ns) | 40.7 | — | 0 | 7 | 0.90 |
| npt_200 (3 ns) | 41.1 | 53.0 | 0 | 0 | 0.90 |

Leaflets are balanced at 139/139 (POPC+CHL). The receptor TM core sits 3 Å below the
midplane, and the bundle axis is tilted ~22°. APL is still falling, as intended for
the deliberately under-packed build. It must plateau in the unbiased leg before
biasing; the QC is rerun over both 100 ns legs at the §6.4 decision point.

**4 fs + Lipid21:** grompp warns that the oleoyl C=C bond's period (20 fs) is exactly
5 × dt. HMR does not touch heavy–heavy bonds; Lipid21 is validated with HMR at 4 fs,
and Verlet is stable well past this ratio. `equilibrate.sh` accepts **only** this
"oscillational period" warning and fails on any other; there is no blanket `-maxwarn`.

### Unbiased legs crashed at 4 fs — methionine methyls, fixed by LINCS accuracy (Day 1)

**What happened.** Both systems completed all six restrained NPT stages, then aborted
in unrestrained production with "Too many LINCS warnings (1000)": active at 28 ps,
inactive at 60 ps. **Every failing atom was a methionine CE** (the S–CH₃ methyl):
Met74, 92, 101, 132, 153, 163, 205, 245, 266 and 283, the same set in both systems.
Constraint deviations were tiny (rms ~4×10⁻⁵), but CE–H bonds rotated > 30° per step.
That is ~25× the thermal angular speed of an HMR methyl, so it was an integration
artefact, not physics. Warnings began in the first 4 fs NPT stage. Solute restraints
kept them under the per-run abort limit until production.

**Test** (`env/test_timestep.sh`; 100 ps each from the same equilibrated `npt_10` state;
warnings counted in stderr, since GROMACS does not write them to md.log):

| setup | LINCS warnings / 100 ps | ns/day (A30) | T Solute/MEMB/SOLV (K) |
|---|---|---|---|
| 4 fs, lincs-order 4 / iter 1 (as run) | aborted at 28 ps | — | — |
| **4 fs, lincs-order 6 / iter 2** | **0** | **393** | 310.7 / 310.8 / 308.7 |
| 3 fs, order 4 / iter 1 | 0 | 359 | 311.3 / 311.2 / 309.4 |
| 2 fs, order 4 / iter 1 | 0 | 255 | 310.0 / 310.1 / 309.5 |

**Decision:** 4 fs with **lincs-order 6, lincs-iter 2** for every stage and for OPES
(`systems/mdp/_common.mdp`). Default LINCS is too coarse for coupled CH₃ constraints with
3 amu hydrogens at 4 fs, and the GROMACS manual recommends higher order/iterations for
large time steps. The equilibrated `npt_10` states are kept: deviations stayed ~10⁻⁴
throughout and the bilayer QC was clean. The failed production legs are kept under
`systems/<start>/failed_prod_lincs4/` and are **not** used for anything. The production
legs are watched continuously for any recurrence.

### Repository visibility

2026-10-07: the repository owner approved pushing to the **public** GitHub repo
(`Kuhfeldrf/oprm1-opes-activation`). This supersedes the runbook's §12.5 private-first
instruction. Pushes follow each stage gate from here on.

## Day 1 → 2 — §6.4 / §8 decision point: PASSED (2026-10-07, 23:30)

1. **Does it match what the runbook assumed?** Yes. Both endpoints are stable for 100 ns
   unbiased. The active state held without its G protein (CV1 1.20 nm, drift −1.3 Å;
   threshold 3 Å). The inactive state kept the R3.50–T6.34 contact and its Na⁺. DAMGO
   stayed salt-bridged to D149 in both. Full table in `cv/cv_definition.md`.
2. **Is the next step still right?** Yes. The 2D space separates the states by 10 SD (CV1)
   and 9 SD (CV2), and the deposited active structures cluster (Day 1). No §8 failure
   signature: no CV1 drift > 3 Å, CV2 settled (trend ≤ 0.02 nm/100 ns), helicity 0.89–0.90.
3. **What would failure look like?** Endpoint overlap, drift, bilayer collapse, or lipid
   entering the bundle. None are seen. The bilayer QC over both legs is P–P 43.6 Å, APL
   settled at ~48 Å² (receptor footprint removed), 0 lipid atoms in the bundle in every
   frame, and core water ≤ 8.

Note: the active leg's NPxxY (CV2) sits at 0.32 nm against the equilibrated inactive
reference, below the deposited-structure gap (~0.38 nm). The active TM7 relaxes
partway without a G protein; TM6 does not. Recorded now because it bears on whether the
active basin found by OPES will match the deposited active geometry.

## Day 2 — OPES launch, first parameter revision (2026-10-08 00:05)

**prod1** (4 walkers × 2 starts; SIGMA 0.025/0.027 nm; walls at CV1 0.53–1.50 nm and
CV2 ≤ 0.47 nm; KAPPA 2000; BARRIER 50 kJ/mol). Checked at 10 ns/walker (40 ns per
start):
- The shared bias was building normally. nker reached 1100 (active start) and 486
  (inactive start).
- The active start was already exploring toward inactive: one walker reached CV1
  0.79 nm. There were no committed crossings yet, as expected at this stage.
- **But the walkers were pushed through the walls.** In the active start, walker 2 spent
  44% of frames beyond the CV2 upper wall (max 0.58 nm vs wall 0.47), and another
  walker reached CV1 1.595 nm (wall 1.50). The inactive start went to CV1 0.49 nm
  (wall 0.53).
- **Diagnosis:** the walls were too soft relative to BARRIER. At KAPPA 2000 kJ/mol/nm², a
  wall reaches 50 kJ/mol only after 0.16 nm of penetration, about 6 SIGMA. OPES fills
  the basin up to BARRIER and then simply climbs the soft wall.
- **Fix:** KAPPA 50000 (45 kJ/mol at 0.03 nm penetration, about 1 SIGMA). Wall positions are
  unchanged, since they come from the measured endpoints. Per §9, the parameters were
  revised instead of spending more wall-clock. prod1 is kept as
  `runs/prod1_softwalls/` and is **not** used in any result.
- **prod2** relaunched from scratch with the same walker starting frames.

## Day 2 — prod2 monitoring snapshot at ~71 ns/walker (2026-10-08, 6.5 h wall-clock)

284 ns aggregate per start; 260 ns/day per walker, so 375 ns/walker lands in about
29 h more. Snapshot: `results/monitor_prod2_70ns.json`.

- **The stiff walls hold.** Wall contact is ≤ 3.2% of frames in the active start. It is
  ≤ 1.5% in the inactive start except walker 1 (8.6%, lower CV1 wall at 0.53 nm; it
  sits at TM3–TM6 0.51–0.55 nm, tighter than the deposited inactive state).
- **Ligand and Na⁺ are stable.** DAMGO contacts are ≥ 304 and d_salt ≤ 0.50 nm in all 8
  walkers. Na⁺ stays seated in the inactive start (0.87–0.99) and absent in the active
  start.
- **Bias.** The nker growth is sublinear: active 1216→2496 and inactive 657→1104
  across the quarters, so deposition is slowing. rct is still small (−2…+2 kJ/mol), so
  the bias is far from filling the 50 kJ/mol BARRIER.
- **Crossings: 0 in each direction.** CV1 and CV2 move independently:
  - **CV1 crosses.** Active-start walker 0 spent 21% of its frames with TM3–TM6 inside
    the inactive CV1 range (minimum 0.537 nm). Inactive-start walkers reached 1.05–1.08 nm,
    near the active box edge (1.089).
  - **CV2 does not follow.**
    - In the active start, NPxxY RMSD never drops below 0.11 nm (inactive box ≤ 0.118).
      Walker 0, with TM6 closed, has NPxxY at 0.42 nm, even further from inactive.
    - In the inactive start, NPxxY never exceeds 0.27 nm (active box ≥ 0.265), and only
      when CV1 is low.
  - The CV1–CV2 correlation per walker is −0.46…+0.21. The walkers are filling the
    off-diagonal corners (TM6 closed with NPxxY out; TM6 open with NPxxY in), not the
    coupled path.

Three questions:
1. **Does it match what the runbook assumed?** Partly. Transitions in CV1 alone are
   fast, as expected for a biased distance. But the TM7/NPxxY rearrangement is
   hysteretic within the runtime: each start keeps its own NPxxY state. This is the
   early form of the §8 "start-dependent FES / missing slow DOF" signature, not yet a
   verdict. The bias is far below BARRIER, and OPES fills orthogonal directions before
   committing to a coupled path.
2. **Is the next step still right?** Yes, continue to the planned 375 ns/walker; there is no
   reason to stop (no wall pinning, no ligand loss, no fold loss). If CV2 is still
   start-locked at the end, the verdict records that NPxxY/TM7 (and possibly the Na⁺
   site, which correlates with it in the inactive start) needs an extra or better CV.
   I will not change parameters mid-run.
3. **What would failure look like at the end?**
   - Zero committed crossings in either direction.
   - Two FESs that differ by > 2 kcal/mol in the CV2 direction, each with a minimum at
     its own starting NPxxY value.
   - rct still rising in the last quarter.

## Day 2 — prod2 snapshot at ~182 ns/walker (2026-10-08 19:05, 19 h wall-clock)

Snapshot: `results/monitor_prod2_182ns.json`. Current speed is 259 ns/day per walker
(measured over 2 min). The average since launch is 230 ns/day, because of a slower
stretch between 15:00 and 19:00 with no error in the logs. 375 ns/walker lands in
about 18 h more.

- **Still 0 committed crossings in each direction.** Each CV now reaches the other
  state on its own, but never together:
  - In the active start, walker 1 brought NPxxY to 0.068 nm (inactive box ≤ 0.118)
    while TM3–TM6 stayed ≥ 1.07 nm. Walkers 0 and 2 closed TM3–TM6 to 0.52–0.56 nm
    with NPxxY ≥ 0.16 nm.
  - In the inactive start, walker 3 reached NPxxY 0.367 nm and CV1 1.083 nm, 0.006 nm
    short of the active box. These were not simultaneous.
- **Bias.** nker growth keeps slowing: active +730 and inactive +176 in the last quarter.
  rct is flat at −0.4…+1.5 kJ/mol, so the deposited bias stays far below
  BARRIER 50 kJ/mol. OPES is refining the explored region, not pushing outward.
- **Walls.** Walker 1 of the inactive start is at the lower CV1 wall 8.9% of the time.
  It prefers TM3–TM6 at about 0.50–0.51 nm, so the FES for CV1 < 0.55 nm is
  wall-limited and will be masked in the verdict. All other walkers are ≤ 2.9%.
- **Ligand, salt bridge, Na⁺, fold.** Unchanged and stable.

Three questions:
1. **Does it match what the runbook assumed?** No, not yet. The two-step pattern
   persists: TM6 and NPxxY each move separately, and the coupled activation transition
   has not happened in 730 ns per start.
2. **Is the next step still right?** Yes, finish the planned 375 ns/walker and then
   judge. Nothing is breaking (§8), so stopping early would only lose data.
3. **What would failure look like?** The same picture at the end: off-diagonal
   sampling only, and FES minima at each start’s own state.

## Day 3 — overnight throughput stall (2026-10-09 06:45)

Both OPES jobs (orcaga14, orcaga15) slowed about 20× at the same times:
- partly from 17:00 on Oct 8;
- continuously from 21:35 Oct 8 to 05:50 Oct 9 (about 0.4 ns per 45 min instead of 8).

Both recovered between 05:50 and 06:35; the current speed is 230 ns/day per walker.
The checkpoint timestamps in `opes.log` show the stall.

There are no errors, LINCS warnings or energy anomalies, and constraint RMSD is
unchanged. The stall hit two separate nodes simultaneously, so the shared resource is
the likely cause: `/scratch` is an NFS mount at 94% full. PLUMED appends COLVAR and
KERNELS every 500 steps on every walker, which is sensitive to NFS latency.

Effect: lost wall-clock only, not data. Position at 06:45 is 204 ns (active) and 211 ns
(inactive) per walker, about 17 h from completion if there are no further stalls. The
job time limit (mdrun -maxh 70) is not at risk.

If it recurs: restart from checkpoint with the run directory on node-local disk,
copying back at checkpoints.

## Day 3 — interim FES at ~205 ns/walker (2026-10-09 morning)

Run with `sbatch analysis/fes_job.sh runs/prod2 results/interim_fes_prod2`:
- grid CV1 0.48–1.56 nm, CV2 0–0.52 nm;
- the measured basin boxes;
- reweighting with the total bias.

The pipeline runs end to end. Its first output was all NaN ΔG, which is correct rather
than a bug: each starts FES has zero sampled grid points in the other starts basin
box (`basins_sampled`: active start {inactive 0.0, active 1.0}; inactive start the
reverse). `fes.py` now reports basin coverage explicitly. The two maps share only 12
grid points (out of about 6000 sampled by either), so the start-difference RMS (0.16
kcal/mol on those 12 points) means nothing yet.

Shape of the maps (`results/interim_fes_prod2/figures/fes_both_starts.png`):
- **Active start.**
  - Global minimum at CV1 1.25, CV2 0.32 nm.
  - A shoulder at CV1 ≈ 0.93 nm (about 3–4 kcal/mol).
  - A distinct secondary minimum at **CV1 0.67, CV2 0.42 nm** (about 5–6 kcal/mol):
    TM6 closed to inactive-like distance, NPxxY *more* displaced than active.
  - The low-CV2 region (< 0.14 nm) is unsampled at < 12 kcal/mol.
- **Inactive start.** A single basin at CV1 0.64, CV2 0.05 nm. Inside 12 kcal/mol it
  reaches only CV1 ≤ 0.96 nm and CV2 ≤ 0.30 nm.

Reading: TM6 closure from the active state is reversible within the CVs. The NPxxY
return to its inactive conformation is not, and from the inactive side TM6 opening is
uphill by > 12 kcal/mol in this map. The (CV1 0.67, CV2 0.42) state is a TM6-in /
TM7-out off-pathway (or intermediate) state, not the inactive basin. So far this
supports the hypothesis that a degree of freedom missing from the CVs (TM7 / NPxxY
rotamers / Y7.53 or Na⁺) gates the transition.

Final judgement waits for 375 ns/walker; parameters unchanged.

## Day 4 — §10 verdict (2026-10-10)

prod2 finished: 2 starts × 4 walkers × 375 ns = 3.0 μs biased, both jobs exiting 0. The
analysis is `sbatch analysis/final_analysis_job.sh runs/prod2 results/fes_prod2`
(FES on all data and on the second half, convergence, PLUMED cross-check, bilayer QC).

| §10 criterion | Measured | Verdict |
|---|---|---|
| 1 reversible crossing | 0 / 0 committed crossings in each start | FAIL |
| 2 start-independence | each start sampled 0% of the other basin; in the second half, the shared region below 8 kcal/mol (338 points) differs by 2.2 kcal/mol RMS, 5.9 max | FAIL |
| 3 no hysteresis | forward and reverse paths visit different regions (`figures/paths_overlay.png`) | FAIL |
| 4 convergence | walker SD active 2.2 → 2.2, inactive 0.75 → 1.9 kcal/mol (grows); Q3 vs Q4 RMS 2.1 / 0.9 kcal/mol | FAIL |

**Supporting checks:**
- DAMGO stayed bound (coordination ≥ 283; D149 salt bridge ≤ 0.50 nm).
- Walls were touched in ≤ 2.4% of frames, except inactive walker 1 at 7.2% (lower CV1
  wall; FES below about 0.55 nm masked).
- The PLUMED `REWEIGHT_BIAS` cross-check matches `fes.py` to 0.004–0.007 kcal/mol RMS,
  with identical minima.

**Plausibility checks on the analysis itself:**
- **The first membrane QC on the biased walkers flagged three walkers.** It showed
  "lipids in the bundle", tilt up to 90° and helicity down to 0.79. All of it was a **PBC
  artefact**: single-rank GROMACS wraps atoms into the box at neighbour-search steps, so
  the receptor is split once it diffuses across a box edge. An axis-free test also
  counted lipids at the "centre" from 250 ns on in inactive walker 3, again because of
  the split receptor. `membrane_qc.py` now makes the solute whole, centres it and wraps
  the rest by residue. After that fix, all 8 walkers show helicity ≥ 0.84, tilt ≤ 26°,
  and a single lipid atom in 1 of 608 frames. The CVs were never affected, because
  `WHOLEMOLECULES` is in `plumed_opes.dat`.
- **PLUMED 2.10 template differences from runbook §7.6:**
  - `READ VALUES=opes.bias` creates component `<label>.bias`;
  - `HISTOGRAM` with `LOGWEIGHTS` needs `NORMALIZATION=true`.
  Both are fixed in `cv/plumed/plumed_reweight.dat`.

Three questions:
1. **Does it match what the runbook assumed?** No. The runbook assumed the two
   coordinates capture activation. Instead:
   - TM6 (CV1) moves freely in both directions;
   - NPxxY (CV2) moves separately, sometimes against CV1;
   - the coupled transition never happened in 3 μs.
   Na⁺ at D2.50 never left in the inactive start and never bound in the active start.
   That is an unbiased slow variable differing between the starts, named as a suspect
   on Day 1.
2. **Is the next step still right?** Per §10, this is a clean negative and is reported as
   the finding, not tuned away. Option B was conditional on Option A converging early,
   so it is not started. The next experiments, out of this week's scope, are:
   - a third CV for Na⁺ occupancy and/or Y7.53 χ1;
   - raising BARRIER to 80–100 kJ/mol;
   - training a learned CV on these 3 μs.
3. **What would a false negative look like?** A force-field or setup error that freezes
   TM7 or sodium. Two checks argue against it:
   - the unbiased endpoints are stable and match deposited geometry;
   - both CVs individually traverse their full range.
   The remaining risk is the fixed protonation of D2.50, which propka predicts would
   change on activation. This is recorded as a model limitation.

## Day 5 — re-examining the diagnosis; follow-up runs (2026-10-11)

**Prompt.** I was asked to reconsider the Day 4 interpretation rather than build on it.
I had prepared an Na⁺-biased third-CV rerun. It was **not** launched.

**Hysteresis diagnostic.** I compared frames from the two starts at matched CV values,
using `results/fes_prod2/colvar_strided`:

| At CV1 0.55–0.72, CV2 0.18–0.32 nm | active start (w0, w2) | inactive start (all walkers) |
|---|---|---|
| Na⁺ at D2.50 | 0.00 | 0.91 ± 0.06 |
| R167 CZ – T281 OG1 | 0.78 ± 0.10 nm | 0.47 ± 0.08 nm |
| W295 χ1 | −1.24 | −1.26 |
| lig_cont / d_salt | 388 / 0.33 nm | 412 / 0.33 nm |

At CV1 0.72–0.95 nm, CV2 0.10–0.22 nm the picture is the same: Na⁺ 0.00 vs 0.89 and
side-chain contact 0.95 vs 0.59 nm.

**Active-start walkers that closed TM6 were still relaxing at the end.**
Means over 0–50 / 100–200 / 250–375 ns:

| Walker | CV2 (nm) | side-chain contact (nm) |
|---|---|---|
| w0 | 0.40 / 0.36 / 0.31 | 0.92 / 0.77 / 0.66 |
| w2 | 0.31 / 0.28 / 0.24 | 1.05 / 0.92 / 0.84 |

**Revised reading:**
1. **The hysteresis is real and involves at least two variables outside the CVs:** Na⁺
   occupancy and the R3.50–T6.34 side-chain contact. Na⁺ was over-weighted on Day 4.
2. **Insufficient time is not excluded.** The active-start deactivation was incomplete
   but still moving toward inactive in every variable.
3. **Na⁺ is not needed for that direction.** The active start relaxed toward inactive
   without it.
4. **No transducer is present.** A G-protein-free active state decaying (TM6 first,
   NPxxY lagging) may be correct physics. The Option A expectation may not hold without
   a transducer.

The §10 FAIL stands. The README diagnosis is corrected in place, with a dated
correction section. Tag v0.1-poc is left as it is so the history stays honest.

**Why not the Na⁺-biased run:**
- It addresses one of two hidden variables.
- It does not test the time explanation.
- It risks dragging an ion into the collapsed active pocket past DAMGO.

**Follow-ups launched instead.** They cost the same (two 4-GPU nodes for about
1.5 days) and give two interpretable answers:
1. **`opes_active_ext`** (job 204211). The prod2 active start continues from its
   checkpoints and STATE for +375 ns/walker. `opes.tpr` was extended with
   `gmx_mpi convert-tpr -extend 375000`, and the bias is identical. The rendered
   restart input was diffed against the original: it differs only in the RESTART and
   STATE_RFILE lines. COLVAR/KERNELS/STATE/cpt/tpr at 375 ns are backed up in
   `runs/prod2_backup_375ns/`. The v0.1 analysis is unchanged, since it is time ≤ 375 ns.
2. **`inactive_nona`** (job 204212, equilibration + 100 ns unbiased).
   - Built by `systems/make_nona.py` from the end of the inactive unbiased leg: the
     D2.50 Na⁺ (atom 36847, 0.34 nm from D116 CG) is swapped with a bulk water (atoms
     39580–39583). The ion ends up 4.55 nm from the bilayer midplane, ≥ 1.84 nm from any
     non-solvent atom.
   - Topology, index, composition and PLUMED indices are unchanged.
   - `plumed driver` on the swapped frame gives identical CVs and observables, except
     na_site 0.88 → 0.00.
   - Same equilibration protocol, seed 20261131.
   - The unbiased leg also shows whether Na⁺ re-enters spontaneously. If it does, the
     sodium-free condition is not stable and that is itself the answer.
   - OPES on it follows with the prod2 parameters unchanged.

Three questions:
1. **Does it match what was assumed?** The Day 4 diagnosis assumed a single hidden
   variable. The data show two, plus an incomplete relaxation.
2. **Is the next step still right?** The two runs separate time from sodium. Neither
   needs new CV definitions, so the comparison with prod2 is direct.
3. **What would failure look like?**
   - Extension: the active-start walkers stay where they ended. That points to a
     coordinate problem, not time.
   - inactive_nona: the start-dependence persists with sodium removed, which rules
     sodium out.
   - inactive_nona unbiased leg: Na⁺ returns to the site.

### Day 5 (2026-10-11, late): inactive_nona unbiased leg done, OPES launched

- **Unbiased leg (job 204212, 100 ns): Na⁺ did not return.** na_site > 0.1 in 0 of
  10,001 frames.
- **The receptor stayed inactive without sodium.** Second-half means:

  | | CV1 d_tm36 (nm) | CV2 rmsd_npxxy (nm) | d_hbond (nm) | na_site |
  |---|---|---|---|---|
  | inactive_nona | 0.638 | 0.065 | 0.519 | 0.00 |
  | inactive (with Na⁺) | 0.633 | 0.125 | 0.515 | 0.89 |

  - 96% of second-half frames lie in the inactive box. Maxima over the run: CV1 0.776,
    CV2 0.145.
  - Removing Na⁺ did not destabilise the inactive state on 100 ns. CV2 sits closer to
    the reference without it.
- **OPES on inactive_nona (job 209397, orcaga13).**
  - Command: `systems/opes.sh inactive_nona prod2 375 SIGMA=0.025,0.027 UW=1.50,0.47
    LW=0.53,0.0 BARRIER=50`. The parameters are prod2's, unchanged.
  - Output goes to `runs/prod2/inactive_nona`.
  - Walkers start from t = 62.5, 75, 87.5 and 100 ns of the unbiased leg, with seeds
    20261100+k (the same seeds as prod2).
  - Startup check: nker is identical across walkers (shared bias), and the walls are
    rendered as in prod2.
- **Active extension (job 204211) at about 422 ns per walker**, averaged over the 47 ns
  since the restart:
  - w0 and w2 are still partly deactivated. w2's CV2 continues to fall (0.24 → 0.20).
  - w1 and w3 stay active-like.
  - No committed crossings yet.

### Day 7 (2026-10-13): active-start extension finished (750 ns per walker)

- **Run.** Job 204211 ended normally at 750 ns per walker, 3.0 μs from the active
  start in total.
- **Analysis.** Job 213803 ran `final_analysis_job.sh` on `runs/v02_ext`.
  - `runs/v02_ext` is a directory of symlinks: `active` → `prod2/active`,
    `inactive` → `prod2/inactive`.
  - Output: `results/fes_v02_ext/`.
- **Committed crossings: still 0 in both directions.**
  - Closest approach to the inactive box over the whole run is 1.23 bandwidths, by
    walker w2 at 369 ns, i.e. *before* the extension.
  - Closest approach during the extension, per walker, in bandwidths:

    | walker | first 375 ns | extension |
    |---|---|---|
    | w0 | 5.4 | 3.8 |
    | w1 | 14.5 | 3.9 |
    | w2 | 1.2 | 1.4 |
    | w3 | 19.6 | 16.9 |

- **The relaxation seen on Day 5 did not continue.**
  - w0 and w2 had CV1 in the inactive range from about 30 to 330 ns. They then drifted
    back to CV1 0.8–1.0 nm.
  - w3 went the other way: CV2 has been inactive-like (0.05–0.12 nm) since about 450 ns,
    while CV1 stays at 1.1–1.4 nm.
  - At the end, w0 and w1 return to CV2 0.3–0.45 nm (active-like).
- **New observation: the two CVs move independently from the active start.**
  - Walkers reach the inactive value of each CV separately, but never both together.
  - The `paths_overlay` plot shows an L shape:
    - (CV1 ≈ 0.6, CV2 ≈ 0.25–0.35): TM6 closed, NPxxY still active-like;
    - (CV1 ≈ 1.2, CV2 ≈ 0.05–0.15): TM6 open, NPxxY inactive-like.
  - Neither corner connects to the inactive box in the sampling from this start.
  - The inactive start, by contrast, does not visit the (1.2, 0.1) corner.
- **§10 numbers**, prod2 at 375 ns → extension at 750 ns:

  | | prod2 (375 ns) | extension (750 ns) |
  |---|---|---|
  | Walker SD, active start, first → second half (kcal/mol) | 2.2 → 2.2 | 1.9 → 2.7 |
  | Grid points sampled from both starts | 706 | 1280 |
  | RMS difference on those points (kcal/mol) | 2.7 | 3.4 |
  | Second-half low-FE region: RMS (kcal/mol) | 2.2 | 1.7 |
  | Second-half low-FE region: max (kcal/mol) | 5.9 | 4.5 |
  | Second-half low-FE region: points | 338 | 328 |

  - The walker SD *grows* in the extension, so §10.4 still fails.
  - The second-half comparison sets the active start's 375–750 ns against the inactive
    start's 188–375 ns.
  - Basins sampled are still disjoint: active start {inactive 0, active 1}.
- **QC.**
  - PLUMED reweight cross-check: 0.004 / 0.007 kcal/mol RMS.
  - Membrane QC: no flags on any walker.
- **Interpretation.**
  - Doubling the simulated time did not complete the active→inactive transition in
    both CVs at once, and the walkers did not converge. Time alone is unlikely to be
    the limitation.
  - The decoupling of CV1 and CV2 is consistent with a missing coordinate, one that
    links TM6 and the NPxxY motif, but does not prove it.
  - The sodium-free comparison (job 209397) still decides whether Na⁺ contributes.

### Day 7 (2026-10-13, later): sodium-free inactive start finished; v0.2

- **Run.** Job 209397 ended normally at 375 ns per walker.
- **Analysis.** Job 213829 ran on `runs/v02_nona`.
  - The active start is its first 375 ns, from `prod2_backup_375ns`, so sodium is the
    only difference from v0.1. The active rows reproduce v0.1 exactly.
  - Output: `results/fes_v02_nona/`.
- **Results for the sodium-free start:**
  - 0 crossings toward active.
  - Closest approach to the active box: 4.7 kernel widths, at walker 2, 328 ns,
    (1.01, 0.24). With Na⁺ it was 2.5.
  - Na⁺ never re-entered the site.
  - Second-half difference vs the active start: 2.5 kcal/mol RMS, 6.1 max, over 410
    shared low-FE points. v0.1 had 2.2 RMS and 5.9 max over 338.
  - Walker SD 2.1 → 2.8 kcal/mol (grows). Q3 vs Q4 RMS 0.7.
  - rct at the end is 0.04 kJ/mol.
  - DAMGO contacts ≥ 305; d_salt ≤ 0.52 nm.
  - Reweight cross-check 0.006 kcal/mol; membrane QC has no flags.
- **paths_overlay.**
  - The sodium-free walkers populate a lobe at CV1 0.8–1.05 nm, CV2 0.05–0.25 nm, which
    the Na⁺-bound start barely visited. That is TM6 partly open with NPxxY inactive.
  - The lobe points at the active start's (1.2, 0.1) corner but does not join the
    active basin.
- **Conclusion.**
  - Na⁺ loosens TM6 in the inactive state but is not the rate-limiting difference.
  - Together with the extension, time and sodium are both ruled out at this scale.
  - TM6 and NPxxY move one at a time in every data set.
  - Two explanations remain: a missing coupling coordinate, or no stable active basin
    without a transducer.
  - README updated (section "Follow-ups (v0.2)"); tagged `v0.2-followups`.
- **Next options (not started):**
  - an Nb39-bound active start, to test the transducer explanation;
  - a learned or path CV trained on the 6 μs here, to test the coordinate explanation.

## v0.3 — agonist vs antagonist efficacy test (PLAN_v03.md)

### Day 8 (2026-10-13/14): G1–G3

- **G1 (commit e0541c1).** Four complexes on one set of 5C1M receptor coordinates;
  naloxone GAFF2/AM1-BCC (`NLX`, +1). Rebuilding both v0.1 starts with the generalised
  `build_topology.py` (job 214938) reproduces every file byte-for-byte bar the GROMACS
  header; every PLUMED index is unchanged.
- **G2 (commit e87a2f3).** Composition identical within each pair, all neutral:

  | pair | atoms | POPC | CHL | WAT | Na⁺ | Cl⁻ | box (Å) |
  |---|---|---|---|---|---|---|---|
  | `free_*` | 84.5k | 194 | 84 | 11,913 | 20 | 33 | 95 × 95 × 86 |
  | `nb39_*` | 162k | 194 | 84 | 30,844 | 73 | 84 | 95 × 95 × 150 |

  Two defects in new code paths, both fixed there:
  - packmol-memgen sized the z box from the extracellular side (±43 Å) and ignored Nb39
    reaching z = −54.3; packmol stopped with "Fixed molecule are outside the PBC box".
    `pack.sh` takes a z override (`--dims`), used for `nb39_*` only.
  - The Nb39 disulfide used a computed residue *index*; LEaP resolves `unit.N` by PDB
    residue *number*. Harmless for the receptor (numbered from 1), fatal for Nb39
    (numbered from 282). The cysteine numbers are now read from the residues.
- **G3 equilibration** (jobs 216031/216032/216713/216714, `long`, 1 × A30 each).
  - Seeds set explicitly, 20261401–04. `equilibrate.sh` derives its default from the
    *first letter* of the system name, which would give both members of each pair the
    same seed.
  - Unbiased leg 10 ns, as the plan says; the script's default is 100.
  - Measured single-A30 throughput at 4 fs: `free_*` 282 / 278 ns/day, `nb39_*`
    165 / 165 ns/day.

  Gate, over the 10 ns unbiased leg (active = CV1 ≥ 1.00 nm and CV2 ≥ 0.25 nm):

  | system | active | CV1 min (nm) | CV2 min (nm) | d_salt ≤ 0.6 nm | d_salt max | lig_cont 0–1 → 9–10 ns |
  |---|---|---|---|---|---|---|
  | v0.1 `active` (ref.) | 99.9 % | 1.133 | 0.250 | 100 % | 0.391 | 354 → 363 |
  | `free_dam` | 100 % | 1.147 | 0.359 | 100 % | 0.495 | 360 → 355 |
  | `free_nlx` | 100 % | 1.205 | 0.272 | 100 % | 0.399 | 283 → 286 |
  | `nb39_dam` | 100 % | 1.209 | 0.328 | 100 % | 0.420 | 345 → 363 |
  | `nb39_nlx` | 100 % | 1.192 | 0.369 | 100 % | 0.398 | 281 → 260 |

  Membrane QC (`analysis/membrane_qc.py`, v0.1 thresholds, stride 20; job 219526,
  `results/g3_qc/`): no flags on any system.

  | system | P–P (Å) | APL (Å²) | TM centre − midplane (Å) | TM tilt (°) | core waters | TM helicity |
  |---|---|---|---|---|---|---|
  | `free_dam` | 42.7 | 50.4 | −3.5 | 24.0 | 1.5 | 0.89 |
  | `free_nlx` | 42.1 | 51.1 | −3.0 | 24.0 | 0.2 | 0.89 |
  | `nb39_dam` | 41.5 | 52.4 | −1.9 | 23.5 | 0.5 | 0.89 |
  | `nb39_nlx` | 41.4 | 52.7 | −2.5 | 23.9 | 0.5 | 0.89 |

  Nb39–receptor interface (heavy-atom pairs within 4.5 Å, `analysis/nb_contacts.py`):
  `nb39_dam` 290 → 248 (retention 0.86, minimum 243), `nb39_nlx` 259 → 274 (1.06,
  minimum 222); 28–30 receptor residues in contact throughout. The interface holds in
  both.

  **G3 passes for all four systems.** Watch-items, not failures:
  `free_nlx` has the lowest CV2 minimum (0.272 nm, nearest the threshold), and
  `nb39_nlx` lig_cont fell 8 % with the D149 salt bridge intact throughout.

### Deviations from PLAN_v03.md

- **`lig_rmsd` is not logged by PLUMED.** The plan lists it as an added observable; the
  measurement templates carry `nb_cont` but not `lig_rmsd`. It will be computed from the
  trajectories (full system, every 100 ps) after scaffold alignment, as the plan defines
  it. The 10 ps observables are unaffected.
- **Nb39 production walltime.** At 165 ns/day a 500 ns replica needs about 73 h, past
  `replicas.sh`'s former fixed `-maxh 70`. `MAXH` is now a parameter; `nb39_*` run with
  `--time=5-00:00:00 MAXH=118` instead of finishing through a restart.
- **`nb_cont` during G3.** `equilibrate.sh` renders the plain template, so the G3 leg has
  no `nb_cont`; the interface was measured from the trajectory instead
  (`analysis/nb_contacts.py`, same 4.5 Å heavy-atom definition).
