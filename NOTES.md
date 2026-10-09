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
The checkpoint timestamps in  show the stall.

There are no errors, LINCS warnings or energy anomalies, and constraint RMSD is
unchanged. The stall hit two separate nodes simultaneously, so the shared resource is
the likely cause:  is an NFS mount at 94% full. PLUMED appends COLVAR and
KERNELS every 500 steps on every walker, which is sensitive to NFS latency.

Effect: lost wall-clock only, not data. Position at 06:45 is 204 ns (active) and 211 ns
(inactive) per walker, about 17 h from completion if there are no further stalls. The
job time limit (mdrun -maxh 70) is not at risk.

If it recurs: restart from checkpoint with the run directory on node-local disk,
copying back at checkpoints.
