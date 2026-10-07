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
