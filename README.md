# oprm1-opes-activation

OPES enhanced-sampling molecular dynamics of the **human μ-opioid receptor** (OPRM1,
UniProt P35372) with the agonist **DAMGO** bound, in a POPC/cholesterol bilayer. The
goal is the conformational free-energy surface (FES) of receptor activation along two
hand-picked coordinates:
- the TM3–TM6 Cα distance;
- the NPxxYA backbone RMSD.

The run is designed so that the FES computed from an **active** start and from an
**inactive** start can be compared. If they agree, the landscape is real. If they
disagree, the two coordinates are missing a slow motion. This is a week-one
proof of concept for a methods proposal, not a biology result.

**Verdict (v0.1-poc): FAIL on all four §10 criteria. This is a clean negative.**
3 μs of multi-walker OPES (2 starts × 4 walkers × 375 ns) produced:
- no reversible activation transition;
- two start-dependent landscapes.

The two coordinates are not sufficient for μOR activation at this sampling budget.
Details and the diagnosis are [below](#verdict).

Quickest way to check that the setup works: [`examples/two-cv-smoke/`](examples/two-cv-smoke/)
(about 5 s on 4 CPU cores, no GPU and no allocation needed).

## Collective variables (human numbering)

| CV | Definition | Atoms |
|---|---|---|
| **CV1 `d_tm36`** | Cα–Cα distance **Arg167 (3.50) – Thr281 (6.34)** | `@CA-167`, `@CA-281` via `MOLINFO` |
| **CV2 `rmsd_npxxy`** | Cα RMSD of **Asn334–Ala339** (NPVLYA, 7.49–7.54) against the equilibrated inactive structure, after superposition on the TM1–TM5 scaffold (134 Cα) | dual-weight reference `cv/ref_inactive_npxxy.pdb`: occupancy = alignment, B-factor = displacement |

> ⚠️ **Numbering.** The published Filizola-group papers use **rodent (mouse/rat)
> numbering, which is human − 2**. Rodent R165 / T279 / N332 are human R167 / T281 /
> N334, so add 2 to a published residue number to get the human one. Every number in
> this repository is human P35372. This is asserted on every deposited structure
> (`cv/assert_numbering.py`) and again at topology build time
> (`systems/build_topology.py`).

Secondary observables are recorded but never biased:
- R167–T281 side-chain H-bond;
- the W295 χ1 toggle;
- the DAMGO–D149 salt bridge;
- DAMGO–receptor contacts;
- Na⁺ occupancy at D2.50 (Asp116).

Full definitions, measured endpoint values and atom indices are in
[`cv/cv_definition.md`](cv/cv_definition.md).

## Structures

Option A of the runbook: the agonist DAMGO bound, in both directions. Full record:
[`structures/structures.md`](structures/structures.md).

| Start | Source | Removed | Reverted |
|---|---|---|---|
| active | **10TM** chain R + DAMGO (chain S), 3.0 Å | second receptor copy (F) and its DAMGO (H); 2 resolved cholesterols; residues 66–67 | **W158F** (thermostabilising F158W) |
| inactive | **9PXU** chain R, 3.4 Å; DAMGO transplanted from 10TM by scaffold superposition | nanobody (C, K), Fab (H, L), naloxone; residues 348–354 (BRIL is in DBREF but has no modelled coordinates) | **L266M**, **R271K** (SEQADV "conflict") |

Common to both:
- receptor residues 69–346, capped ACE/NME;
- identical protonation (propka 3.5.1, unified across starts: Asp116 deprotonated;
  His173, His225, His299 and His321 as HID);
- the Cys142–Cys219 disulfide;
- no loops modelled.

The 9PXU Na⁺ at D2.50 is kept. 10TM has none, so in the active start that ion starts in
bulk solvent.

## System and method

- **Force fields:** Amber ff19SB + Lipid21 + OPC water, GAFF2/AM1-BCC for DAMGO (net +1).
- **Bilayer:** POPC:CHL1 7:3 built with packmol-memgen, oriented with OPM's 8F7Q frame.
  194 POPC, 84 CHL, 150 mM NaCl; 84,511 atoms in a 95 × 95 × 86 Å box. Both starts
  have identical composition and atom order.
- **Simulation:** GROMACS 2025 (PLUMED-patched 2025.0 for mdrun, 2025.3 for grompp),
  HMR with a 4 fs step, 310 K, 1 bar semi-isotropic.
- **Endpoints:** each start equilibrated and then run 100 ns unbiased. The endpoint CV
  values are measured from these runs, not looked up.
- **Bias:** PLUMED 2.10.1 `OPES_METAD`, built from source with OPES enabled.
  - BARRIER 50 kJ/mol, PACE 500, SIGMA 0.025 / 0.027 nm (from the endpoint spread).
  - 4 walkers per start sharing one bias (`WALKERS_MPI`).
  - Upper/lower walls (KAPPA 50,000 kJ mol⁻¹ nm⁻²) at the measured endpoints plus margin.
- **Reweighting:** w = exp(V_total/kT), with V the total bias (OPES plus walls).
  Implemented in `analysis/fes.py`. It agrees with PLUMED `REWEIGHT_BIAS` + `HISTOGRAM`
  to 0.004–0.007 kcal/mol RMS, with identical minima
  (`results/fes_prod2/plumed_reweight/comparison.json`).

Versions and build recipes are in [`env/`](env/versions.md). Every decision and
deviation from the runbook, stage by stage, is in [`NOTES.md`](NOTES.md).

**Relation to the precedent.** These are the same two order parameters used by the
Filizola group (Kapoor, Provasi & Filizola, *Mol Pharmacol* 2020;98:475,
doi:10.1124/mol.119.119339; Kapoor et al., *Sci Rep* 2017;7:11255,
doi:10.1038/s41598-017-11483-8). Those landscapes came from **adaptive sampling plus
Markov state models**. This pilot uses **OPES biased sampling** with reweighting: the
same coordinates, but a different method. It is not a reproduction of their method.
- Meral, Provasi & Filizola (*J Chem Phys* 2018;149:224101, doi:10.1063/1.5060960)
  used metadynamics, but **biased path CVs (S, Z) built from contact maps**. These two
  descriptors were only the projection space, so that paper is precedent for the
  coordinates, not for biasing them directly.
- The review by Marino, Shang & Filizola (*Br J Pharmacol* 2018;175:2834,
  doi:10.1111/bph.13774) is cited for context only. Its structural claims, such as the
  size of the TM6 displacement, could not be verified and are not relied on here.

## Verdict

Production run `runs/prod2`: 2 starts × 4 walkers × 375 ns = **3.0 μs biased**.

| §10 criterion | Result | |
|---|---|---|
| 1. Reversible crossing (≥ 3–5 each way, both CVs) | **0 / 0** in both starts | FAIL |
| 2. Start-independence (≤ 1–2 kcal/mol) | Basins are disjoint: neither start's FES has any sampled point in the other start's basin. Where both sampled below 8 kcal/mol (second half, 338 grid points), the two maps differ by 2.2 kcal/mol RMS (5.9 kcal/mol max) | FAIL |
| 3. No hysteresis | Forward and reverse paths occupy different regions; they overlap only near CV1 ≈ 0.6 nm, CV2 ≈ 0.2–0.35 nm | FAIL |
| 4. Convergence | Walker-to-walker SD is 2.2 → 2.2 kcal/mol (active start) and 0.75 → 1.9 kcal/mol (inactive start; it grows). Q3-vs-Q4 RMS is 2.1 / 0.9 kcal/mol | FAIL |

Supporting checks, all of which pass:
- DAMGO stayed bound throughout (DAMGO–receptor heavy-atom coordination ≥ 283, against 300–550 typical; DAMGO N–D149 CG ≤ 0.50 nm).
- The fold and bilayer stayed intact on all 8 biased walkers (`results/fes_prod2/membrane_qc_*`):
  TM helicity ≥ 0.84 in every frame, bundle tilt ≤ 26°, P–P thickness 43.8–44.7 Å,
  APL ≈ 48 Å², and ≤ 5 waters in the hydrophobic core. Lipid heavy atoms appeared inside
  the bundle in only 1 of 608 frames, as a single atom.
- Walls were touched in ≤ 2.4% of frames in 7 of 8 walkers. The exception is inactive-start walker 1 at 7.2%, on the lower CV1 wall at 0.53 nm: that walker prefers 0.50–0.51 nm, so the FES below about 0.55 nm on CV1 is wall-limited.

**What happened** (`results/fes_prod2/figures/`):
- **CV1 alone moves freely in both directions.**
  - Two of the four active-start walkers closed TM3–TM6 to inactive distances
    (0.55–0.65 nm) within 100 ns, but NPxxY stayed at 0.2–0.35 nm. The inactive value is
    ≤ 0.12 nm.
  - The other two active-start walkers did the reverse: NPxxY dropped to about 0.1 nm
    while TM6 stayed open.
  - Inactive-start walkers reached TM3–TM6 up to 1.22 nm, and NPxxY up to 0.48 nm, but
    never both at once. The closest approach to the active box was (1.11, 0.24) nm,
    2.5 kernel widths outside it.
- **The two CVs never changed together in one trajectory.** The closest the active start
  came to the inactive box was (0.636, 0.130) nm, 1.2 kernel widths outside it.
- **The active start's FES has a secondary basin at CV1 0.6–0.7 nm, CV2 0.3–0.4 nm**
  (TM6 closed, NPxxY displaced). The inactive start never visits it.
- **Na⁺ at D2.50 was a slow, unbiased difference between the starts.** It stayed bound
  in every inactive-start walker (≤ 0.1% of frames unbound) and never bound in the
  active start.
- **The bias reached its BARRIER range without driving a coupled transition.** opes.bias
  spans −50 to about +9 kJ/mol, and rct stays flat near −1 kJ/mol.

**Diagnosis**, using the runbook §10 table rows "no crossings", "landscapes disagree"
and "CV1 reached but CV2 flat":
- Rate-limiting degrees of freedom lie outside the 2D space. The prime suspects are:
  - Na⁺ occupancy (and the D2.50 protonation that a fixed-charge model cannot change);
  - the TM7 / Y7.53 side-chain rearrangement that NPxxY backbone RMSD does not resolve.
- This is the outcome the runbook calls a passing week as a *clean negative*: two
  hand-picked coordinates are not enough, which is the argument for learning the
  coordinate.
- The 3 μs here is a ready-made training set for that. It holds 8 biased walkers with
  all secondary observables recorded.
- A cheaper next step: repeat with Na⁺ occupancy as a third biased CV, or with BARRIER
  raised to 80–100 kJ/mol.

**External validation is out of scope.** Everything above is an internal check, and
internal convergence detects variance, not force-field bias. Agreement with measured
function (cAMP, BRET, GTPγS, Kᵢ) was not tested this week.

## Repository layout

| Path | Contents |
|---|---|
| `env/` | build scripts for PLUMED and patched GROMACS, environment files, versions |
| `structures/` | fetch, numbering assertion, cleanup; cleaned complexes; DAMGO parameters |
| `cv/` | CV definitions, scaffold, RMSD references, PLUMED inputs (`cv/plumed/`) |
| `systems/` | membrane build, topology, equilibration and OPES launchers, mdp files, built systems |
| `analysis/` | endpoint measurement, membrane QC, daily OPES monitor, FES/reweighting, convergence |
| `results/` | endpoint values; daily snapshots; `fes_prod2/` (FES grids, summaries, figures, strided COLVARs, QC) |
| `examples/two-cv-smoke/` | self-contained runnable check of both CVs |
| `NOTES.md` | stage-by-stage plausibility log and every deviation |

Reproduce the analysis from a finished run:
`sbatch analysis/final_analysis_job.sh runs/prod2 results/fes_prod2`.

## Where the trajectories live

Not in git; they total about 9 GB of OPES xtc plus checkpoints and OPES KERNELS/STATE. They are on
ORCA (Portland State University):

```
/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation/
  systems/{active,inactive}/prod.*          100 ns unbiased endpoints
  runs/prod2/{active,inactive}/w{0..3}/      OPES production (opes.xtc, COLVAR.k, KERNELS, STATE)
  runs/prod1_softwalls/                      superseded first launch (soft walls), not used
```

ORCA `/scratch` is not backed up and is subject to purge. Copy these elsewhere before
relying on them. Strided COLVARs (every 100 ps) for all 8 walkers are committed in
`results/fes_prod2/colvar_strided/`.
