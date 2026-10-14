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

**Why it failed (v0.2 follow-ups, [below](#follow-ups-v02-2026-10-13)).** Two
controlled follow-ups narrowed the explanation:
- **Not sampling time.** Doubling the active start to 750 ns per walker gave still no
  crossing, and convergence got worse.
- **Not the sodium ion.** A sodium-free inactive start was no closer to the active
  start's landscape.

What remains:
- **A missing coordinate.** TM6 and the NPxxY motif move independently in every run,
  and the bias on two separate CVs never makes them move together.
- **No transducer.** The receptor has no G protein or nanobody, so a stable active
  basin may not exist in this system.

These two are not yet distinguished. Details of the original verdict are
[below](#verdict).

## What the pilot establishes

The verdict above is about one landscape. The framework that produced it, and that the
12-receptor panel would reuse, is working and measured. Each item below is a file in this
repository, not an expectation.

**The toolchain runs.** GROMACS 2025.3 with PLUMED 2.10.1 (`OPES_METAD`, built from
source), multi-walker shared bias across 4 walkers, on a 150k-atom scale membrane system.
PLUMED's own reweighting agrees with the OPES estimate to 0.004-0.007 kcal/mol RMS
(`results/fes_prod2/plumed_reweight/comparison.json`).

**The acceptance criteria were written before the data, and they fired.** The four §10
criteria - reversible crossing, start-independence, no hysteresis, replica convergence -
rejected the v0.1 landscape. The supporting checks (ligand retention, bilayer integrity,
wall contact) passed at the same time, so the failure was localised to sampling rather
than to a broken system. A pipeline that reports a landscape only when it converges is
the point of the §10 table.

**The diagnosis was narrowed by controlled runs, not by argument.** Doubling the active
start to 750 ns/walker and removing the sodium ion each changed the answer not at all,
which excludes simulation length and Na+ coordination. TM6 and NPxxY reach their inactive
values separately in every dataset and never together - the signature of a missing
coupling coordinate, which is a protocol change rather than a cost problem.

**The systems are reproducible to the byte.** `systems/build_topology.py` rebuilt both
v0.1 endpoints after being generalised for a second ligand and a transducer chain: every
output file is identical bar the GROMACS header, and every PLUMED atom index is unchanged
(job 214938).

**Per-landscape cost is measured, not estimated.** On one A30 at 4 fs with HMR:
239 ns/day per walker for the 84.5k-atom system (4 walkers, one node), 282 ns/day for a
single replica, and 165 ns/day for the 162k-atom Nb39 system. These are the numbers that
size a panel campaign.

**The efficacy test is built and pre-registered.** `PLAN_v03.md` was committed before any
v0.3 data existed. Its four systems - DAMGO and naloxone, each with and without the Nb39
transducer, on one set of receptor coordinates with composition matched within each pair -
are built, equilibrated, and past every validity gate (`NOTES.md`, Day 8).


Quickest way to check that the setup works: [`examples/two-cv-smoke/`](examples/two-cv-smoke/)
(about 5 s on 4 CPU cores, no GPU and no allocation needed).

## Relation to the 12-receptor panel protocol

This repository is the μ-opioid prototype for a panel protocol: twelve peptide-activated
GPCRs of the enteroendocrine membrane, each solved in both the active and the inactive
conformation, each given an apo reference landscape and one landscape per bound peptide,
with efficacy read as the signed shift between them. The table states, element by element,
what the prototype has demonstrated and what it has not. It is the honest version of that
mapping; the gaps are the prototype's findings, not omissions.

| Panel protocol element | Prototype status |
|---|---|
| 150k-atom receptor + bilayer system, GROMACS + PLUMED | **Demonstrated.** 84.5k atoms without transducer, 162k with Nb39; build is byte-reproducible |
| Multi-walker OPES along activation coordinates | **Demonstrated.** 4 walkers, shared bias, ~6 μs total; reweighting agrees to 0.004–0.007 kcal/mol |
| CV1 TM6 displacement, CV2 NPxxY | **Demonstrated as biased coordinates.** Both reach their inactive values |
| CV3 connector region | **Not yet biased.** Logged as `chi_w295` only. The prototype's central finding — TM6 and NPxxY never move together — is what motivates promoting it to a biased coordinate |
| Reversible driving between experimental endpoints | **Not achieved on 6 μs with two CVs.** 0 committed crossings in either direction |
| ΔG between endpoints | **Not obtained.** The two starts' basins are disjoint, so no endpoint difference can be read |
| Apo reference landscape per receptor | **Not run.** Every prototype run is DAMGO-bound |
| Peptide pose carried from the co-folding stage | **Not used.** Ligand poses come from experimental structures (5C1M, 10TM, 9PXU) |
| Efficacy as signed liganded − unliganded shift | **Not computed.** Requires the two landscapes above |
| Convergence from both starts and across replicas | **Implemented and enforced.** The four §10 criteria rejected the v0.1 landscape; the supporting physical checks passed at the same time, localising the failure to sampling |
| Recovering documented agonist/antagonist efficacy | **In progress, pre-registered.** v0.3: DAMGO vs naloxone, each ± Nb39, built and past every validity gate before production |
| Regression model over accumulated landscapes | **Not started.** Needs converged landscapes as input |
| Per-landscape cost for campaign sizing | **Measured.** 239 ns/day per walker (84.5k atoms, A30, 4 fs HMR); 165 ns/day at 162k atoms |

Two consequences worth stating plainly, because they shape what a panel campaign should
buy:

- **The limitation the prototype found is a coordinate problem, not a hardware problem.**
  Doubling simulation time made convergence worse, not better. Adding GPUs to a
  two-coordinate OPES run on this system would buy more unconverged sampling. The panel
  protocol therefore needs the third coordinate and the transducer-stabilised reference,
  both of which this repository is now set up to test.
- **The prototype's throughput requirement is parallel, not coupled.** OPES walkers
  exchange only a small shared-bias file, and the v0.3 replicas are independent by
  construction — coupling them through GROMACS's `-multidir` cost 91% of wall time in
  inter-simulation signalling and was removed. Activation MD fills a large machine by
  running many landscapes at once, not by tightly coupling one.

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
  (TM6 closed, NPxxY displaced). In the inactive-start map that region lies about
  8 kcal/mol or more above the minimum.
- **Na⁺ at D2.50 was a slow, unbiased difference between the starts.** It stayed bound
  in every inactive-start walker (≤ 0.1% of frames unbound) and never bound in the
  active start.
- **The bias reached its BARRIER range without driving a coupled transition.** opes.bias
  spans −50 to about +9 kJ/mol, and rct stays flat near −1 kJ/mol.

**Diagnosis**, using the runbook §10 table rows "no crossings", "landscapes disagree"
and "CV1 reached but CV2 flat":
- **There is a slow variable outside the 2D space.** At the same CV values (CV1
  0.55–0.72 nm, CV2 0.18–0.32 nm), frames from the two starts still differ in two
  recorded observables:
  - Na⁺ at D2.50: 0.00 in the active start vs 0.91 in the inactive start;
  - the R167–T281 side-chain contact: 0.78 nm (broken) vs 0.47 nm (formed).

  W295 χ1 and the ligand observables do not differ.
- **Which variable is rate-limiting is not established.** See the correction below.

### Correction (2026-10-11): diagnosis not established

The first version of this section said the two coordinates are insufficient and named
Na⁺ as the prime suspect. On re-analysis that goes beyond the data:
- **The active-start relaxation had not finished.** The two walkers that closed TM6
  were still drifting toward the inactive state when the run stopped. Their means over
  0–50, 100–200 and 250–375 ns were:
  - walker 0: CV2 0.40 → 0.36 → 0.31 nm; side-chain contact 0.92 → 0.77 → 0.66 nm;
  - walker 2: CV2 0.31 → 0.28 → 0.24 nm; side-chain contact 1.05 → 0.92 → 0.84 nm.

  So too little sampling time explains the active-start side as well as a missing
  coordinate does.
- **Na⁺ is not shown to be causal.** The active start was relaxing toward inactive
  without any Na⁺ at the site. Na⁺ may matter mainly in the activating direction, where
  the inactive start kept its ion throughout.
- **No transducer.** Neither system has a G protein or nanobody. A G-protein-free,
  agonist-bound receptor is not generally expected to hold a fully active conformation
  stably. The active-start walkers decaying toward inactive (TM6 first, NPxxY lagging)
  may therefore be correct physics, and a start-independent landscape for this system
  may have the inactive state lowest. The runbook's Option A expectation that "agonist
  favours active" may not hold without a transducer.

The §10 verdict (FAIL) is unchanged; only the explanation is withdrawn. Follow-up runs
designed to separate these explanations (NOTES.md, Day 5):
1. **The active start, extended** from its checkpoints by 375 ns per walker with the
   same bias. Does deactivation complete in both CVs?
2. **A sodium-free inactive start** (`systems/inactive_nona`): the D2.50 Na⁺ swapped
   with a bulk water, so both starts are chemically identical. Is the start-dependence
   caused by Na⁺?

### Follow-ups (v0.2, 2026-10-13)

Both follow-ups used the prod2 parameters unchanged and went through the same analysis
(`analysis/final_analysis_job.sh`). The PLUMED reweight cross-check agrees to
≤ 0.007 kcal/mol in both, and the membrane QC raised no flags on any walker.

**1. Active start extended to 750 ns per walker** (`results/fes_v02_ext/`, 3.0 μs from
the active start).
- **Still 0 / 0 committed crossings.** The closest approach to the inactive box
  (1.23 kernel widths) happened at 369 ns, before the extension. During the extension,
  no walker got closer than 1.4 kernel widths.
- **The Day-5 relaxation did not continue.**
  - Walkers 0 and 2 drifted back from CV1 ≈ 0.6 nm to 0.8–1.0 nm.
  - Walker 3 moved its NPxxY to inactive values (0.05–0.12 nm) at about 450 ns,
    but kept TM6 open (CV1 1.1–1.4 nm).
- **The sampled region is L-shaped.** Walkers reach the inactive value of one CV or the
  other, never both:
  - (CV1 ≈ 0.6, CV2 ≈ 0.3): TM6 closed, NPxxY active-like;
  - (CV1 ≈ 1.2, CV2 ≈ 0.1): TM6 open, NPxxY inactive-like.
- **Convergence got worse.** Walker-to-walker SD rose from 1.9 to 2.7 kcal/mol (first
  vs second half).
- **The two starts still disagree.** The second-half difference against the inactive
  start, over their shared low-FE region, is 1.7 kcal/mol RMS (4.5 max). Their basins
  are still disjoint.

**2. Sodium-free inactive start** (`systems/inactive_nona`, `runs/prod2/inactive_nona`,
`results/fes_v02_nona/`).
- **Setup.** The D2.50 Na⁺ was swapped with a bulk water.
  - 100 ns unbiased: Na⁺ never returned, and the inactive state stayed stable (96% of
    frames in the box).
  - Then 4 × 375 ns OPES, compared against the active start's first 375 ns. This is
    the same comparison as v0.1, with sodium the only change.

  | | inactive start with Na⁺ (v0.1) | sodium-free |
  |---|---|---|
  | Crossings toward active | 0 | 0 |
  | Closest approach to the active box (kernel widths) | 2.5 | 4.7 |
  | Second-half Δ vs active start, RMS / max (kcal/mol) | 2.2 / 5.9 (338 pts) | 2.5 / 6.1 (410 pts) |
  | Walker-to-walker SD, first → second half (kcal/mol) | 0.75 → 1.9 | 2.1 → 2.8 |

- **Removing Na⁺ loosened the inactive state but did not open a path to the active
  state.**
  - Sodium-free walkers populate a new lobe at CV1 0.8–1.05 nm, CV2 0.05–0.25 nm (TM6
    partly open, NPxxY inactive).
  - That lobe points toward the active start's TM6-open/NPxxY-inactive corner, but
    never joins the active basin.
  - Start-independence did not improve.
- **So Na⁺ is not the rate-limiting difference** on this timescale.

**Interpretation.**
- In all three data sets, TM6 and NPxxY each reach both of their values, but only one
  at a time.
- Biasing them as two independent CVs does not supply the coupling between them.
- Two explanations remain, and these runs cannot tell them apart:
  - **a missing coordinate** that couples TM6 and NPxxY, e.g. the R3.50 / Y7.53
    microswitch region or a path CV between the endpoints;
  - **no stable active basin** without a transducer.

  A G-protein-mimetic nanobody (Nb39) bound to the active start would test the second:
  if the active start then holds, the coordinate question can be asked cleanly.

The 6 μs of OPES and 300 ns unbiased here is a ready-made training set for a learned CV,
such as Deep-TDA or Deep-LDA on the L-shaped region. It holds 12 biased walker trajectories with all
secondary observables recorded.

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
| `results/` | endpoint values; daily snapshots; `fes_prod2/` (v0.1), `fes_v02_ext/` and `fes_v02_nona/` (follow-ups): FES grids, summaries, figures, strided COLVARs, QC |
| `examples/two-cv-smoke/` | self-contained runnable check of both CVs |
| `NOTES.md` | stage-by-stage plausibility log and every deviation |

Reproduce the analysis from a finished run:
`sbatch analysis/final_analysis_job.sh runs/prod2 results/fes_prod2`.
The follow-ups use directories of symlinks that pair the starts:
- `runs/v02_ext` links `{active,inactive}` → `prod2/{active,inactive}`.
- `runs/v02_nona` links `active/w*/COLVAR.*` → `prod2_backup_375ns`, i.e. the first
  375 ns, and `inactive` → `prod2/inactive_nona`.

## Where the trajectories live

Not in git; they total about 20 GB of OPES xtc plus checkpoints and OPES KERNELS/STATE. They are on
ORCA (Portland State University):

```
/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation/
  systems/{active,inactive,inactive_nona}/prod.*   100 ns unbiased endpoints
  runs/prod2/{active,inactive}/w{0..3}/            OPES production (opes.xtc, COLVAR.k, KERNELS, STATE);
                                                   active extended to 750 ns
  runs/prod2/inactive_nona/w{0..3}/                sodium-free inactive start, 375 ns
  runs/prod2_backup_375ns/active/                  active-start COLVAR/KERNELS/STATE at 375 ns
  runs/prod1_softwalls/                            superseded first launch (soft walls), not used
```

ORCA `/scratch` is not backed up and is subject to purge. Copy these elsewhere before
relying on them. Strided COLVARs (every 100 ps) for every walker are committed in
`results/fes_*/colvar_strided/`.
