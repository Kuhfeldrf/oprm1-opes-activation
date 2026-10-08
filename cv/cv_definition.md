# Collective variables

All residue numbers are **human P35372** (GPCRdb `oprm_human`). The published Filizola-group
papers use **rodent numbering, which is human −2**; for example, rodent R165/T279/N332
are human R167/T281/N334. Human residue 165 is Val.

| Role | Ballesteros–Weinstein | Human residue |
|---|---|---|
| D(RY) | 3.49 | Asp166 |
| **R**(DRY) | **3.50** | **Arg167** |
| Y(DRY) | 3.51 | Tyr168 |
| L (no ionic-lock Glu/Asp in μOR) | 6.30 | Leu277 |
| **T** | **6.34** | **Thr281** |
| C W x P | 6.47–6.50 | Cys294-Trp295-Thr296-Pro297 |
| **NPxxYA** | **7.49–7.54** | **Asn334 Pro335 Val336 Leu337 Tyr338 Ala339** (NPVLYA) |
| D2.50 (Na⁺ site) | 2.50 | Asp116 |
| D3.32 (ligand amine) | 3.32 | Asp149 |

Verified three ways: GPCRdb residue service; `cv/assert_numbering.py` on all 8 deposited
entries; and the build-time assertion in `systems/build_topology.py` on the simulated
topology.

## Biased CVs

**CV1 `d_tm36`** — Cα(Arg167)–Cα(Thr281) distance, in nm.

**CV2 `rmsd_npxxy`** — RMSD of the six NPxxYA Cα atoms against the **inactive** reference,
in nm. It uses a dual-weight reference (`cv/build_rmsd_ref.py`, runbook §7.2):

| Atom set | Residues | Count | occupancy (align) | B-factor (displace) |
|---|---|---|---|---|
| TM1–TM5 scaffold Cα | 72–96, 105–131, 140–170, 184–206, 228–255 | 134 | 1.00 | 0.00 |
| NPxxYA Cα | 334–339 | 6 | 0.00 | 1.00 |

The scaffold is defined in `cv/scaffold.py` from GPCRdb segment boundaries (TM1 67–98,
TM2 103–133, TM3 138–173, TM4 182–208, TM5 226–264). It excludes:
- helix ends;
- TM5's cytoplasmic end (256–264), which moves with TM6 on activation;
- ICL2, ICL3 and ECL2 (the 9MQx ECL2 gap is at 225–226);
- all of TM6 and TM7.

On deposited structures the scaffold superposes to 0.9–1.4 Å across active and inactive
states, so it is a stable frame (NOTES.md).

Reference coordinates: initially the built 9PXU-derived inactive system (`system_ref.pdb`
of `systems/inactive`). They are replaced by the **equilibrated** inactive structure before
any OPES run (§6.2); the provenance line in the PDB REMARK records which one is in use.
PDB serials in the reference equal topology indices, which are identical in both
systems by construction.

## Atom indices and chain IDs

- **Deposited files (author chain IDs, used only by `cv/*.py` on raw mmCIF):** 10TM `R`
  (copy `F`); 8F7Q `R`, `M`; 8Y72, 9PPQ, 9PXU `R`; 9MQH, 9MQI, 9MQJ `A`.
- **Simulated systems:** single receptor chain `R`, ligand chain `L`. Topology order is
  receptor (ACE68…NME347), then DAMGO, then lipids, water and ions. The two starts have
  identical composition and therefore **identical atom indices**.
- PLUMED selects by residue through `MOLINFO STRUCTURE=system_ref.pdb`. That file is
  written by `systems/build_topology.py` with canonical numbering restored. tleap and
  packmol-memgen both renumber the receptor from 1, so canonical R167 would otherwise
  be residue 100. Explicit 1-based indices (DAMGO, Na⁺, receptor heavy atoms) are
  substituted by `systems/render_plumed.py` from `build_report.json`, never typed by
  hand. Resolved indices for the production systems are listed below once the systems
  are built.

## Secondary observables (recorded, never biased)

| Name | Definition | Why |
|---|---|---|
| `d_hbond` | Arg167 CZ – Thr281 OG1 | the real inactive-state TM3–TM6 restraint in μOR (no ionic lock: 6.30 is Leu277) |
| `chi_w295` | Trp295 χ1 | CWxP toggle switch |
| `d_salt` | DAMGO N1 (amine) – Asp149 CG | ligand anchored by the conserved D3.32 salt bridge |
| `lig_cont` | DAMGO heavy atoms × receptor heavy atoms, switching r₀ 0.45 nm | ligand retention, not dragged out by the bias |
| `na_site` | Asp116 CG × all Na⁺, r₀ 0.40 nm | sodium-site occupancy: a slow DOF outside the CV space |

## Endpoint values

### Deposited coordinates (Å; `cv/deposited_cv_values.txt`)

| | CV1 Cα–Cα | CV2 vs 9PXU |
|---|---|---|
| inactive (9MQH, 9MQJ, 9MQI, 9PXU) | 6.58–6.75 | 0.00–1.15 |
| active (10TM, 8F7Q×2, 8Y72) | 12.24–13.18 | 3.63–4.10 |
| 9PPQ (GPCRdb: active) | 13.00 | 4.06 |

### Equilibrated systems (unbiased 100 ns each, first 10 ns discarded; `results/endpoints.json`)

CV2 is measured against the **equilibrated** inactive reference: the scaffold-aligned mean
over 50–100 ns of the inactive leg (`cv/rebuild_ref_equilibrated.py`; crystal-derived
original kept as `cv/ref_inactive_npxxy.crystal.pdb`).

| | CV1 d_tm36 (nm) | CV2 rmsd_npxxy (nm) | R167 CZ–T281 OG1 (nm) | Na⁺ site | DAMGO–D149 (nm) |
|---|---|---|---|---|---|
| active leg | 1.199 ± 0.055 | 0.320 ± 0.027 | 1.16 ± 0.17 | 0.00 | 0.321 ± 0.015 |
| inactive leg | 0.633 ± 0.025 | 0.060 ± 0.029 | 0.52 ± 0.05 | 0.92 ± 0.06 | 0.319 ± 0.017 |
| separation | 0.565 (10 SD) | 0.260 (9 SD) | | | |

The active state held without its G protein for 100 ns: CV1 drifted −1.3 Å between
the first and last 10 ns, against a 3 Å failure threshold. The inactive state kept the
R167–T281 contact and its Na⁺.

**OPES parameters derived from these** (`opes_proposal`):
- SIGMA 0.025, 0.027 nm: the smaller within-well SD per CV.
- Lower walls 0.53 nm (CV1) and 0 (CV2, no-op).
- Upper walls 1.50 nm (CV1) and 0.47 nm (CV2).
- BARRIER 50 kJ/mol (≈ 12 kcal/mol).

Basin boxes for crossing counts (nm; d_tm36 range, rmsd_npxxy range):
- inactive: [0.558, 0.684] × [0, 0.118]
- active: [1.089, 1.364] × [0.265, 0.401]
