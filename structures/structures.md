# Structures: selection, cleanup, every deletion and reversion

Produced by `structures/fetch.sh` → `cv/assert_numbering.py` → `structures/clean.py`.
Machine-readable record per start: `structures/clean/<start>/clean_report.json`.
The cleaned complexes are committed as `structures/clean/<start>/complex_opm.pdb`.

## Final choices (Option A: DAMGO-bound μOR, two starts)

| Start | Source | Res. | Receptor kept | Why |
|---|---|---|---|---|
| **active** | **10TM** chain R + DAMGO chain S | 3.0 Å | 69–346 (+ACE68, NME347) | No G protein or nanobody in the model. One reversion. DAMGO is experimental. |
| **inactive** | **9PXU** chain R + DAMGO transplanted from 10TM | 3.4 Å | 69–346 (+ACE68, NME347) | Fully resolved 68–354 with no ECL2 gap, the Na⁺ site resolved, and an orthosteric pocket already occupied. Higher resolution than 9MQH (3.9 Å). |

Both starts have **identical sequence, caps, protonation and disulfide**. They differ only
in conformation and in where one Na⁺ starts (see below).

## Deletions

| Start | Deleted | Detail |
|---|---|---|
| active (10TM) | auth chain **F** | second receptor copy (identical CVs to R; one copy kept) |
| | auth chain **H** | DAMGO bound to copy F |
| | CLR401, CLR402 (chain R) | resolved cholesterols; the membrane adds 30 mol% CHL uniformly to both starts |
| | residues 66–67 | outside the common range |
| inactive (9PXU) | auth chains **C, K** | nanobody chains |
| | auth chains **H, L** | Fab heavy/light |
| | naloxone (A1APV501) | replaced by DAMGO |
| | residues 348–354 | outside the common range (10TM ends at 347) |
| | BRIL (P0ABE7, auth −108…−8) | in DBREF only; **no BRIL coordinates are modelled** in the deposited file, so nothing to delete or cap |
| both | all waters, all hydrogens, alternative conformers | — |

## Reverted mutations (checked against SEQADV, not trusted from the runbook)

| Start | Engineered | Reverted to | SEQADV label |
|---|---|---|---|
| active (10TM) | Trp158 | **Phe158** | "engineered mutation" |
| inactive (9PXU) | Leu266 | **Met266** | "conflict" |
| inactive (9PXU) | Arg271 | **Lys271** | "conflict" |

pdbfixer `applyMutations` + `addMissingAtoms`. No structure needed more than two
reversions, so §4.5's rigidity flag is not triggered. **8F7Q also carries F158W**
(SEQADV "conflict"); it is not a start, but it would need reverting if it were.

## Completed side chains (pdbfixer, heavy atoms only)

- active: 10 residues, including the reverted F158.
- inactive: 34 residues, mostly surface Arg/Lys/Glu side chains truncated in the 3.4 Å
  map; plus the reverted M266/K271. The full list is in `clean_report.json`.

## Loops

No loops were modelled. 69–346 is continuous in both sources. ICL3 (265–269 GPCRdb)
is resolved in both. The 9MQH/9MQI/9MQJ ECL2 gap (225–226) is avoided by not using them
as starts.

## Caps

Residue 68 becomes **ACE** (its CA → CH3, C, O kept) and residue 347 becomes **NME** (N,
CA → CH3 kept), from the experimental backbone atoms.

## Protonation (pH 7.4, propka 3.5.1 via pdb2pqr 3.6.1), unified across starts

propka run separately on each start **disagreed** at two residues:

| Residue | active (10TM) | inactive (9PXU) | Used for both |
|---|---|---|---|
| Asp116 (D2.50) | ASH (protonated) | ASP | **ASP**: the Na⁺ it coordinates is kept |
| His173 | HIE | HID | **HID** (inactive assignment) |

Both starts are the same chemical system, so one pattern is applied to both. A
protonation difference between starts would appear as a fake free-energy difference.
Final non-default states: **HID173, HID225, HID299 (H6.52), HID321**. All other
titratable residues are in their standard states.

> The active-start propka result (D2.50 protonated when Na⁺ is absent) matches the
> hypothesis that D2.50 protonation accompanies activation. A fixed-charge simulation
> cannot capture this, and it is recorded as a model limitation rather than a choice.

## Disulfide

**Cys142–Cys219** (ECL2–TM3), detected geometrically (SG–SG < 2.5 Å) in both starts,
named CYX and bonded explicitly in tleap.

## DAMGO

- Chemistry: Tyr-D-Ala-Gly-N-MePhe-Gly-ol (CCD TYR-DAL-GLY-MEA-ETA), 37 heavy atoms, net
  **+1**. Parameters are GAFF2/AM1-BCC, one residue `DAM`, reused from the parent repo
  (`structures/params/damgo.mol2`, `.frcmod`).
- Coordinates come from **10TM chain S**. Template atom names are assigned by
  **molecular-graph isomorphism** (element-matched, 37 atoms / bonds identical; ring-flip
  symmetric maps resolved by RMSD). Hydrogens are placed by local-frame superposition
  from the template.
- Inactive start: DAMGO is moved into 9PXU by superposing 10TM onto 9PXU on the TM1–TM5
  scaffold (scaffold RMSD 1.39 Å). The closest heavy-atom contact is **2.32 Å
  (Gln126 NE2)**, with no other contacts under 2.5 Å; minimisation resolves it. DAMGO
  N1 to Asp149 (D3.32) carboxylate is 4.66 Å, versus 3.86 Å in 10TM. The salt bridge is
  expected to re-form during restrained equilibration; `d_salt` tracks it.

## Na⁺ (sodium site)

9PXU resolves a Na⁺ **2.41 Å** from Asp116 (D2.50) carboxylate oxygen. It is kept, as §8
requires. 10TM has none (active state, collapsed pocket). The topology builder balances
total Na⁺/Cl⁻ so both systems contain the same number of ions; in the active start the
corresponding Na⁺ starts in bulk. **Sodium-site occupancy is a slow degree of freedom
outside the two CVs** and is tracked as `na_site`. If the two starts disagree, it is a
named suspect.

## Membrane frame

Both starts are superposed on the TM1–TM5 scaffold onto **OPM's 8F7Q** (human, canonical
numbering; 9PXU and 10TM are not in OPM). Scaffold RMSD to 8F7Q is 0.51 Å (active) and
1.34 Å (inactive). OPM half-thickness is 17.8 Å.

## §6.5 assertion output

`structures/assert_numbering.out`: all 8 entries canonical P35372 with zero offset,
target residues correct, author chains R,F / R,M / R / A as specified. Exit 0.
