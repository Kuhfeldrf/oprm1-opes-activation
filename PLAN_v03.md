# v0.3 plan: agonist vs antagonist efficacy test with Nb39

Written 2026-10-13, before any v0.3 data existed. Deviations go in NOTES.md.

## Question

Does an agonist (DAMGO) hold the active receptor state better than an antagonist
(naloxone)? This tests efficacy, not the full activation landscape. v0.1/v0.2 showed
that the full landscape is out of reach with the current CVs and system.

## Systems

All four systems share **one set of receptor coordinates**. Within each pair, the
ligand is the only difference.

| System | Receptor | Nb39 | Ligand | Role |
|---|---|---|---|---|
| `nb39_dam` | 5C1M, humanised | yes | DAMGO | agonist, transducer-stabilised |
| `nb39_nlx` | same | yes | naloxone | antagonist, transducer-stabilised |
| `free_dam` | same | removed before packing | DAMGO | agonist, no transducer |
| `free_nlx` | same | removed before packing | naloxone | antagonist, no transducer |

**Receptor: 5C1M chain A** (mouse μOR, BU72 agonist, Nb39, 2.07 Å). It replaces 10TM
for v0.3.
- **Why not 10TM.** 10TM's ICL3 sits about 6 Å from 5C1M's. Placing Nb39 on it gives
  19 heavy atoms within 2 Å of the receptor.
- **Numbering.** Human = mouse + 2, constant over the whole chain (checked by
  alignment).
- **Humanising.** Residues 68–347 contain three mouse/human differences; they are
  mutated to human: N139T, V189I, I308V (human numbers). Residue 158 is already Phe,
  matching the human native residue.
- **Completeness.** No gaps in 68–347, and YCM (59) is outside the range.
- **Termini.** Caps are ACE68 and NME347, as in v0.1.
- **Protonation.** The v0.1 states are reused unchanged, so the receptor's atom order
  matches v0.1 and the CV2 reference indices stay valid. This is asserted at build
  time.
- **Active state.** CV1 is 12.9 Å and CV2 4.6 Å against 9PXU. There is no Na⁺ at
  D2.50.

**Nb39: 5C1M chain B**, VHH residues 3–124, with the AAA tag removed.
- CDR3 residues 103–109 (GQSSSPY) are disordered and get modelled by pdbfixer.
- The Cys24–Cys96 disulfide is kept.
- Termini are free and charged.
- Protonation comes from propka at pH 7.4.
- In `system_ref.pdb`, Nb39 residues are numbered +1000 (chain N) so they cannot
  collide with receptor numbers in `MOLINFO`.

**Ligands**, each placed by TM1–TM5 scaffold superposition. Neither has heavy-atom
contacts under 2.5 Å with the 5C1M receptor.
- **DAMGO** from 10TM. Closest contact 2.67 Å. Existing GAFF2 parameters (`DAM`).
- **Naloxone** from 9PXU. Closest contact 2.93 Å.
  - New GAFF2/AM1-BCC parameters (`NLX`) with the amine protonated (net +1).
  - Hydrogens are added to the crystal heavy atoms, so stereochemistry comes from the
    structure.

Both ligands carry +1, so each pair has identical composition: the same lipids,
waters and ions.

## Protocol

1. **Build.**
   - Pack into POPC:CHL1 7:3 at 150 mM NaCl, using the v0.1 packing settings.
   - Match composition within each pair (`nb39_*`, `free_*`).
   - Same force fields and HMR as v0.1.
2. **Equilibrate** with `systems/equilibrate.sh`: the staged restraint release, then a
   short unbiased leg (10 ns, QC only).
3. **Production.**
   - 4 independent replicas per system, each with fresh velocities and its own seed.
   - 500 ns each, unbiased. Total 4 × 4 × 0.5 = 8 μs.
   - One 4-GPU node per system; the four run concurrently.
   - Observables are logged every 10 ps.
4. **Observables.**
   - The v0.1 set: CV1, CV2, d_hbond, chi_w295, d_salt, na_site, lig_cont.
   - Added for v0.3: Nb39–receptor heavy-atom contacts (`nb_cont`, in the Nb39 systems)
     and ligand RMSD in the pocket after scaffold alignment (`lig_rmsd`).

## Pre-registered readouts and decision rule

Active-state definition: **CV1 ≥ 1.00 nm and CV2 ≥ 0.25 nm**. This is looser than the
v0.1 active box, so that small breathing is not counted as deactivation.

**Primary readout: the no-transducer pair (`free_*`).**
- Per replica, measure the fraction of 100–500 ns spent active and the time of first
  persistent exit. A persistent exit means 10 ns or more outside the active state.
- Compare DAMGO and naloxone by the mean fraction, with SE taken over replicas.
- **Discrimination** if DAMGO − naloxone > 2 SE *and* DAMGO is higher in at least 3 of
  4 replica-rank pairings (Mann–Whitney on the 4 v 4 replica means, one-sided).
- **Inconclusive** if the difference is under 2 SE. Inconclusive is not the same as
  "no effect".

**Secondary readout: the Nb39 pair (`nb39_*`).**
- Nb39 contact retention (`nb_cont` relative to its 0–10 ns mean).
- Fraction of time active, compared in the same way.

The expectation going in: Nb39 holds both systems active on this timescale. A
difference is therefore expected mainly in the interface (`nb_cont`) and in
fluctuations, not in a full exit. So `nb39_*` is mostly a positive control for the
active state.

**Validity gates.** A replica that fails one is reported and excluded; nothing is
silently dropped.
- The ligand stays bound: d_salt ≤ 0.6 nm in 95% of frames, and lig_cont is not
  collapsing.
- The membrane QC passes, with the same thresholds as v0.1.
- Equilibration: the system holds active (both CVs inside the active definition) at
  the end of the restrained stages.

## Out of scope for v0.3

These are recorded so they are not lost:
- **The full activation landscape** (Nb39, plus a path or Deep-TDA CV trained on the
  v0.1/v0.2 data). It needs about 10 μs or more and builds on the `nb39_dam` system
  made here.
- **Binding affinity** (funnel metadynamics or ABFE). It is a different method, and
  affinity alone does not separate agonist from antagonist.
- **External validation** against cAMP, BRET or GTPγS data. It is not tested.

## Stage gates (each committed and pushed)

| Gate | What it covers |
|---|---|
| G1 | structures + naloxone parameters |
| G2 | four built systems with composition matched within each pair |
| G3 | equilibration QC |
| G4 | production done, with monitoring snapshots along the way |
| G5 | analysis + README section + tag `v0.3-efficacy` |
