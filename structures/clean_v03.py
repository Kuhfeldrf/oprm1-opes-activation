"""v0.3 start complexes (PLAN_v03.md): one active receptor, two ligands, ± Nb39.

    receptor  5C1M chain A (mouse; human = mouse + 2), residues 68-347 human, humanised
              N139T, V189I, I308V; ACE68/NME347 caps; v0.1 protonation pattern
    Nb39      5C1M chain B, VHH 3-124 (AAA tag dropped); CDR3 103-109 (GQSSSPY, disordered
              in the crystal) built by pdbfixer, then relaxed with everything else frozen
    DAMGO     from 10TM by TM1-TM5 scaffold superposition (parent template damgo.mol2)
    naloxone  from 9PXU by TM1-TM5 scaffold superposition (structures/params/naloxone.mol2)

Writes structures/clean/{nb39_dam,nb39_nlx,free_dam,free_nlx}/complex_opm.pdb and a
clean_report.json each. The receptor coordinates are identical in all four; within each
pair only the ligand differs.

Order in complex_opm.pdb: receptor (chain R), Na+ none, ligand (chain L), Nb39 (chain N).
Run with the oprm1-opes env python (pdbfixer, openmm); pdb2pqr from mor-pilot.
Usage: python structures/clean_v03.py
"""
import json
import subprocess
import sys
from pathlib import Path

import gemmi
import numpy as np
import openmm as mm
from openmm import app, unit
from pdbfixer import PDBFixer

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "cv"))
from clean import (CANON, FIRST, LAST, OPM, PDB2PQR, ca_of, damgo_onto,  # noqa: E402
                   protonation_states, read_mol2, scaffold_fit)
from scaffold import CV1_PAIR, NPXXYA  # noqa: E402

RAW, OUT, PAR = HERE / "raw", HERE / "clean", HERE / "params"
SRC, REC, NB = "5C1M", "A", "B"
MOUSE_TO_HUMAN = 2
HUMANISE = ["ASN-139-THR", "VAL-189-ILE", "ILE-308-VAL"]
NB_RANGE = (3, 124)                 # VHH domain; 125-127 (AAA) are tag
NB_CDR3_GAP = range(103, 110)       # missing in 5C1M, modelled here
LIGANDS = {"dam": {"pdb": "10TM", "rec": "R", "lig": ("chain", "S"), "mol2": "damgo.mol2",
                   "resn": "DAM", "cation_N": "N1"},
           "nlx": {"pdb": "9PXU", "rec": "R", "lig": ("resname", "A1APV"),
                   "mol2": "naloxone.mol2", "resn": "NLX", "cation_N": None}}
SYSTEMS = {"nb39_dam": ("dam", True), "nb39_nlx": ("nlx", True),
           "free_dam": ("dam", False), "free_nlx": ("nlx", False)}


def read(pid):
    st = gemmi.read_structure(str(RAW / f"{pid}.cif"))
    st.remove_alternative_conformations()
    st.remove_hydrogens()
    st.setup_entities()
    return st


def fixer_from(chain, path, seqres=None):
    """Write one chain to PDB (optionally with SEQRES) and return a PDBFixer on it."""
    s = gemmi.Structure()
    m = gemmi.Model("1")
    m.add_chain(chain)
    s.add_model(m)
    s.setup_entities()
    if seqres:
        s.entities[0].full_sequence = seqres
    s.write_pdb(str(path))
    return PDBFixer(filename=str(path))


def mindist(A, B):
    return float(np.min(np.linalg.norm(np.asarray(A)[:, None] - np.asarray(B)[None], axis=2)))


def main():
    work = OUT / "v03_work"
    work.mkdir(parents=True, exist_ok=True)
    rep0 = {"source": SRC}
    src = read(SRC)

    # ---------------- receptor: renumber to human, cut 68-347, humanise, complete
    rec = gemmi.Chain("R")
    for r in src[0][REC]:
        if r.het_flag != "A":
            continue
        r = r.clone()
        r.seqid.num += MOUSE_TO_HUMAN
        if FIRST <= r.seqid.num <= LAST:
            rec.add_residue(r)
    nums = [r.seqid.num for r in rec]
    gaps = [n for n in range(FIRST, LAST + 1) if n not in nums]
    if gaps:
        raise SystemExit(f"FATAL: receptor gaps {gaps}")
    if any(r.name == "YCM" for r in rec):
        raise SystemExit("FATAL: modified residue (YCM) in range")
    fx = fixer_from(rec, work / "receptor_raw.pdb")
    fx.findMissingResidues()
    fx.missingResidues = {}
    before = {int(r.id): r.name for r in fx.topology.residues()}
    for m_ in HUMANISE:
        wt, n, _ = m_.split("-")
        if before[int(n)] != wt:
            raise SystemExit(f"FATAL: residue {n} is {before[int(n)]}, expected mouse {wt}")
    fx.applyMutations(HUMANISE, "R")
    fx.findMissingAtoms()
    rep0["receptor_completed_sidechains"] = sorted(f"{r.name}{r.id}:{len(v)}"
                                                   for r, v in fx.missingAtoms.items())
    fx.addMissingAtoms()
    with open(work / "receptor_fixed_native.pdb", "w") as f:
        app.PDBFile.writeFile(fx.topology, fx.positions, f, keepIds=True)
    rep0["humanised"] = [f"{before[int(m_.split('-')[1])]}{m_.split('-')[1]}"
                         f"{m_.split('-')[2][:3]}" for m_ in HUMANISE]

    # ---------------- Nb39: VHH 3-124, model the CDR3 gap, relax it alone
    nbent = [e for e in src.entities if src[0][NB].get_polymer().subchain_id()
             in e.subchains][0]
    seq = list(nbent.full_sequence)          # numbering 1..127 in 5C1M
    nb = gemmi.Chain("N")
    for r in src[0][NB]:
        if r.het_flag == "A" and NB_RANGE[0] <= r.seqid.num <= NB_RANGE[1]:
            nb.add_residue(r.clone())
    present = {r.seqid.num for r in nb}
    missing = [n for n in range(NB_RANGE[0], NB_RANGE[1] + 1) if n not in present]
    if missing != list(NB_CDR3_GAP):
        raise SystemExit(f"FATAL: Nb39 gaps {missing}, expected {list(NB_CDR3_GAP)}")
    # author numbering = SEQRES position + constant offset (checked residue by residue)
    offs = {r.seqid.num - r.label_seq for r in nb}
    if len(offs) != 1:
        raise SystemExit(f"FATAL: Nb39 author/SEQRES offset not constant {offs}")
    off = offs.pop()
    for r in nb:
        if seq[r.label_seq - 1] != r.name:
            raise SystemExit(f"FATAL: Nb39 {r.name}{r.seqid.num} does not match SEQRES")
    sub = seq[NB_RANGE[0] - off - 1:NB_RANGE[1] - off]
    rep0["nb39_author_minus_seqres"] = off
    fxn = fixer_from(nb, work / "nb39_raw.pdb", seqres=sub)
    fxn.findMissingResidues()
    nres = len(list(list(fxn.topology.chains())[0].residues()))
    internal = {k: v for k, v in fxn.missingResidues.items() if 0 < k[1] < nres}
    if sum(len(v) for v in internal.values()) != len(NB_CDR3_GAP):
        raise SystemExit(f"FATAL: pdbfixer sees missing residues {fxn.missingResidues}")
    fxn.missingResidues = internal
    rep0["nb39_terminal_missing_ignored"] = {str(k): v for k, v in fxn.missingResidues.items()
                                             if k not in internal} or None
    fxn.findMissingAtoms()
    rep0["nb39_completed_sidechains"] = sorted(f"{r.name}{r.id}:{len(v)}"
                                               for r, v in fxn.missingAtoms.items())
    fxn.addMissingAtoms()
    # pdbfixer renumbers inserted residues sequentially from the chain start
    ids = [int(r.id) for r in fxn.topology.residues()]
    if ids != list(range(NB_RANGE[0], NB_RANGE[1] + 1)):
        # assign SEQRES positions explicitly
        for k, r in enumerate(fxn.topology.residues()):
            r.id = str(NB_RANGE[0] + k)
    names = [r.name for r in fxn.topology.residues()]
    if names != sub:
        raise SystemExit("FATAL: Nb39 sequence after modelling does not match SEQRES")

    # combine receptor + Nb39 to relax the loop against the receptor surface
    mod = app.Modeller(fx.topology, fx.positions)
    mod.add(fxn.topology, fxn.positions)

    class _F:  # minimal object with topology/positions for relax_loop
        topology, positions = mod.topology, mod.positions
    loop = set(range(NB_CDR3_GAP.start - 6, NB_CDR3_GAP.stop + 3))   # CDR3 97-111 incl. Arg99
    # receptor ids 102-110 also exist: freeze by chain, not by residue id
    relaxed = relax_loop_chain(_F, "N", loop, work / "rec_nb39_relaxed_H.pdb")
    nbxyz, (e0, e1), nfree = relaxed
    rep0["nb39_loop_relax"] = {"residues": f"{min(loop)}-{max(loop)}", "free_heavy_atoms": nfree,
                               "E_kJ_before": round(e0, 1), "E_kJ_after": round(e1, 1)}

    # ---------------- OPM frame (receptor scaffold onto OPM 8F7Q)
    ast = gemmi.read_structure(str(work / "receptor_fixed_native.pdb"))
    opm = gemmi.read_structure(str(RAW / OPM[0]))
    R, t, r_opm, n_opm = scaffold_fit(ca_of(ast[0]["R"]), ca_of(opm[0][OPM[1]]))
    rep0["opm_scaffold_rmsd_A"] = round(r_opm, 3)
    tr = gemmi.Transform(gemmi.Mat33(R.tolist()), gemmi.Vec3(*t))
    ast[0].transform_pos_and_adp(tr)
    recf = ast[0]["R"]
    seqh = {r.seqid.num: r.name for r in recf}
    bad = {n: seqh.get(n) for n, v in CANON.items() if seqh.get(n) != v}
    if bad or seqh[158] != "PHE":
        raise SystemExit(f"FATAL: canonical residue check failed {bad} / F158 {seqh[158]}")
    nb_atoms = [(rid, nm, R @ x + t) for (rid, nm), x in nbxyz.items()]
    rec_xyz = np.array([a.pos.tolist() for r in recf for a in r])
    rep0["nb39_receptor_min_heavy_dist_A"] = round(mindist([x for *_, x in nb_atoms], rec_xyz), 2)

    # CVs of the start (Angstrom; CV2 against deposited 9PXU as in v0.1)
    ref = ca_of(read("9PXU")[0]["R"])
    ca = ca_of(recf)
    Rr, tr_, _, _ = scaffold_fit(ca, ref)
    X = np.array([ca[n] for n in NPXXYA]) @ Rr.T + tr_
    Y = np.array([ref[n] for n in NPXXYA])
    rep0["CV1_CA_A"] = round(float(np.linalg.norm(ca[CV1_PAIR[0]] - ca[CV1_PAIR[1]])), 2)
    rep0["CV2_vs_9PXU_A"] = round(float(np.sqrt(((X - Y) ** 2).sum(1).mean())), 2)

    # ---------------- protonation: receptor = v0.1 pattern; Nb39 = propka
    v01 = json.loads((OUT / "active" / "clean_report.json").read_text())
    rec_states = {}
    for k, v in v01["final_states_nondefault"].items():
        rec_states[int("".join(c for c in k if c.isdigit()))] = v
    rep0["receptor_states_from_v01"] = v01["final_states_nondefault"]
    nb_only = work / "nb39_frame.pdb"
    lines = []
    for i, (rid, nm, x) in enumerate(sorted(nb_atoms, key=lambda a: a[0]), 1):
        resn = sub[rid - NB_RANGE[0]]
        lines.append(f"ATOM  {i:5d} {nm if len(nm) == 4 else ' ' + nm:<4s} {resn:>3s} N"
                     f"{rid:4d}    {x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}  1.00  0.00"
                     f"          {nm[0]:>2s}")
    nb_only.write_text("\n".join(lines) + "\nTER\nEND\n")
    r_ = subprocess.run([PDB2PQR, "--ff", "AMBER", "--titration-state-method", "propka",
                         "--with-ph", "7.4", "--keep-chain", "--pdb-output",
                         str(work / "nb39_pqr.pdb"), str(nb_only), str(work / "nb39.pqr")],
                        capture_output=True, text=True)
    (work / "nb39_pdb2pqr.log").write_text(r_.stdout + r_.stderr)
    if r_.returncode != 0:
        raise SystemExit("FATAL: pdb2pqr failed on Nb39")
    nb_states = protonation_states(work / "nb39_pqr.pdb")
    nb_cys = {rid: x for rid, nm, x in nb_atoms if nm == "SG"}
    nb_ss = [(i, j) for i in nb_cys for j in nb_cys if i < j
             and np.linalg.norm(nb_cys[i] - nb_cys[j]) < 2.5]
    if nb_ss != [(24, 96)]:
        raise SystemExit(f"FATAL: Nb39 disulfide {nb_ss}, expected [(24, 96)]")
    rep0["nb39_disulfides"] = nb_ss
    rep0["nb39_states_nondefault"] = {f"{sub[n - NB_RANGE[0]]}{n}": v for n, v in nb_states.items()
                                      if v not in ("CYS", sub[n - NB_RANGE[0]])
                                      or sub[n - NB_RANGE[0]] == "HIS"}
    # receptor disulfide
    sg = {r.seqid.num: np.array(r["SG"][0].pos.tolist()) for r in recf if r.name == "CYS"}
    rec_ss = [(i, j) for i in sg for j in sg if i < j and np.linalg.norm(sg[i] - sg[j]) < 2.5]
    if rec_ss != [(142, 219)]:
        raise SystemExit(f"FATAL: receptor disulfide {rec_ss}")

    # ---------------- ligands into the final receptor frame
    ligs = {}
    for key, cfg in LIGANDS.items():
        s = read(cfg["pdb"])
        mode, sel = cfg["lig"]
        if mode == "chain":
            heavy = [(a.element.name.upper(), np.array(a.pos.tolist())) for r in s[0][sel]
                     for a in r]
        else:
            heavy = [(a.element.name.upper(), np.array(a.pos.tolist()))
                     for r in s[0][cfg["rec"]] if r.name == sel for a in r]
        Rl, tl, r_tx, _ = scaffold_fit(ca_of(s[0][cfg["rec"]]), ca)
        tgt = [(el, Rl @ x + tl) for el, x in heavy]
        atoms, bonds = read_mol2(PAR / cfg["mol2"])
        placed, map_rmsd, _ = damgo_onto(atoms, bonds, tgt)
        lh = np.array([x for _, el, x in placed if el != "H"])
        labs = [f"{r.name}{r.seqid.num}:{a.name}" for r in recf for a in r]
        D = np.linalg.norm(lh[:, None] - rec_xyz[None], axis=2)
        close = sorted({labs[j] for i, j in zip(*np.where(D < 2.5))})
        # cationic N: the only N bonded to 4 atoms (amine, protonated) in the template
        deg = np.bincount(np.array(bonds).ravel(), minlength=len(atoms))
        catN = [a["name"] for i, a in enumerate(atoms) if a["el"] == "N" and deg[i] == 4]
        if len(catN) != 1:
            raise SystemExit(f"FATAL: {key}: cationic N ambiguous {catN}")
        n_xyz = [x for nm, _, x in placed if nm == catN[0]][0]
        d149 = min(np.linalg.norm(np.array(a.pos.tolist()) - n_xyz) for r in recf
                   if r.seqid.num == 149 for a in r if a.name.startswith("OD"))
        ligs[key] = {"atoms": placed, "resn": cfg["resn"], "report": {
            "source": cfg["pdb"], "scaffold_rmsd_A": round(r_tx, 3),
            "template_map_rmsd_A": round(map_rmsd, 3),
            "min_heavy_dist_to_receptor_A": round(float(D.min()), 2),
            "receptor_contacts_lt_2.5A": close, "cationic_N": catN[0],
            "cationic_N_to_D149_OD_min_A": round(float(d149), 2),
            "min_heavy_dist_to_nb39_A": round(mindist(lh, [x for *_, x in nb_atoms]), 2)}}

    # ---------------- write the four complexes
    for name, (lk, with_nb) in SYSTEMS.items():
        out = OUT / name
        out.mkdir(parents=True, exist_ok=True)
        lines, serial = [], 0

        def atom(nm, resn, chn, resi, xyz, el, het=False):
            nonlocal serial
            serial += 1
            n4 = nm if len(nm) == 4 else f" {nm:<3}"
            return ("HETATM" if het else "ATOM  ") + f"{serial:5d} {n4:4s} {resn:>3s} {chn}" \
                   f"{resi:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00" \
                   f"          {el:>2s}"
        for r in recf:
            n = r.seqid.num
            if n == FIRST:
                for a in r:
                    if a.name in ("CA", "C", "O"):
                        lines.append(atom("CH3" if a.name == "CA" else a.name, "ACE", "R", n,
                                          a.pos.tolist(), a.element.name))
                continue
            if n == LAST:
                for a in r:
                    if a.name in ("N", "CA"):
                        lines.append(atom("C" if a.name == "CA" else a.name, "NME", "R", n,
                                          a.pos.tolist(), a.element.name))
                continue
            resn = "CYX" if n in (142, 219) else rec_states.get(n, r.name)
            for a in r:
                lines.append(atom(a.name, resn, "R", n, a.pos.tolist(), a.element.name))
        lines.append("TER")
        for nm, el, x in ligs[lk]["atoms"]:
            lines.append(atom(nm, ligs[lk]["resn"], "L", 1, x, el, het=True))
        lines.append("TER")
        if with_nb:
            for rid in range(NB_RANGE[0], NB_RANGE[1] + 1):
                resn = sub[rid - NB_RANGE[0]]
                resn = "CYX" if rid in (24, 96) else nb_states.get(rid, resn)
                if resn == "CYS":
                    resn = "CYX" if rid in (24, 96) else "CYS"
                for r_id, nm, x in nb_atoms:
                    if r_id == rid:
                        lines.append(atom(nm, resn, "N", rid, x, nm[0]))
            lines.append("TER")
        lines.append("END")
        (out / "complex_opm.pdb").write_text("\n".join(lines) + "\n")
        rep = dict(rep0, system=name, ligand=ligs[lk]["report"], nb39=with_nb,
                   ligand_resname=ligs[lk]["resn"])
        z = np.array([float(l[46:54]) for l in lines if l.startswith(("ATOM", "HETATM"))])
        rep["z_range_A"] = [round(float(z.min()), 1), round(float(z.max()), 1)]
        (out / "clean_report.json").write_text(json.dumps(rep, indent=1, default=str))
        print(name, json.dumps({k: rep[k] for k in ("ligand", "z_range_A")}, default=str))
    print(json.dumps(rep0, indent=1, default=str))


def relax_loop_chain(fobj, chain_id, loop_resids, out_pdb):
    """Relax the modelled Nb39 CDR3 with everything else frozen (massless).

    pdbfixer places inserted loops crudely; here the loop was threaded through Arg99, so
    a plain minimisation stalls. Stage 1 replaces all non-bonded terms by a bounded soft
    repulsion (atoms can pass through each other); stage 2 is full amber14 + GBn2. Free:
    heavy atoms and H of `loop_resids` on the Nb39 chain, plus every hydrogen (only heavy
    atoms are kept downstream). Afterwards L-chirality at every CA and trans peptide
    bonds in the loop are asserted, because stage 1 can invert centres.
    Returns Nb39 heavy-atom positions (Angstrom) keyed by (resid, name), (E_start,
    E_end) in kJ/mol for the full potential, and the number of free heavy atoms."""
    mod = app.Modeller(fobj.topology, fobj.positions)
    ff = app.ForceField("amber14-all.xml", "implicit/gbn2.xml")
    mod.addHydrogens(ff, pH=7.4)
    nb_chain = list(mod.topology.chains())[-1]
    atoms = list(mod.topology.atoms())
    free = [at.residue.chain is nb_chain and int(at.residue.id) in loop_resids
            or at.element.symbol == "H" for at in atoms]
    nfree_heavy = sum(f and at.element.symbol != "H" for f, at in zip(free, atoms))

    # L-handedness (improper N-CA-C-CB) and trans omega restraints on the free residues:
    # the soft stage otherwise inverts centres (seen at 105, 108, 111). The improper target
    # is measured on the frozen crystal residues, which are all L.
    P0 = np.array(mod.positions.value_in_unit(unit.nanometer))
    res_at = {}
    for at in atoms:
        if at.residue.chain is nb_chain:
            res_at.setdefault(int(at.residue.id), {})[at.name] = at.index

    def dihedral(i, j, k, l):
        b0, b1, b2 = P0[i] - P0[j], P0[k] - P0[j], P0[l] - P0[k]
        b1 = b1 / np.linalg.norm(b1)
        v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
        return float(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))
    imp = lambda a: (a["N"], a["CA"], a["C"], a["CB"])
    ref_imp = [dihedral(*imp(a)) for rid, a in res_at.items()
               if rid not in loop_resids and "CB" in a]
    if np.ptp(np.sign(ref_imp)) != 0:
        raise SystemExit("FATAL: frozen Nb39 residues have mixed improper sign")
    imp0 = float(np.arctan2(np.mean(np.sin(ref_imp)), np.mean(np.cos(ref_imp))))

    def restraints():
        f = mm.CustomTorsionForce("k*(1-cos(theta-theta0))")
        f.addPerTorsionParameter("k")
        f.addPerTorsionParameter("theta0")
        for rid in sorted(loop_resids):
            a = res_at.get(rid, {})
            if "CB" in a:
                f.addTorsion(*imp(a), [500.0, imp0])
            b = res_at.get(rid + 1, {})
            if a and b and "N" in b:
                f.addTorsion(a["CA"], a["C"], b["N"], b["CA"], [500.0, np.pi])
        return f

    def build(soft):
        system = ff.createSystem(mod.topology, nonbondedMethod=app.NoCutoff, constraints=None)
        system.addForce(restraints())
        if soft:
            nbf = [f for f in system.getForces() if isinstance(f, mm.NonbondedForce)][0]
            rep = mm.CustomNonbondedForce("eps*step(sig-r)*(1-r/sig)^2; sig=0.25; eps=500")
            rep.setNonbondedMethod(mm.CustomNonbondedForce.NoCutoff)
            for _ in atoms:
                rep.addParticle([])
            for k in range(nbf.getNumExceptions()):
                i, j, *_ = nbf.getExceptionParameters(k)
                rep.addExclusion(i, j)
            for k in reversed(range(system.getNumForces())):
                if isinstance(system.getForce(k), (mm.NonbondedForce, mm.CustomGBForce,
                                                    mm.GBSAOBCForce)):
                    system.removeForce(k)
            system.addForce(rep)
        for at, f in zip(atoms, free):
            if not f:
                system.setParticleMass(at.index, 0.0)
        return system

    full = build(False)
    ctx_full = mm.Context(full, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("CPU"))
    ctx_full.setPositions(mod.positions)
    e0 = ctx_full.getState(getEnergy=True).getPotentialEnergy()
    soft = build(True)     # reported energies include the (small) restraint term
    ctx = mm.Context(soft, mm.VerletIntegrator(0.001), mm.Platform.getPlatformByName("CPU"))
    ctx.setPositions(mod.positions)
    mm.LocalEnergyMinimizer.minimize(ctx, 1.0, 10000)
    ctx_full.setPositions(ctx.getState(getPositions=True).getPositions())
    mm.LocalEnergyMinimizer.minimize(ctx_full, 1.0, 20000)
    st = ctx_full.getState(getEnergy=True, getPositions=True)
    with open(out_pdb, "w") as f:
        app.PDBFile.writeFile(mod.topology, st.getPositions(), f, keepIds=True)
    pos = st.getPositions(asNumpy=True).value_in_unit(unit.angstrom)

    # sanity: CA chirality (L) everywhere on Nb39, trans omega inside the loop
    byres = {}
    for at in atoms:
        if at.residue.chain is nb_chain:
            byres.setdefault(int(at.residue.id), {})[at.name] = pos[at.index]
    # handedness of every CA, compared with the frozen crystal residues (all L)
    sgn = {rid: np.sign(np.dot(np.cross(a["N"] - a["CA"], a["C"] - a["CA"]), a["CB"] - a["CA"]))
           for rid, a in byres.items() if "CB" in a}
    ref = {v for rid, v in sgn.items() if rid not in loop_resids}
    if len(ref) != 1:
        raise SystemExit(f"FATAL: frozen Nb39 residues have mixed CA handedness {ref}")
    flipped = [rid for rid, v in sgn.items() if v not in ref]
    if flipped:
        raise SystemExit(f"FATAL: Nb39 CA chirality inverted by loop relax at {flipped}")

    def dih(p0, p1, p2, p3):
        b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
        b1 /= np.linalg.norm(b1)
        v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
        return np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))
    for rid in sorted(loop_resids):
        if rid + 1 in byres:
            a, b = byres[rid], byres[rid + 1]
            om = dih(a["CA"], a["C"], b["N"], b["CA"])
            if abs(om) < 150:
                raise SystemExit(f"FATAL: Nb39 omega {rid}-{rid + 1} = {om:.0f} deg after relax")
    xyz = {(int(a.residue.id), a.name): pos[a.index] for a in atoms
           if a.residue.chain is nb_chain and a.element.symbol != "H"}
    return xyz, (e0.value_in_unit(unit.kilojoule_per_mole),
                 st.getPotentialEnergy().value_in_unit(unit.kilojoule_per_mole)), nfree_heavy


if __name__ == "__main__":
    main()
