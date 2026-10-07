"""Runbook §4.4: build the two Option-A start complexes (DAMGO-bound μOR).

    active    10TM chain R + DAMGO (chain S)              revert F158W
    inactive  9PXU chain R + DAMGO transplanted from 10TM  revert L266M, R271K
              (+ the resolved Na+ at the D2.50 site)

Both starts carry the IDENTICAL sequence 69-346, with residue 68 -> ACE and 347 -> NME
caps, so the two systems differ only in conformation. Everything not listed is
deleted: BRIL, nanobody, Fab, naloxone, cholesterol, the second 10TM receptor copy.

Outputs per start (in --out/<start>/):
    receptor_fixed.pdb   pdbfixer-completed, reverted, OPM frame, no H, canonical numbering
    receptor_pqr.pdb     pdb2pqr/propka pH 7.4 (hydrogens used only to read states)
    complex_opm.pdb      ACE/NME capped, Amber residue names, no protein H, + Na+ + DAM (with H)
    clean_report.json    every deletion, reversion, protonation state, clash and CV value

Run with the oprm1-opes env python; pdb2pqr is called from the mor-pilot env.
Usage: python structures/clean.py --raw structures/raw --out structures/clean
"""
import argparse
import itertools
import json
import subprocess
import sys
from pathlib import Path

import gemmi
import networkx as nx
import numpy as np
from networkx.algorithms import isomorphism
from openmm.app import PDBFile
from pdbfixer import PDBFixer

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cv"))
from scaffold import SCAFFOLD, NPXXYA, CV1_PAIR  # noqa: E402

PDB2PQR = "/scratch/kuhfeldr-Kuhfeld_temp/miniforge3/envs/mor-pilot/bin/pdb2pqr"
FIRST, LAST = 68, 347          # 68 -> ACE, 347 -> NME, modelled residues 69-346
STARTS = {
    "active":   {"pdb": "10TM", "chain": "R", "mutations": ["TRP-158-PHE"], "ions": []},
    "inactive": {"pdb": "9PXU", "chain": "R", "mutations": ["LEU-266-MET", "ARG-271-LYS"],
                 "ions": ["NA"]},
}
DAMGO_SRC = ("10TM", "R", "S")   # pdb, receptor chain it binds, ligand chain
OPM = ("8f7q_opm.pdb", "R")      # membrane frame: OPM-oriented 8F7Q, canonical numbering
CANON = {166: "ASP", 167: "ARG", 168: "TYR", 281: "THR", 334: "ASN", 335: "PRO",
         336: "VAL", 337: "LEU", 338: "TYR", 339: "ALA"}


def kabsch(P, Q):
    pc, qc = P.mean(0), Q.mean(0)
    U, _, Vt = np.linalg.svd((P - pc).T @ (Q - qc))
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1.0, 1.0, d]) @ U.T
    return R, qc - R @ pc


def ca_of(chain):
    return {r.seqid.num: np.array(r["CA"][0].pos.tolist())
            for r in chain if r.het_flag == "A" and r.find_atom("CA", "*")}


def scaffold_fit(mobile_ca, target_ca):
    common = [n for n in SCAFFOLD if n in mobile_ca and n in target_ca]
    P = np.array([mobile_ca[n] for n in common]); Q = np.array([target_ca[n] for n in common])
    R, t = kabsch(P, Q)
    rmsd = float(np.sqrt((((P @ R.T + t) - Q) ** 2).sum(1).mean()))
    return R, t, rmsd, len(common)


def read_mol2(path):
    txt = Path(path).read_text()
    atoms_blk = txt.split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0].strip().splitlines()
    bonds_blk = txt.split("@<TRIPOS>BOND")[1].split("@<TRIPOS>")[0].strip().splitlines()
    atoms = []
    for l in atoms_blk:
        f = l.split()
        el = "".join(c for c in f[1] if c.isalpha())[:1].upper()
        atoms.append({"name": f[1], "xyz": np.array([float(x) for x in f[2:5]]),
                      "type": f[5], "el": el, "q": float(f[8])})
    bonds = [(int(l.split()[1]) - 1, int(l.split()[2]) - 1) for l in bonds_blk]
    return atoms, bonds


def damgo_onto(mol2_atoms, mol2_bonds, target):
    """Map the GAFF2 template onto target heavy atoms by graph isomorphism.

    target: list of (element, xyz). Returns list of (name, xyz) for ALL template atoms,
    heavy atoms at target coordinates, hydrogens placed by local-frame superposition."""
    G = nx.Graph()
    heavy = [i for i, a in enumerate(mol2_atoms) if a["el"] != "H"]
    for i in heavy:
        G.add_node(i, el=mol2_atoms[i]["el"])
    G.add_edges_from((i, j) for i, j in mol2_bonds if i in G and j in G)
    T = nx.Graph()
    X = np.array([x for _, x in target])
    for k, (el, _) in enumerate(target):
        T.add_node(k, el=el)
    for a, b in itertools.combinations(range(len(target)), 2):
        if np.linalg.norm(X[a] - X[b]) < 1.75:
            T.add_edge(a, b)
    if G.number_of_nodes() != T.number_of_nodes() or G.number_of_edges() != T.number_of_edges():
        raise SystemExit(f"FATAL: DAMGO graphs differ: template {G.number_of_nodes()} atoms/"
                         f"{G.number_of_edges()} bonds, structure {T.number_of_nodes()}/"
                         f"{T.number_of_edges()}")
    gm = isomorphism.GraphMatcher(G, T, node_match=lambda a, b: a["el"] == b["el"])
    best = None
    tmpl = np.array([mol2_atoms[i]["xyz"] for i in heavy])
    for m in gm.isomorphisms_iter():          # ring flips give several equivalent maps
        Y = np.array([X[m[i]] for i in heavy])
        R, t = kabsch(tmpl, Y)
        r = float(np.sqrt((((tmpl @ R.T + t) - Y) ** 2).sum(1).mean()))
        if best is None or r < best[0]:
            best = (r, m)
    if best is None:
        raise SystemExit("FATAL: no element-preserving isomorphism between template and DAMGO")
    rmsd, m = best
    pos = {i: X[m[i]] for i in heavy}
    full = nx.Graph(); full.add_edges_from(mol2_bonds)
    for h, a in enumerate(mol2_atoms):
        if a["el"] != "H":
            continue
        (p,) = list(full.neighbors(h))
        frame = [p] + [n for n in full.neighbors(p) if n in pos]
        for n in list(frame[1:]):
            frame += [k for k in full.neighbors(n) if k in pos and k not in frame]
        frame = frame[:6]
        R, t = kabsch(np.array([mol2_atoms[i]["xyz"] for i in frame]),
                      np.array([pos[i] for i in frame]))
        pos[h] = mol2_atoms[h]["xyz"] @ R.T + t
    return [(a["name"], a["el"], pos[i]) for i, a in enumerate(mol2_atoms)], rmsd, m


def protonation_states(pqr_pdb):
    """Read states from the hydrogens pdb2pqr placed (names alone are not diagnostic)."""
    st = gemmi.read_structure(str(pqr_pdb))
    states = {}
    for r in st[0][0] if len(st[0]) == 1 else [r for ch in st[0] for r in ch]:
        names = {a.name for a in r}
        n = r.seqid.num
        if r.name in ("HIS", "HID", "HIE", "HIP"):
            hd1, he2 = "HD1" in names, "HE2" in names
            states[n] = "HIP" if hd1 and he2 else "HID" if hd1 else "HIE" if he2 else None
            if states[n] is None:
                raise SystemExit(f"FATAL: His{n} has neither HD1 nor HE2")
        elif r.name in ("ASP", "ASH"):
            states[n] = "ASH" if "HD2" in names else "ASP"
        elif r.name in ("GLU", "GLH"):
            states[n] = "GLH" if "HE2" in names else "GLU"
        elif r.name in ("LYS", "LYN"):
            states[n] = "LYS" if "HZ3" in names else "LYN"
        elif r.name in ("CYS", "CYX", "CYM"):
            states[n] = "CYS"
    return states


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--raw", type=Path, default=HERE / "raw")
    ap.add_argument("--out", type=Path, default=HERE / "clean")
    ap.add_argument("--mol2", type=Path, default=HERE / "params" / "damgo.mol2")
    a = ap.parse_args()

    opm = gemmi.read_structure(str(a.raw / OPM[0]))
    opm_ca = ca_of(opm[0][OPM[1]])
    mol2_atoms, mol2_bonds = read_mol2(a.mol2)

    # DAMGO heavy atoms in the 10TM frame, plus 10TM's receptor CA (to transplant)
    src = gemmi.read_structure(str(a.raw / f"{DAMGO_SRC[0]}.cif"))
    src.remove_alternative_conformations()
    dam_src = [(at.element.name.upper(), np.array(at.pos.tolist()))
               for r in src[0][DAMGO_SRC[2]] for at in r if at.element.name != "H"]
    src_ca = ca_of(src[0][DAMGO_SRC[1]])

    work = {}
    for start, cfg in STARTS.items():
        out = a.out / start
        out.mkdir(parents=True, exist_ok=True)
        rep = {"start": start, "source": cfg["pdb"], "chain": cfg["chain"]}
        st = gemmi.read_structure(str(a.raw / f"{cfg['pdb']}.cif"))
        st.remove_alternative_conformations()
        st.remove_hydrogens()
        model = st[0]
        rep["deleted_chains"] = sorted({ch.name for ch in model if ch.name != cfg["chain"]})
        rec = model[cfg["chain"]]
        # the receptor chain also holds hetero groups (BRIL is polymer at negative numbers)
        het = sorted({f"{r.name}{r.seqid.num}" for r in rec if r.het_flag == "H"})
        outside = [r.seqid.num for r in rec if r.het_flag == "A"
                   and not FIRST <= r.seqid.num <= LAST]
        ions = [np.array(r[0].pos.tolist()) for r in rec
                if r.het_flag == "H" and r.name in cfg["ions"]]
        rep["deleted_hetero_in_receptor_chain"] = [h for h in het
                                                   if not any(h.startswith(i) for i in cfg["ions"])]
        rep["deleted_residues_outside_range"] = (f"{min(outside)}..{max(outside)} "
                                                 f"({len(outside)} residues)" if outside else "none")
        keep = gemmi.Structure()
        keep.cell = st.cell
        m = gemmi.Model("1")
        ch = gemmi.Chain("R")
        for r in rec:
            if r.het_flag == "A" and FIRST <= r.seqid.num <= LAST:
                ch.add_residue(r)
        m.add_chain(ch)
        keep.add_model(m)
        nums = [r.seqid.num for r in ch]
        gaps = [n for n in range(FIRST, LAST + 1) if n not in nums]
        if gaps:
            raise SystemExit(f"FATAL: {start}: residues missing in {FIRST}-{LAST}: {gaps}")
        raw_pdb = out / "receptor_raw.pdb"
        keep.setup_entities()
        keep.write_pdb(str(raw_pdb))

        # ---- revert engineered mutations, complete side chains (no hydrogens)
        fx = PDBFixer(filename=str(raw_pdb))
        fx.findMissingResidues()
        fx.missingResidues = {}
        before = {int(r.id): r.name for r in fx.topology.residues()}
        for mut in cfg["mutations"]:
            wt_now, num, _ = mut.split("-")
            if before[int(num)] != wt_now:
                raise SystemExit(f"FATAL: {start} residue {num} is {before[int(num)]}, "
                                 f"expected engineered {wt_now}")
        fx.applyMutations(cfg["mutations"], "R")
        fx.findMissingAtoms()
        rep["completed_sidechains"] = sorted(
            f"{r.name}{r.id}:{len(v)}" for r, v in fx.missingAtoms.items())
        fx.addMissingAtoms()
        after = {int(r.id): r.name for r in fx.topology.residues()}
        rep["reverted"] = [f"{before[int(m_.split('-')[1])]}{m_.split('-')[1]}"
                           f"->{after[int(m_.split('-')[1])]}" for m_ in cfg["mutations"]]
        fixed_pdb = out / "receptor_fixed_native.pdb"
        with open(fixed_pdb, "w") as f:
            PDBFile.writeFile(fx.topology, fx.positions, f, keepIds=True)

        # ---- OPM frame by scaffold superposition onto OPM-oriented 8F7Q
        fst = gemmi.read_structure(str(fixed_pdb))
        R, t, r_opm, n_opm = scaffold_fit(ca_of(fst[0]["R"]), opm_ca)
        rep["opm_scaffold_rmsd_A"] = round(r_opm, 3)
        rep["opm_scaffold_atoms"] = n_opm
        tr = gemmi.Transform(gemmi.Mat33(R.tolist()), gemmi.Vec3(*t))
        fst[0].transform_pos_and_adp(tr)
        fst.write_pdb(str(out / "receptor_fixed.pdb"))
        ions = [R @ x + t for x in ions]
        rec_final = fst[0]["R"]
        seq = {r.seqid.num: r.name for r in rec_final}
        bad = {n: seq.get(n) for n, v in CANON.items() if seq.get(n) != v}
        if bad:
            raise SystemExit(f"FATAL: {start}: canonical residue check failed {bad}")

        # ---- DAMGO into this receptor's frame (transplant via scaffold for 9PXU)
        Rl, tl, r_tx, _ = scaffold_fit(src_ca, ca_of(rec_final))
        rep["damgo_transplant_scaffold_rmsd_A"] = round(r_tx, 3)
        dam_target = [(el, Rl @ x + tl) for el, x in dam_src]
        dam, map_rmsd, _ = damgo_onto(mol2_atoms, mol2_bonds, dam_target)
        rep["damgo_template_vs_10TM_heavy_rmsd_A"] = round(map_rmsd, 3)
        rec_heavy = np.array([at.pos.tolist() for r in rec_final for at in r])
        rec_lab = [f"{r.name}{r.seqid.num}:{at.name}" for r in rec_final for at in r]
        dh = np.array([x for _, el, x in dam if el != "H"])
        D = np.linalg.norm(dh[:, None] - rec_heavy[None], axis=2)
        rep["damgo_min_heavy_dist_A"] = round(float(D.min()), 2)
        close = sorted({rec_lab[j] for i, j in zip(*np.where(D < 2.5))})
        rep["damgo_heavy_contacts_lt_2.5A"] = close
        d149 = [k for k, l in enumerate(rec_lab) if l.startswith("ASP149:OD")]
        n1 = [x for (nm, el, x) in dam if nm == "N1"][0]
        rep["damgo_N1_to_D149_OD_min_A"] = round(float(min(
            np.linalg.norm(rec_heavy[k] - n1) for k in d149)), 2)
        if ions:
            d116 = [k for k, l in enumerate(rec_lab) if l.startswith("ASP116:OD")]
            rep["na_to_D116_OD_min_A"] = [round(float(min(np.linalg.norm(rec_heavy[k] - x)
                                                          for k in d116)), 2) for x in ions]

        # ---- protonation at pH 7.4 (propka via pdb2pqr), states read from hydrogens
        pqr = out / "receptor.pqr"
        r_ = subprocess.run([PDB2PQR, "--ff", "AMBER", "--titration-state-method", "propka",
                             "--with-ph", "7.4", "--keep-chain", "--pdb-output",
                             str(out / "receptor_pqr.pdb"), str(out / "receptor_fixed.pdb"),
                             str(pqr)], capture_output=True, text=True)
        (out / "pdb2pqr.log").write_text(r_.stdout + r_.stderr)
        if r_.returncode != 0:
            raise SystemExit(f"FATAL: pdb2pqr failed for {start}; see {out}/pdb2pqr.log")
        states = protonation_states(out / "receptor_pqr.pdb")
        rep["propka_D116_D2.50"] = states.get(116)
        rep["propka_H299_H6.52"] = states.get(299)
        rep["propka_states"] = states
        work[start] = (rep, rec_final, ions, dam, seq, out, fst)

    # ---- one protonation pattern for BOTH starts: they are the same chemical system,
    # and any composition difference between starts becomes a free-energy artefact (§8).
    # Inactive-start (9PXU, Na+-bound) assignments are used, D2.50 forced to charged ASP
    # because the Na+ it coordinates is kept.
    final_states = dict(work["inactive"][0]["propka_states"])
    final_states[116] = "ASP"
    disagree = {n: {k: work[k][0]["propka_states"].get(n) for k in work}
                for n in final_states
                if len({work[k][0]["propka_states"].get(n) for k in work}) > 1}
    for start in STARTS:
        rep, rec_final, ions, dam, seq, out, fst = work[start]
        states = final_states
        rep["propka_disagreements_between_starts"] = {f"{seq[n]}{n}": v for n, v in disagree.items()}
        rep["final_states_nondefault"] = {f"{seq[n]}{n}": v for n, v in final_states.items()
                                          if v not in ("CYS", seq[n]) or seq[n] == "HIS"}
        del rep["propka_states"]

        # ---- disulfides
        sg = {r.seqid.num: np.array(r["SG"][0].pos.tolist()) for r in rec_final
              if r.name == "CYS" and r.find_atom("SG", "*")}
        ss = [(i, j) for i, j in itertools.combinations(sorted(sg), 2)
              if np.linalg.norm(sg[i] - sg[j]) < 2.5]
        rep["disulfides"] = ss
        ss_res = {x for p in ss for x in p}

        # ---- write complex: ACE68, 69-346 Amber names, NME347, Na+, DAM
        lines, serial = [], 0
        def atom(name, resn, chn, resi, xyz, el, het=False):
            nonlocal serial
            serial += 1
            nm = name if len(name) == 4 else f" {name:<3}"
            return ("HETATM" if het else "ATOM  ") + f"{serial:5d} {nm:4s} {resn:>3s} {chn}" \
                   f"{resi:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}  1.00  0.00" \
                   f"          {el:>2s}"
        for r in rec_final:
            n = r.seqid.num
            if n == FIRST:
                for at in r:
                    if at.name in ("CA", "C", "O"):
                        lines.append(atom("CH3" if at.name == "CA" else at.name, "ACE", "R", n,
                                          at.pos.tolist(), at.element.name))
                continue
            if n == LAST:
                for at in r:
                    if at.name in ("N", "CA"):
                        # ff19SB (aminoct12.lib) names the NME methyl carbon "C", not "CH3"
                        lines.append(atom("C" if at.name == "CA" else at.name, "NME", "R", n,
                                          at.pos.tolist(), at.element.name))
                continue
            resn = r.name
            if n in ss_res:
                resn = "CYX"
            elif n in states and resn != "CYS":
                resn = states[n]
            for at in r:
                lines.append(atom(at.name, resn, "R", n, at.pos.tolist(), at.element.name))
        lines.append("TER")
        for k, x in enumerate(ions):
            lines.append(atom("Na+", "Na+", "I", 1 + k, x, "Na", het=True))
            lines.append("TER")
        for nm, el, x in dam:
            lines.append(atom(nm, "DAM", "L", 1, x, el, het=True))
        lines += ["TER", "END"]
        (out / "complex_opm.pdb").write_text("\n".join(lines) + "\n")

        # ---- CVs on the cleaned start (Angstrom; CV2 vs deposited 9PXU, as in NOTES.md)
        ref = gemmi.read_structure(str(a.raw / "9PXU.cif"))[0]["R"]
        ref_ca, ca = ca_of(ref), ca_of(rec_final)
        Rr, tr_, _, _ = scaffold_fit(ca, ref_ca)
        X = np.array([ca[n] for n in NPXXYA]) @ Rr.T + tr_
        Y = np.array([ref_ca[n] for n in NPXXYA])
        rep["CV1_CA_A"] = round(float(np.linalg.norm(ca[CV1_PAIR[0]] - ca[CV1_PAIR[1]])), 2)
        rep["CV2_vs_9PXU_A"] = round(float(np.sqrt(((X - Y) ** 2).sum(1).mean())), 2)
        z = np.array([ca[n][2] for n in ca])
        rep["receptor_CA_z_range_A"] = [round(float(z.min()), 1), round(float(z.max()), 1)]
        (out / "clean_report.json").write_text(json.dumps(rep, indent=1, default=str))
        print(json.dumps(rep, indent=1, default=str))


if __name__ == "__main__":
    main()
