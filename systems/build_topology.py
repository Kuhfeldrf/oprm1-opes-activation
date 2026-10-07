"""Build GROMACS topologies for both starts with IDENTICAL composition (§8).

Input:  systems/<start>/packed.pdb   (packmol-memgen output)
Output: systems/<start>/{system.top, system.gro, system_ref.pdb, index.ndx, posre_*.itp,
        build_report.json}

What this enforces, each of which otherwise produces a system that runs and is wrong:
  * composition: the two starts get exactly the same number of POPC, CHL, water, Na+,
    Cl- (bulk lipids/waters/ions furthest from the receptor are removed, bulk waters are
    converted to ions). Any compositional difference between starts would appear as a
    free-energy difference between them.
  * ordering: receptor, then DAMGO, then everything else, so protein+ligand occupy
    topology indices 1..N contiguously (WHOLEMOLECULES, posres, PLUMED groups).
  * numbering: packmol-memgen and tleap renumber the receptor from 1. system_ref.pdb is
    written with CANONICAL P35372 numbers (ACE68..NME347) and checked against
    D166/R167/Y168/T281/N334-A339 (runbook §7.1).
  * disulfide C142-C219 bonded explicitly (CYX alone is not enough in tleap).
  * neutrality, atom count and per-species charge are verified across Amber -> GROMACS.
  * hydrogen mass repartitioning (3.024 amu, solute + lipids, not water) for dt = 4 fs.

Run in the mor-pilot env (parmed, tleap):
    python systems/build_topology.py --root systems --starts active inactive
"""
import argparse
import collections
import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import parmed

AA = {"ALA", "ARG", "ASN", "ASP", "ASH", "CYS", "CYX", "GLN", "GLU", "GLH", "GLY", "HID",
      "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP",
      "TYR", "VAL", "ACE", "NME"}
LIPRES = {"PA", "PC", "OL", "CHL"}
OFFSET = 67  # packmol/tleap residue 1 == ACE68
CANON = {166: "ASP", 167: "ARG", 168: "TYR", 281: "THR", 334: "ASN", 335: "PRO",
         336: "VAL", 337: "LEU", 338: "TYR", 339: "ALA"}
SS = (142, 219)
MOL2 = Path(__file__).resolve().parent.parent / "structures" / "params"


def xyz(l):
    return np.array([float(l[30:38]), float(l[38:46]), float(l[46:54])])


def parse(packed):
    """Group the packed PDB into molecules: receptor, dam, lipids (POPC = PA+PC+OL,
    CHL), waters, ions. Lipids/waters are delimited by TER records in packmol output."""
    rec, dam, mols = [], [], []
    cur = []
    for l in packed.read_text().splitlines():
        if l.startswith(("TER", "END")):
            if cur:
                mols.append(cur)
            cur = []
            continue
        if not l.startswith(("ATOM", "HETATM")):
            continue
        rn = l[17:21].strip()
        if rn in AA:
            rec.append(l)
        elif rn == "DAM":
            dam.append(l)
        else:
            cur.append(l)
    if cur:
        mols.append(cur)
    groups = collections.defaultdict(list)
    for m in mols:
        names = {l[17:21].strip() for l in m}
        if names <= {"WAT"}:
            # packmol writes all waters of a block without TER: split by residue number
            byres = collections.OrderedDict()
            for l in m:
                byres.setdefault(l[22:27], []).append(l)
            groups["WAT"] += list(byres.values())
        elif names <= {"Na+"} or names <= {"Cl-"}:
            for l in m:
                groups[l[17:21].strip()].append([l])
        elif names == {"PA", "PC", "OL"}:
            groups["POPC"].append(m)
        elif names == {"CHL"}:
            groups["CHL"].append(m)
        else:
            sys.exit(f"FATAL: unrecognised molecule with residues {names}")
    return rec, dam, groups


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", type=Path, default=Path("systems"))
    ap.add_argument("--starts", nargs="+", default=["active", "inactive"])
    ap.add_argument("--hmr", type=float, default=3.024)
    a = ap.parse_args()

    data = {}
    for s in a.starts:
        rec, dam, g = parse(a.root / s / "packed.pdb")
        # the crystal Na+ (inactive) is written by packmol as part of the protein file;
        # it surfaces here as a Na+ residue like any other and is kept (it is at the site).
        data[s] = (rec, dam, g)
        print(s, {k: len(v) for k, v in g.items()}, "receptor atoms", len(rec),
              "DAM atoms", len(dam))

    # ---------------- identical composition: common minimum, converting water->ion if short
    species = ["POPC", "CHL", "Na+", "Cl-"]
    target = {sp: min(len(data[s][2][sp]) for s in a.starts) for sp in species}
    # neutrality: keep Na-Cl difference identical (solute identical), use the max ion count
    for ion in ("Na+", "Cl-"):
        target[ion] = max(len(data[s][2][ion]) for s in a.starts)
    reports = {}
    for s in a.starts:
        rec, dam, g = data[s]
        rep = {"start": s, "removed": {}, "water_to_ion": {}}
        P = np.array([xyz(l) for l in rec])
        def dist_to_rec(m):
            c = np.array([xyz(l) for l in m]).mean(0)
            return float(np.min(np.linalg.norm(P - c, axis=1)))
        for sp in ("POPC", "CHL"):
            excess = len(g[sp]) - target[sp]
            if excess > 0:
                # remove the lipid(s) farthest from the receptor, from the fuller leaflet
                zs = np.array([np.array([xyz(l) for l in m]).mean(0)[2] for m in g[sp]])
                up = zs > 0
                fuller = up if up.sum() >= (~up).sum() else ~up
                cand = sorted((i for i in range(len(g[sp])) if fuller[i]),
                              key=lambda i: -dist_to_rec(g[sp][i]))[:excess]
                g[sp] = [m for i, m in enumerate(g[sp]) if i not in set(cand)]
                rep["removed"][sp] = excess
        for ion, name in (("Na+", "Na+"), ("Cl-", "Cl-")):
            short = target[ion] - len(g[ion])
            if short > 0:
                # convert bulk waters (|z| > 30 A, far from receptor) to ions
                cand = sorted((i for i, m in enumerate(g["WAT"]) if abs(xyz(m[0])[2]) > 30),
                              key=lambda i: -dist_to_rec(g["WAT"][i]))[:short]
                for i in cand:
                    o = g["WAT"][i][0]
                    g[ion].append([o[:12] + f"{name:<4s}{name:>4s}" + o[20:]])
                g["WAT"] = [m for i, m in enumerate(g["WAT"]) if i not in set(cand)]
                rep["water_to_ion"][ion] = short
        reports[s] = rep
    nwat = min(len(data[s][2]["WAT"]) for s in a.starts)
    for s in a.starts:
        g = data[s][2]
        excess = len(g["WAT"]) - nwat
        if excess:
            P = np.array([xyz(l) for l in data[s][0]])
            order = sorted(range(len(g["WAT"])),
                           key=lambda i: -float(np.min(np.linalg.norm(P - xyz(g["WAT"][i][0]),
                                                                      axis=1))))
            drop = set(order[:excess])
            g["WAT"] = [m for i, m in enumerate(g["WAT"]) if i not in drop]
            reports[s]["removed"]["WAT"] = excess
        reports[s]["composition"] = {k: len(v) for k, v in g.items()}
    comps = {s: reports[s]["composition"] for s in a.starts}
    if len({json.dumps(c, sort_keys=True) for c in comps.values()}) != 1:
        sys.exit(f"FATAL: compositions still differ: {comps}")
    print("composition (identical):", comps[a.starts[0]])

    # ---------------- per start: write pieces, tleap, verify, convert, HMR
    for s in a.starts:
        D = a.root / s
        rec, dam, g = data[s]
        rep = reports[s]
        box = None
        for l in (D / "packmol.inp").read_text().splitlines():
            if l.strip().startswith("pbc "):
                v = [float(x) for x in l.split()[1:7]]
                box = np.array(v[3:]) - np.array(v[:3])
        rep["box_A"] = box.round(2).tolist()

        (D / "tl_receptor.pdb").write_text(
            "\n".join(l for l in rec if (l[76:78].strip() or l[12:16].strip()[0]) != "H")
            + "\nTER\nEND\n")
        (D / "tl_dam.pdb").write_text("\n".join(
            l[:17] + "DAM L   1" + l[26:] for l in dam) + "\nTER\nEND\n")
        env = []
        for sp in ("POPC", "CHL", "WAT", "Na+", "Cl-"):
            for m in g[sp]:
                env += m + ["TER"]
        (D / "tl_env.pdb").write_text("\n".join(env) + "\nEND\n")

        i, j = SS[0] - OFFSET, SS[1] - OFFSET
        leap = f"""source leaprc.protein.ff19SB
source leaprc.lipid21
source leaprc.water.opc
source leaprc.gaff2
loadamberparams {MOL2}/damgo.frcmod
DAM = loadmol2 {MOL2}/damgo.mol2
rec = loadpdb tl_receptor.pdb
bond rec.{i}.SG rec.{j}.SG
lig = loadpdb tl_dam.pdb
env = loadpdb tl_env.pdb
sys = combine {{ rec lig env }}
set sys box {{ {box[0]:.3f} {box[1]:.3f} {box[2]:.3f} }}
charge sys
saveamberparm sys system.parm7 system.rst7
quit
"""
        (D / "build.leap").write_text(leap)
        r = subprocess.run(["tleap", "-f", "build.leap"], cwd=D, capture_output=True, text=True)
        (D / "tleap.log").write_text(r.stdout + r.stderr)
        errs = [l for l in (r.stdout + r.stderr).splitlines()
                if re.match(r"\s*(FATAL|Error!?)", l)]
        if errs or not (D / "system.parm7").exists():
            sys.exit(f"FATAL: tleap failed for {s}: {errs[:5]}")
        added = [l for l in r.stdout.splitlines() if "Added missing heavy atom" in l]
        rep["tleap_added_heavy_atoms"] = len(added)

        amb = parmed.load_file(str(D / "system.parm7"), xyz=str(D / "system.rst7"))
        q = sum(x.charge for x in amb.atoms)
        rep["net_charge"] = round(q, 4)
        if abs(q) > 1e-3:
            sys.exit(f"FATAL: {s} not neutral ({q:+.4f})")
        # contiguity: receptor then DAM first
        names = [res.name for res in amb.residues]
        ndam = names.index("DAM")
        if any(n not in AA for n in names[:ndam]) or names.count("DAM") != 1:
            sys.exit(f"FATAL: {s}: receptor/DAM not contiguous at start of topology")
        n_solute = amb.residues[ndam].atoms[-1].idx + 1
        rep["solute_atoms_1_to"] = n_solute
        # disulfide present
        sg = [at for res in amb.residues[:ndam] if res.name == "CYX" for at in res.atoms
              if at.name == "SG"]
        bonded = any(b.atom1 in sg and b.atom2 in sg for b in amb.bonds)
        if len(sg) != 2 or not bonded:
            sys.exit(f"FATAL: {s}: disulfide C142-C219 not bonded")

        # canonical numbering check on the built topology
        for k, res in enumerate(amb.residues[:ndam]):
            canon = k + 1 + OFFSET
            if canon in CANON and res.name.replace("ASH", "ASP") != CANON[canon]:
                sys.exit(f"FATAL: {s}: topology residue {k+1} ({res.name}) is not canonical "
                         f"{CANON[canon]}{canon}")

        if a.hmr:
            parmed.tools.HMassRepartition(amb, a.hmr).execute()
        aq = collections.defaultdict(float); an = collections.Counter()
        for res in amb.residues:
            aq[res.name] += sum(x.charge for x in res.atoms); an[res.name] += 1
        amb.save(str(D / "system.top"), format="gromacs", overwrite=True)
        amb.save(str(D / "system.gro"), overwrite=True)
        gm = parmed.load_file(str(D / "system.top"), xyz=str(D / "system.gro"))
        gq = collections.defaultdict(float); gn = collections.Counter()
        for res in gm.residues:
            gq[res.name] += sum(x.charge for x in res.atoms); gn[res.name] += 1
        bad = [nm for nm in set(aq) | set(gq)
               if an[nm] != gn[nm] or abs(aq[nm] - gq[nm]) > 1e-3]
        if len(amb.atoms) != len(gm.atoms) or bad:
            sys.exit(f"FATAL: {s}: Amber->GROMACS conversion changed the system: {bad}")
        rep["atoms"] = len(gm.atoms)
        hmass = sorted({round(at.mass, 3) for at in gm.atoms if at.element == 1
                        and at.residue.name not in ("WAT",)})
        rep["solute_lipid_H_masses"] = hmass

        # system_ref.pdb: solute only, canonical numbering, serial == topology index
        out = []
        for at in gm.atoms[:n_solute]:
            res = at.residue
            k = res.idx
            resnum = k + 1 + OFFSET if res.name != "DAM" else 1
            chain = "R" if res.name != "DAM" else "L"
            nm = at.name if len(at.name) == 4 else f" {at.name:<3}"
            out.append(f"ATOM  {at.idx+1:5d} {nm:4s} {res.name:>3s} {chain}{resnum:4d}    "
                       f"{at.xx:8.3f}{at.xy:8.3f}{at.xz:8.3f}  1.00  0.00")
        (D / "system_ref.pdb").write_text("\n".join(out) + "\nEND\n")

        # index groups: thermostat (Solute / MEMB / SOLV) and convenience groups
        solute = list(range(1, n_solute + 1))
        memb = [at.idx + 1 for at in gm.atoms if at.residue.name in LIPRES]
        solv = [at.idx + 1 for at in gm.atoms
                if at.residue.name in ("WAT", "Na+", "Cl-")]
        if len(solute) + len(memb) + len(solv) != len(gm.atoms):
            sys.exit(f"FATAL: {s}: index groups do not partition the system")
        lig = [at.idx + 1 for at in gm.atoms if at.residue.name == "DAM"]
        heavy_solute = [at.idx + 1 for at in gm.atoms[:n_solute] if at.element != 1]
        rec_heavy = [i for i in heavy_solute if i < lig[0]]
        with open(D / "index.ndx", "w") as f:
            for name, idx in (("System", list(range(1, len(gm.atoms) + 1))),
                              ("Solute", solute), ("MEMB", memb), ("SOLV", solv),
                              ("DAM", lig), ("Solute_heavy", heavy_solute),
                              ("Receptor_heavy", rec_heavy)):
                f.write(f"[ {name} ]\n")
                for c in range(0, len(idx), 15):
                    f.write(" ".join(f"{x:6d}" for x in idx[c:c + 15]) + "\n")
        rep["groups"] = {"Solute": len(solute), "MEMB": len(memb), "SOLV": len(solv),
                         "DAM": f"{lig[0]}-{lig[-1]}"}
        # 1-based indices consumed by systems/render_plumed.py
        rep["plumed"] = {
            "receptor_last": lig[0] - 1,
            "dam_first": lig[0], "dam_last": lig[-1],
            "dam_heavy": [at.idx + 1 for at in gm.atoms if at.residue.name == "DAM"
                          and at.element != 1],
            "dam_N1": [at.idx + 1 for at in gm.atoms
                       if at.residue.name == "DAM" and at.name == "N1"][0],
            "na": [at.idx + 1 for at in gm.atoms if at.residue.name == "Na+"],
        }

        # position restraints on solute heavy atoms, switched by -DPOSRES -DPOSRES_FC=...
        top = (D / "system.top").read_text()
        blocks = re.split(r"(?=\[ moleculetype \])", top)
        new, done = [], []
        for b in blocks:
            mname = re.search(r"\[ moleculetype \]\s*\n(?:;.*\n)*\s*(\S+)", b)
            if mname:
                atoms_sec = b.split("[ atoms ]")[1].split("[")[0]
                rows = [l.split() for l in atoms_sec.splitlines()
                        if l.strip() and not l.strip().startswith(";")]
                resn = {r[3] for r in rows}
                if resn & (AA | {"DAM"}) and not resn & {"WAT", "PA", "CHL"}:
                    itp = f"posre_{mname.group(1)}.itp"
                    with open(D / itp, "w") as f:
                        f.write("[ position_restraints ]\n")
                        for r in rows:
                            if not r[4].startswith("H"):
                                f.write(f"{int(r[0]):6d}     1  POSRES_FC  POSRES_FC  POSRES_FC\n")
                    b = b.rstrip() + f'\n\n#ifdef POSRES\n#include "{itp}"\n#endif\n\n'
                    done.append(mname.group(1))
            new.append(b)
        (D / "system.top").write_text("".join(new))
        rep["posres_moleculetypes"] = done
        if len(done) != 2:
            sys.exit(f"FATAL: {s}: expected posres for receptor and DAM, got {done}")
        (D / "build_report.json").write_text(json.dumps(rep, indent=1))
        print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
