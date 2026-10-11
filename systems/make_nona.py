"""Derive a sodium-free copy of a built system: the Na+ at the D2.50 (Asp116) site is
swapped with a bulk water, so topology, atom order, composition and every PLUMED index
stay identical; only coordinates change.

Used for the v0.2 follow-up (NOTES.md, Day 5): in prod2 the inactive start kept its
D2.50 Na+ for the whole 1.5 us and the active start never had one, so the two starts
differed in a variable the bias never touched. inactive_nona removes that difference.

Usage: python systems/make_nona.py systems/inactive systems/inactive_nona [frame.gro]
  frame.gro defaults to <src>/prod.gro (end of the 100 ns unbiased leg).
Writes <dst>/system.gro (no velocities) and <dst>/build_report.json (with a provenance
block); topology, index, posres and system_ref.pdb are copied unchanged.
"""
import json
import shutil
import sys
from pathlib import Path

import numpy as np

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
frame = Path(sys.argv[3]) if len(sys.argv) > 3 else src / "prod.gro"
dst.mkdir(parents=True, exist_ok=True)
rep = json.loads((src / "build_report.json").read_text())

lines = frame.read_text().splitlines()
title, n, atoms, box = lines[0], int(lines[1]), lines[2:2 + int(lines[1])], lines[2 + int(lines[1])]
resn = [a[5:10].strip() for a in atoms]
name = [a[10:15].strip() for a in atoms]
xyz = np.array([[float(a[20:28]), float(a[28:36]), float(a[36:44])] for a in atoms])
L = np.array([float(v) for v in box.split()[:3]])

# Asp116 CG: canonical residue 116 in system_ref.pdb (serial == topology index)
cg = None
for l in (src / "system_ref.pdb").read_text().splitlines():
    if l.startswith(("ATOM", "HETATM")) and l[12:16].strip() == "CG" and int(l[22:26]) == 116:
        cg = int(l[6:11]) - 1
assert cg is not None and resn[cg] in ("ASP", "ASH"), "Asp116 CG not found"


def mind(a, b):
    d = b - a
    d -= L * np.round(d / L)
    return np.linalg.norm(d, axis=-1)


na = [i - 1 for i in rep["plumed"]["na"]]
dna = mind(xyz[cg], xyz[na])
site = int(na[int(np.argmin(dna))])
assert dna.min() < 0.45, f"no Na+ within 0.45 nm of D116 CG (closest {dna.min():.2f} nm)"

# bulk water: farthest from the bilayer midplane in z (PBC), > 1.5 nm from any solute
# or lipid atom
P = np.array([i for i, (r, nm) in enumerate(zip(resn, name)) if r == "PC" and nm == "P31"])
mid = xyz[P, 2].mean()
nonsolv = np.array([i for i, r in enumerate(resn) if r not in ("WAT", "Na+", "Cl-")])
ow = np.array([i for i, (r, nm) in enumerate(zip(resn, name)) if r == "WAT" and nm == "O"])
dz = np.abs((xyz[ow, 2] - mid + L[2] / 2) % L[2] - L[2] / 2)
best = None
for i in ow[np.argsort(-dz)][:200]:
    if mind(xyz[i], xyz[nonsolv]).min() > 1.5 and mind(xyz[i], xyz[na]).min() > 0.8:
        best = int(i)
        break
assert best is not None, "no bulk water found"
# the water molecule = consecutive atoms of the same residue starting at its O
w = [best]
while best + len(w) < n and resn[best + len(w)] == "WAT" and \
        atoms[best + len(w)][:5] == atoms[best][:5]:
    w.append(best + len(w))

na_old, o_old = xyz[site].copy(), xyz[best].copy()
xyz[w] += na_old - o_old          # water moves into the site (keeps its geometry)
xyz[site] = o_old                 # Na+ goes to bulk

out = [f"{title} | {Path(sys.argv[0]).name}: D2.50 Na+ (atom {site+1}) swapped with "
       f"water O {best+1}", f"{n:5d}"]
for a, (x, y, z) in zip(atoms, xyz):
    out.append(f"{a[:20]}{x:8.3f}{y:8.3f}{z:8.3f}")
out.append(box)
(dst / "system.gro").write_text("\n".join(out) + "\n")

for f in ("system.top", "index.ndx", "system_ref.pdb") + tuple(p.name for p in src.glob("posre_*.itp")):
    shutil.copy2(src / f, dst / f)
rep["derived"] = {
    "from": str(frame), "operation": "swap D2.50 Na+ with bulk water",
    "na_atom": int(site) + 1, "na_to_d116cg_before_nm": round(float(dna.min()), 3),
    "water_atoms": [int(i) + 1 for i in w],
    "na_to_bilayer_mid_after_nm": round(float(abs((o_old[2] - mid + L[2] / 2) % L[2] - L[2] / 2)), 2),
    "na_min_dist_to_nonsolvent_after_nm": round(float(mind(o_old, xyz[nonsolv]).min()), 2),
}
(dst / "build_report.json").write_text(json.dumps(rep, indent=1))
print(json.dumps(rep["derived"], indent=1))
