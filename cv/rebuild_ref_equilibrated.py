"""Replace the CV2 reference with the EQUILIBRATED inactive structure (runbook §6.2, §7.2).

Coordinates: mean Cα positions of the 140 reference atoms (134 scaffold + 6 NPxxYA) over
the second half of the inactive unbiased leg, after superposing each frame on the
scaffold. The atom serials, the occupancy/B-factor split, and therefore what CV2
means, are unchanged; only the inactive geometry it is measured against moves from the
cleaned 9PXU crystal coordinates to the simulated inactive basin.

Usage: python cv/rebuild_ref_equilibrated.py systems/inactive cv/ref_inactive_npxxy.pdb
       (the crystal-derived file is kept as cv/ref_inactive_npxxy.crystal.pdb)
"""
import shutil
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.analysis import align

D, ref_path = Path(sys.argv[1]), Path(sys.argv[2])
crystal = ref_path.with_name(ref_path.stem + ".crystal.pdb")
if not crystal.exists():
    shutil.copy(ref_path, crystal)
lines = [l for l in crystal.read_text().splitlines() if l.startswith("ATOM")]
serial = [int(l[6:11]) for l in lines]
occ = np.array([float(l[54:60]) for l in lines])
scaf_idx = [s - 1 for s, o in zip(serial, occ) if o > 0]

u = mda.Universe(str(D / "prod.tpr"), str(D / "prod.xtc"))
sel = u.atoms[[s - 1 for s in serial]]
scaf = u.atoms[scaf_idx]
nfr = len(u.trajectory)
frames = range(nfr // 2, nfr)
ref_scaf = None
acc = np.zeros((len(sel), 3))
for i in frames:
    u.trajectory[i]
    X = scaf.positions.copy()
    if ref_scaf is None:
        ref_scaf = X
    R, _ = align.rotation_matrix(X - X.mean(0), ref_scaf - ref_scaf.mean(0))
    acc += (sel.positions - X.mean(0)) @ R.T + ref_scaf.mean(0)
mean = acc / len(frames)
out = [lines[0][:0] + "REMARK dual-weight RMSD reference. occupancy is align weight (TM1-5 scaffold CA), "
       "beta is displace weight (NPxxYA 334-339 CA)",
       f"REMARK EQUILIBRATED inactive: scaffold-aligned mean over frames {frames.start}-{nfr - 1} "
       f"of {D.name}/prod.xtc; serials are topology indices"]
for l, x in zip(lines, mean):
    out.append(l[:30] + f"{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}" + l[54:])
ref_path.write_text("\n".join(out) + "\nEND\n")
print(f"wrote {ref_path} from {len(frames)} frames; crystal-derived kept at {crystal}")
