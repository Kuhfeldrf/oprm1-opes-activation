"""Nb39-receptor interface retention over the G3 unbiased leg.

The equilibration leg logged the plain PLUMED template, which has no nb_cont, so this
measures it from the trajectory with the same definition: receptor heavy atoms
(Receptor_heavy) vs Nb39 heavy atoms (NB39_heavy), 4.5 A cutoff, minimum image.
Usage: python nb_contacts.py SYSTEM_DIR
"""
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import capped_distance

D = Path(sys.argv[1])


def groups(ndx):
    out, cur = {}, None
    for line in ndx.read_text().splitlines():
        line = line.strip()
        if line.startswith("["):
            cur = line.strip("[] ")
            out[cur] = []
        elif line and cur:
            out[cur] += [int(x) - 1 for x in line.split()]
    return out


g = groups(D / "index.ndx")
u = mda.Universe(str(D / "prod.tpr"), str(D / "prod.xtc"))
rec, nb = u.atoms[g["Receptor_heavy"]], u.atoms[g["NB39_heavy"]]
rows = []
for ts in u.trajectory[::10]:
    pairs = capped_distance(nb.positions, rec.positions, 4.5, box=ts.dimensions,
                            return_distances=False)
    rows.append((ts.time / 1000.0, len(pairs), len(np.unique(pairs[:, 0])),
                 len(np.unique(rec[np.unique(pairs[:, 1])].resids))))
a = np.array(rows)
first, last = a[a[:, 0] <= 1.0], a[a[:, 0] >= 9.0]
print(f"{D.name}: frames {len(a)}  "
      f"pairs 0-1 ns {first[:, 1].mean():.0f} -> 9-10 ns {last[:, 1].mean():.0f} "
      f"(retention {last[:, 1].mean() / first[:, 1].mean():.2f}); "
      f"Nb39 atoms in contact {first[:, 2].mean():.0f} -> {last[:, 2].mean():.0f}; "
      f"receptor residues {first[:, 3].mean():.0f} -> {last[:, 3].mean():.0f}; "
      f"min pairs {a[:, 1].min():.0f}")
