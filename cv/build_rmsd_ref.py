"""Runbook §7.2: build the dual-weight PLUMED RMSD reference for CV2.

Input is a PDB written FROM THE SIMULATION SYSTEM (gmx editconf -f topol.tpr),
so atom serials equal PLUMED's 1-based topology indices. Output contains only
CA atoms:
    TM1-TM5 scaffold CA  occupancy 1.00 (align)    B-factor 0.00
    NPxxYA 334-339 CA    occupancy 0.00            B-factor 1.00 (displace)
Coordinates are taken from the input structure (the inactive reference).

Usage: python cv/build_rmsd_ref.py system_ref.pdb out.pdb [--offset N] [--chain X]
--offset: system_resnum = canonical_resnum + offset (0 when numbering is canonical).
"""
import argparse
import sys
from pathlib import Path
import gemmi
sys.path.insert(0, str(Path(__file__).parent))
from scaffold import SCAFFOLD, NPXXYA

EXPECT = {334: "ASN", 335: "PRO", 336: "VAL", 337: "LEU", 338: "TYR", 339: "ALA"}

ap = argparse.ArgumentParser()
ap.add_argument("inp"); ap.add_argument("out")
ap.add_argument("--offset", type=int, default=0)
ap.add_argument("--chain", default=None)
a = ap.parse_args()

st = gemmi.read_structure(a.inp)
ch = st[0][a.chain] if a.chain else st[0][0]
ca = {r.seqid.num: (r.name, r["CA"][0]) for r in ch if r.find_atom("CA", "*")}
for n, nm in EXPECT.items():
    got = ca.get(n + a.offset, ("missing",))[0]
    if got != nm:
        sys.exit(f"FATAL: canonical {n} -> system {n + a.offset} is {got}, expected {nm}")
missing = [n for n in SCAFFOLD if n + a.offset not in ca]
if missing:
    sys.exit(f"FATAL: scaffold residues missing from system: {missing}")

rows = []
for n in SCAFFOLD:
    rows.append((n, 1.0, 0.0))
for n in NPXXYA:
    rows.append((n, 0.0, 1.0))
rows.sort(key=lambda t: ca[t[0] + a.offset][1].serial)
with open(a.out, "w") as f:
    # No '=' in REMARKs: PLUMED parses REMARK key=value pairs as arguments.
    f.write("REMARK dual-weight RMSD reference. occupancy is align weight (TM1-5 scaffold CA), "
            "beta is displace weight (NPxxYA 334-339 CA)\n")
    f.write(f"REMARK source {Path(a.inp).name} offset {a.offset}; serials are topology indices\n")
    for n, occ, b in rows:
        name, at = ca[n + a.offset]
        p = at.pos
        f.write("ATOM  %5d  CA  %3s %1s%4d    %8.3f%8.3f%8.3f%6.2f%6.2f           C\n"
                % (at.serial, name, ch.name[:1] or "A", n + a.offset, p.x, p.y, p.z, occ, b))
    f.write("END\n")
print(f"wrote {a.out}: {len(SCAFFOLD)} align atoms, {len(NPXXYA)} displace atoms")
