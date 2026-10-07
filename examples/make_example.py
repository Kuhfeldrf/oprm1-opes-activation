"""Build examples/two-cv-smoke/ from the production system (runbook §12.3).

The example is the receptor + DAMGO (4,634 atoms) cut out of the built active system,
with the same force-field parameters, in VACUUM with a large box. That is not physics,
and the example does not claim to be: it shows, in seconds on a CPU, that both CVs
resolve to the right atoms, respond to motion, and carry a bias. The production
membrane system is far too large to commit.

The topology is extracted as TEXT from system.top (header + the receptor and DAMGO
moleculetypes, posres includes dropped), so the parameters are byte-identical to
production. Coordinates are the first N atoms of the minimised system.

Usage: python examples/make_example.py systems/active examples/two-cv-smoke
"""
import json
import re
import shutil
import sys
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
dst.mkdir(parents=True, exist_ok=True)
rep = json.loads((src / "build_report.json").read_text())
n = rep["solute_atoms_1_to"]
keep = rep["posres_moleculetypes"]            # receptor ("system1") and "DAM"

top = (src / "system.top").read_text()
blocks = re.split(r"(?=\[ moleculetype \])", top)
header, mols = blocks[0], blocks[1:]
out = [header]
for b in mols:
    name = re.search(r"\[ moleculetype \]\s*\n(?:;.*\n)*\s*(\S+)", b).group(1)
    if name in keep:
        b = re.sub(r"#ifdef POSRES\n#include \"[^\"]+\"\n#endif\n", "", b)
        b = b.split("[ system ]")[0]
        out.append(b)
out.append("[ system ]\nDAMGO-muOR receptor+ligand, vacuum smoke example\n\n"
           "[ molecules ]\n; Compound       #mols\n" + "".join(f"{m:<16s} 1\n" for m in keep))
(dst / "solute.top").write_text("".join(out))

coords = src / "min.gro" if (src / "min.gro").exists() else src / "system.gro"
g = coords.read_text().splitlines()
atoms = g[2:2 + n]
(dst / "solute.gro").write_text(f"muOR+DAMGO from {coords.name}\n{n:5d}\n" + "\n".join(atoms)
                                + "\n  15.00000  15.00000  15.00000\n")

shutil.copy(src / "system_ref.pdb", dst / "system_ref.pdb")
shutil.copy(src.parent.parent / "cv" / "ref_inactive_npxxy.pdb", dst / "ref_inactive_npxxy.pdb")
print(f"example: {n} atoms from {coords.name}; molecules {keep}; files in {dst}")
