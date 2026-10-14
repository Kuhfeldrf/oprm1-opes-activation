#!/usr/bin/env bash
#SBATCH --job-name=oprm1_nlx_params
#SBATCH --partition=short
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --output=logs/nlx_params_%j.log
# GAFF2/AM1-BCC parameters for naloxone (residue NLX, net +1: tertiary amine protonated).
# Heavy atoms come from 9PXU (CCD A1APV), so stereochemistry is the crystal's; hydrogens
# are added by Open Babel at pH 7.4 and the formula/charge are asserted before sqm runs.
# Usage: sbatch structures/params/make_naloxone.sh
set -euo pipefail
ROOT=/scratch/kuhfeldr-Kuhfeld_temp/oprm1-opes-activation
source /scratch/kuhfeldr-Kuhfeld_temp/miniforge3/etc/profile.d/conda.sh
conda activate mor-pilot
P=$ROOT/structures/params; W=$P/nlx_build; mkdir -p $W; cd $W
python - <<'PY'
import gemmi
st = gemmi.read_structure("../../raw/9PXU.cif"); st.remove_alternative_conformations()
res = [r for r in st[0]["R"] if r.name == "A1APV"]
assert len(res) == 1, res
lines = []
for i, a in enumerate(res[0], 1):
    if a.element.name == "H":
        continue
    p = a.pos
    lines.append(f"HETATM{i:5d} {a.name:<4s} NLX L   1    {p.x:8.3f}{p.y:8.3f}{p.z:8.3f}  1.00  0.00"
                 f"          {a.element.name:>2s}")
open("nlx_heavy.pdb", "w").write("\n".join(lines) + "\nEND\n")
print("heavy atoms", len(lines))
PY
obabel nlx_heavy.pdb -O nlx_h.mol2 -p 7.4 2>&1 | tail -1
python - <<'PY'
from collections import Counter
t = open("nlx_h.mol2").read().split("@<TRIPOS>ATOM")[1].split("@<TRIPOS>")[0].split("\n")
el = Counter("".join(c for c in l.split()[5].split(".")[0] if c.isalpha()) for l in t if l.split())
q = sum(float(l.split()[8]) for l in t if l.split())
print("formula", dict(el), "obabel charge sum %.2f" % q)
assert el == Counter({"C": 19, "H": 22, "N": 1, "O": 4}), f"unexpected formula {el} (want C19H22NO4+)"
PY
antechamber -i nlx_h.mol2 -fi mol2 -o nlx_ac.mol2 -fo mol2 -c bcc -nc 1 -at gaff2 -rn NLX -s 2 -pf y
parmchk2 -i nlx_ac.mol2 -f mol2 -o ../naloxone.frcmod -s gaff2
cp nlx_ac.mol2 ../naloxone.mol2
python ../renormalize_charge.py ../naloxone.mol2 1
grep -c "ATTN" ../naloxone.frcmod || true
echo DONE
