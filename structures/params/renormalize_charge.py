"""Make the DAMGO mol2 partial charges sum to exactly the formal charge (+1).

The AM1-BCC charges inherited from the parent repo sum to +1.0020 (rounding in the
6-decimal mol2). That residual would leave every system 0.002 e non-neutral. The
residual is spread evenly over all atoms (shift ~2.7e-5 e/atom); the original is kept
as damgo.orig.mol2.  Usage: python renormalize_charge.py damgo.mol2 1"""
import shutil
import sys
from pathlib import Path

p, target = Path(sys.argv[1]), int(sys.argv[2])
orig = p.with_suffix(".orig.mol2")
if not orig.exists():
    shutil.copy(p, orig)
txt = orig.read_text()
head, rest = txt.split("@<TRIPOS>ATOM\n")
atoms, tail = rest.split("@<TRIPOS>BOND")
rows = atoms.rstrip("\n").splitlines()
q = [float(r.split()[8]) for r in rows]
shift = (target - sum(q)) / len(q)
newq = [round(x + shift, 6) for x in q]
newq[-1] = round(newq[-1] + (target - sum(newq)), 6)
out = []
for r, x in zip(rows, newq):
    f = r.rsplit(None, 1)
    out.append(f"{f[0]} {x:10.6f}")
p.write_text(head + "@<TRIPOS>ATOM\n" + "\n".join(out) + "\n@<TRIPOS>BOND" + tail)
print(f"sum {sum(q):+.6f} -> {sum(newq):+.6f}; shift {shift:+.2e} e/atom")
