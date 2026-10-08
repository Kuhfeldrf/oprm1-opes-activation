"""Fill cv/plumed/*.dat templates for one system from its build_report.json.

Usage: python systems/render_plumed.py <system_dir> <template> <out> [KEY=VALUE ...]
Indices are never typed by hand: they come from the topology the run uses."""
import json
import re
import sys
from pathlib import Path

D, tmpl, out = Path(sys.argv[1]).resolve(), Path(sys.argv[2]), Path(sys.argv[3])
rep = json.loads((D / "build_report.json").read_text())["plumed"]
root = D.parent.parent
subs = {
    "SYSTEM_REF": str(D / "system_ref.pdb"),
    "RMSD_REF": str(root / "cv" / "ref_inactive_npxxy.pdb"),
    "RECEPTOR_LAST": rep["receptor_last"],
    "DAM_FIRST": rep["dam_first"], "DAM_LAST": rep["dam_last"],
    "DAM_N1": rep["dam_N1"],
    "DAM_HEAVY": ",".join(map(str, rep["dam_heavy"])),
    "NA_LIST": ",".join(map(str, rep["na"])),
    "RECEPTOR_HEAVY_GROUP": "rec_heavy",
    "STRIDE": 2500, "SMOKE_D_AT": 1.00, "WALL_KAPPA": "50000.0",
}
for kv in sys.argv[4:]:
    k, v = kv.split("=", 1)
    subs[k] = v
txt = tmpl.read_text()
if "@RECEPTOR_HEAVY_GROUP@" in txt:
    # receptor heavy atoms as a PLUMED GROUP read from the index file
    txt = txt.replace("MOLINFO", f"rec_heavy: GROUP NDX_FILE={D/'index.ndx'} NDX_GROUP=Receptor_heavy\nMOLINFO", 1)
for k, v in subs.items():
    txt = txt.replace(f"@{k}@", str(v))
left = sorted(set(re.findall(r"@[A-Z][A-Z_]*@", txt)))
if left:
    sys.exit(f"FATAL: unfilled placeholders {left}")
out.write_text(txt)
print(f"wrote {out}")
