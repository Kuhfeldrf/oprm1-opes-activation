"""Daily OPES monitor (runbook §9): recrossings, bias growth, walls, ligand retention.

Usage: python analysis/monitor_opes.py runs/prod1 [--boxes results/endpoints.json]
Prints, per start: simulated time per walker, CV ranges, opes.rct / nker / neff trend,
fraction of frames touching a wall, ligand retention extremes, and committed basin
crossings (both CVs inside the measured basin boxes) in each direction.
"""
import argparse
import glob
import json
from pathlib import Path

import numpy as np


def load(p):
    with open(p) as f:
        cols = f.readline().split()[2:]
    d = np.loadtxt(p, comments="#", ndmin=2)
    return {c: d[:, i] for i, c in enumerate(cols)}


def crossings(x, y, bi, ba):
    def inside(b):
        return (x >= b[0]) & (x <= b[1]) & (y >= b[2]) & (y <= b[3])
    st, i2a, a2i = None, 0, 0
    for i, a in zip(inside(bi), inside(ba)):
        s = "I" if i else "A" if a else None
        if s and st and s != st:
            i2a += s == "A"; a2i += s == "I"
        if s:
            st = s
    return i2a, a2i


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--boxes", default="results/endpoints.json")
    a = ap.parse_args()
    boxes = json.load(open(a.boxes))["opes_proposal"]["basin_boxes_for_crossing_count"]
    out = {}
    for start in ("active", "inactive"):
        files = sorted(glob.glob(f"{a.run}/{start}/w*/COLVAR.*"))
        if not files:
            continue
        rep = {"walkers": []}
        tot = [0, 0]
        for f in files:
            c = load(f)
            # drop restart duplicates (PLUMED appends on restart)
            _, idx = np.unique(c["time"], return_index=True)
            c = {k: v[np.sort(idx)] for k, v in c.items()}
            i2a, a2i = crossings(c["d_tm36"], c["rmsd_npxxy"], boxes["inactive"], boxes["active"])
            tot[0] += i2a; tot[1] += a2i
            wall = (c["uw.bias"] > 0.1) | (c["lw.bias"] > 0.1)
            rep["walkers"].append({
                "file": Path(f).name, "ns": round(float(c["time"][-1]) / 1000, 2),
                "d_tm36_range": [round(float(c["d_tm36"].min()), 3), round(float(c["d_tm36"].max()), 3)],
                "rmsd_range": [round(float(c["rmsd_npxxy"].min()), 3), round(float(c["rmsd_npxxy"].max()), 3)],
                "frac_at_wall": round(float(wall.mean()), 4),
                "lig_cont_min": round(float(c["lig_cont"].min()), 1),
                "d_salt_max_nm": round(float(c["d_salt"].max()), 3),
                "na_site_last": round(float(c["na_site"][-1]), 2),
                "i2a": i2a, "a2i": a2i,
            })
        c0 = load(files[0])
        n = len(c0["time"])
        q = [max(0, int(n * f) - 1) for f in (0.25, 0.5, 0.75, 1.0)]
        rep["rct_kJ_at_quarters"] = [round(float(c0["opes.rct"][i]), 2) for i in q]
        rep["nker_at_quarters"] = [int(c0["opes.nker"][i]) for i in q]
        rep["neff_last"] = round(float(c0["opes.neff"][-1]), 1)
        rep["aggregate_ns"] = round(sum(w["ns"] for w in rep["walkers"]), 1)
        rep["crossings_total_i2a_a2i"] = tot
        out[start] = rep
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
