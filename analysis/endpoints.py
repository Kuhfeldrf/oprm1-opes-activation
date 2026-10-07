"""Endpoint statistics from the unbiased legs (runbook §6.4, §8) and OPES parameters.

Reads COLVAR files written by plumed_measure.dat (columns d_tm36, rmsd_npxxy, d_hbond,
chi_w295, d_salt, na_site, lig_cont; time in ps) for the active and inactive legs.

Reports per endpoint: mean/SD of both CVs after equilibration (--skip ns), the §8
failure signatures (CV1 drift > 0.3 nm between first and last 10 ns; CV2 climbing
monotonically, tested as a significant linear trend over the last half), and the
secondary observables. Then proposes OPES parameters:
    SIGMA  = within-well SD, the smaller of the two endpoints per CV (OPES convention:
             the kernel width should not exceed the fluctuation in the narrowest well)
    walls  = 0.10 nm below the inactive mean (CV1) / 0 (CV2), and 0.30 nm / 0.15 nm
             above the active mean; tied to measured values, not copied.

Usage: python analysis/endpoints.py --active A/COLVAR --inactive I/COLVAR [--skip 10]
"""
import argparse
import json

import numpy as np


def load(p):
    with open(p) as f:
        cols = f.readline().split()[2:]
    d = np.loadtxt(p, comments="#")
    return {c: d[:, i] for i, c in enumerate(cols)}


def stats(cv, skip_ps):
    t = cv["time"]
    m = t >= skip_ps
    out = {"t_end_ns": float(t[-1] / 1000), "n": int(m.sum())}
    for c in ("d_tm36", "rmsd_npxxy", "d_hbond", "d_salt", "na_site", "lig_cont"):
        if c in cv:
            out[c] = {"mean": float(cv[c][m].mean()), "sd": float(cv[c][m].std()),
                      "min": float(cv[c][m].min()), "max": float(cv[c][m].max())}
    first = cv["d_tm36"][t < t[0] + 10000].mean()
    last = cv["d_tm36"][t > t[-1] - 10000].mean()
    out["cv1_drift_first_vs_last_10ns_nm"] = float(last - first)
    h = t > t[-1] / 2
    slope = np.polyfit(t[h] / 1000, cv["rmsd_npxxy"][h], 1)[0]
    out["cv2_trend_last_half_nm_per_100ns"] = float(slope * 100)
    # W295 chi1: circular mean
    if "chi_w295" in cv:
        a = cv["chi_w295"][m]
        out["chi_w295_circmean_deg"] = float(np.degrees(np.angle(np.exp(1j * a).mean())))
    flags = []
    if abs(out["cv1_drift_first_vs_last_10ns_nm"]) > 0.3:
        flags.append("CV1 drift > 3 A (runbook §8 failure signature)")
    if abs(out["cv2_trend_last_half_nm_per_100ns"]) > 0.05:
        flags.append("CV2 still trending in the second half (not settled)")
    if "lig_cont" in cv and cv["lig_cont"][m].min() < 0.5 * cv["lig_cont"][m].mean():
        flags.append("ligand contacts dropped below half their mean (check for unbinding)")
    out["flags"] = flags
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--active", required=True)
    ap.add_argument("--inactive", required=True)
    ap.add_argument("--skip", type=float, default=10.0, help="ns discarded as relaxation")
    ap.add_argument("--out", default=None)
    a = ap.parse_args()
    res = {s: stats(load(p), a.skip * 1000) for s, p in
           (("active", a.active), ("inactive", a.inactive))}
    A, I = res["active"], res["inactive"]
    sep1 = A["d_tm36"]["mean"] - I["d_tm36"]["mean"]
    sep2 = A["rmsd_npxxy"]["mean"] - I["rmsd_npxxy"]["mean"]
    res["separation"] = {
        "cv1_nm": sep1, "cv2_nm": sep2,
        "cv1_in_sd": sep1 / max(A["d_tm36"]["sd"], I["d_tm36"]["sd"]),
        "cv2_in_sd": sep2 / max(A["rmsd_npxxy"]["sd"], I["rmsd_npxxy"]["sd"]),
        "cv1_gt_0.2nm": sep1 > 0.2,
    }
    s1 = min(A["d_tm36"]["sd"], I["d_tm36"]["sd"])
    s2 = min(A["rmsd_npxxy"]["sd"], I["rmsd_npxxy"]["sd"])
    res["opes_proposal"] = {
        "SIGMA": f"{s1:.3f},{s2:.3f}",
        "LW": f"{I['d_tm36']['mean'] - 0.10:.2f},0.0",
        "UW": f"{A['d_tm36']['mean'] + 0.30:.2f},{A['rmsd_npxxy']['mean'] + 0.15:.2f}",
        "BARRIER_kJ": 50,
        "basin_boxes_for_crossing_count": {
            "inactive": [I["d_tm36"]["mean"] - 3 * I["d_tm36"]["sd"],
                         I["d_tm36"]["mean"] + 2 * I["d_tm36"]["sd"],
                         0.0, I["rmsd_npxxy"]["mean"] + 2 * I["rmsd_npxxy"]["sd"]],
            "active": [A["d_tm36"]["mean"] - 2 * A["d_tm36"]["sd"],
                       A["d_tm36"]["mean"] + 3 * A["d_tm36"]["sd"],
                       A["rmsd_npxxy"]["mean"] - 2 * A["rmsd_npxxy"]["sd"],
                       A["rmsd_npxxy"]["mean"] + 3 * A["rmsd_npxxy"]["sd"]]},
    }
    txt = json.dumps(res, indent=1)
    print(txt)
    if a.out:
        open(a.out, "w").write(txt)


if __name__ == "__main__":
    main()
