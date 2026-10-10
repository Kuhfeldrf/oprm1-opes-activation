"""§10 criteria 3–4 and the §13 time-series / convergence figures for one OPES run.

- Walker-to-walker error: each walker reweighted alone with the shared bias gives an
  independent FES estimate; the error map is their SD (after removing each walker's
  constant). Reported for the first and second half separately — §10.4 requires that it
  does not grow.
- Block stability: FES from quarter 3 vs quarter 4 (second half), RMS difference.
- Path overlay (§10.3): the 2D regions visited from each start, with basin boxes.
- CV time series per walker, basin boxes shaded, committed crossings marked.
- Strided COLVAR copies for the repository (§12.2).

Usage: python analysis/convergence.py runs/prod2 --out results/fes_prod2 --stride 50
"""
import argparse
import glob
import json
from pathlib import Path

import numpy as np

from fes import KJ2KCAL, crossings, fes2d, read_colvar, total_bias

GRID = (np.linspace(0.48, 1.56, 109), np.linspace(0.0, 0.52, 109))
BW = (0.02, 0.01)


def aligned_sd(Fs, ref, lo=6.0, min_est=2):
    """Mean over grid points of the SD (kcal/mol) across per-walker FES estimates, on
    points where ref < lo kcal/mol and at least min_est walkers sampled (2: the active
    start's walkers split into two disjoint regions in the second half). Each walker is
    first shifted to match ref over the points it shares with ref."""
    m = np.isfinite(ref) & (ref * KJ2KCAL < lo)
    S = []
    for F in Fs:
        ok = m & np.isfinite(F)
        S.append(np.where(ok, F - (np.mean(F[ok] - ref[ok]) if ok.any() else 0), np.nan))
    S = np.array(S) * KJ2KCAL
    n = np.isfinite(S).sum(0)
    use = m & (n >= min_est)
    if use.sum() < 5:
        return None, int(use.sum())
    return float(np.nanmean(np.nanstd(S[:, use], axis=0, ddof=1))), int(use.sum())


def block_rms(Fa, Fb, lo=6.0):
    m = np.isfinite(Fa) & np.isfinite(Fb) & (Fa * KJ2KCAL < lo) & (Fb * KJ2KCAL < lo)
    if m.sum() < 5:
        return None, int(m.sum())
    D = (Fa - Fb)[m] * KJ2KCAL
    return float(np.sqrt(np.mean((D - D.mean()) ** 2))), int(m.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--boxes", default="results/endpoints.json")
    ap.add_argument("--stride", type=int, default=50)
    a = ap.parse_args()
    boxes = json.load(open(a.boxes))["opes_proposal"]["basin_boxes_for_crossing_count"]
    (a.out / "figures").mkdir(parents=True, exist_ok=True)
    cdir = a.out / "colvar_strided"
    cdir.mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    report, data = {}, {}
    for start in ("active", "inactive"):
        files = sorted(glob.glob(f"{a.run}/{start}/w*/COLVAR.*"))
        cv = read_colvar(files)
        data[start] = cv
        b = total_bias(cv)
        x, y, wk, t = cv["d_tm36"], cv["rmsd_npxxy"], cv["walker"], cv["time"]
        F_all = fes2d(x, y, b, GRID, BW)
        half = {}
        for name, sel in (("first_half", t <= t.max() / 2), ("second_half", t > t.max() / 2)):
            ref = fes2d(x[sel], y[sel], b[sel], GRID, BW)
            Fw = [fes2d(x[sel & (wk == k)], y[sel & (wk == k)], b[sel & (wk == k)], GRID, BW)
                  for k in np.unique(wk)]
            sd, n = aligned_sd(Fw, ref)
            half[name] = {"walker_sd_kcal": sd, "grid_points": n}
        q = [(t > t.max() * i / 4) & (t <= t.max() * (i + 1) / 4) for i in range(4)]
        Fq3, Fq4 = (fes2d(x[m], y[m], b[m], GRID, BW) for m in (q[2], q[3]))
        rms34, n34 = block_rms(Fq3, Fq4)
        i2a, a2i, ev = crossings(cv, boxes["inactive"], boxes["active"])
        # deepest excursion toward the other basin, in both CVs at once
        other = boxes["inactive"] if start == "active" else boxes["active"]
        dx = np.maximum(0, np.maximum(other[0] - x, x - other[1]))
        dy = np.maximum(0, np.maximum(other[2] - y, y - other[3]))
        dist = np.hypot(dx / BW[0], dy / BW[1])
        j = int(np.argmin(dist))
        sec = {}
        for c in ("d_hbond", "na_site", "chi_w295", "lig_cont", "d_salt"):
            sec[c] = {"corr_with_CV1": round(float(np.corrcoef(cv[c], x)[0, 1]), 3),
                      "corr_with_CV2": round(float(np.corrcoef(cv[c], y)[0, 1]), 3),
                      "min": round(float(cv[c].min()), 3), "max": round(float(cv[c].max()), 3)}
        report[start] = {
            "walker_error": half,
            "block_q3_vs_q4_rms_kcal": rms34, "block_q3_vs_q4_points": n34,
            "crossings_i2a_a2i": [i2a, a2i],
            "closest_approach_to_other_basin": {
                "walker": int(wk[j]), "time_ns": round(float(t[j]) / 1000, 2),
                "d_tm36": round(float(x[j]), 3), "rmsd_npxxy": round(float(y[j]), 3),
                "distance_in_bandwidths": round(float(dist[j]), 2)},
            "frac_frames_with_CV1_in_other_box_range": round(float(
                ((x >= other[0]) & (x <= other[1])).mean()), 4),
            "frac_frames_with_CV2_in_other_box_range": round(float(
                ((y >= other[2]) & (y <= other[3])).mean()), 4),
            "secondary_observables": sec,
            "rct_last_kJ": round(float(cv["opes.rct"][-1]), 2),
        }
        for k, f in enumerate(files):
            m = np.where(wk == k)[0][::a.stride]
            cols = [c for c in cv if c != "walker"]
            np.savetxt(cdir / f"{start}_w{k}.dat", np.column_stack([cv[c][m] for c in cols]),
                       fmt="%.5g", header="FIELDS " + " ".join(cols)
                       + f"  (every {a.stride}th line of {f})")

        fig, ax = plt.subplots(2, 1, figsize=(11, 5.5), sharex=True, constrained_layout=True)
        for k in np.unique(wk):
            m = wk == k
            ax[0].plot(t[m][::10] / 1000, x[m][::10], lw=0.4, label=f"walker {k}")
            ax[1].plot(t[m][::10] / 1000, y[m][::10], lw=0.4)
        for i, (lab, bx) in enumerate((("inactive", boxes["inactive"]), ("active", boxes["active"]))):
            col = "tab:blue" if lab == "inactive" else "tab:red"
            ax[0].axhspan(bx[0], bx[1], color=col, alpha=0.12, label=f"{lab} box")
            ax[1].axhspan(bx[2], bx[3], color=col, alpha=0.12)
        for e in ev:
            for axx in ax:
                axx.axvline(e["time_ps"] / 1000, color="k", lw=0.8)
        ax[0].set_ylabel("CV1 R167–T281 Cα (nm)"); ax[1].set_ylabel("CV2 NPxxYA RMSD (nm)")
        ax[1].set_xlabel("time per walker (ns)")
        ax[0].set_title(f"{start} start — {i2a} inactive→active, {a2i} active→inactive "
                        "committed crossings (both CVs in box)")
        ax[0].legend(ncol=6, fontsize=7, loc="upper right")
        fig.savefig(a.out / "figures" / f"timeseries_{start}.png", dpi=130)
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 5.5), constrained_layout=True)
    for start, col in (("active", "tab:red"), ("inactive", "tab:blue")):
        cv = data[start]
        ax.scatter(cv["rmsd_npxxy"][::40], cv["d_tm36"][::40], s=0.3, alpha=0.25, c=col,
                   label=f"{start} start", rasterized=True)
    for lab, bx in boxes.items():
        ax.add_patch(Rectangle((bx[2], bx[0]), bx[3] - bx[2], bx[1] - bx[0], fill=False,
                               ec="k", lw=1.2))
        ax.text(bx[2], bx[1] + 0.01, f"{lab} basin", fontsize=8)
    ax.set_xlabel("CV2 NPxxYA RMSD (nm)"); ax.set_ylabel("CV1 R167–T281 Cα (nm)")
    ax.set_title("Visited CV space by start (every 40th frame, all walkers)")
    ax.legend(markerscale=20)
    fig.savefig(a.out / "figures" / "paths_overlay.png", dpi=150)
    (a.out / "convergence.json").write_text(json.dumps(report, indent=1))
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
