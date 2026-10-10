"""OPES reweighting, start-independence, convergence and crossing analysis (runbook §7.6, §10).

Weights are w_t = exp(+V_t / kT), with V the TOTAL bias (opes + walls) at frame t. The
FES is -kT ln P on a grid (Gaussian-kernel weighted histogram), shifted so its minimum
is 0. Values are reported in kcal/mol; PLUMED's units are kJ/mol.

Usage:
  python analysis/fes.py --active runs/active/walker*/COLVAR --inactive runs/inactive/walker*/COLVAR \
      --out results --burn 0
Outputs: results/fes/{active,inactive}.dat, results/fes/quarters_<start>_<k>.dat,
         results/fes/summary.json, results/figures/*.png
"""
import argparse
import json
from pathlib import Path

import numpy as np

KB = 0.0083144626  # kJ/mol/K
T = 310.0
KT = KB * T
KJ2KCAL = 1 / 4.184


def read_colvar(paths):
    """Concatenate walkers. Returns dict of column -> array, plus 'walker' and 'time'."""
    cols, parts = None, []
    for k, p in enumerate(sorted(paths)):
        with open(p) as f:
            hdr = f.readline().split()[2:]
        d = np.loadtxt(p, comments="#")
        if d.ndim == 1:
            d = d[None]
        if cols is None:
            cols = hdr
        elif hdr != cols:
            raise SystemExit(f"FATAL: {p} has different columns")
        parts.append(np.column_stack([d, np.full(len(d), k)]))
    d = np.vstack(parts)
    out = {c: d[:, i] for i, c in enumerate(cols)}
    out["walker"] = d[:, -1].astype(int)
    return out


def total_bias(cv):
    b = np.zeros_like(cv["time"])
    for c in ("opes.bias", "uw.bias", "lw.bias"):
        if c in cv:
            b = b + cv[c]
    return b


def fes2d(x, y, bias, grid, bw):
    """Weighted Gaussian KDE on a grid. Returns F (kJ/mol, min 0) with NaN where unsampled."""
    gx, gy = grid
    if len(x) == 0:
        return np.full((len(gx), len(gy)), np.nan)
    logw = bias / KT
    logw -= logw.max()
    w = np.exp(logw)
    X, Y = np.meshgrid(gx, gy, indexing="ij")
    P = np.zeros_like(X)
    # chunk to bound memory
    for s in range(0, len(x), 20000):
        xs, ys, ws = x[s:s + 20000], y[s:s + 20000], w[s:s + 20000]
        P += np.einsum("n,in,jn->ij", ws,
                       np.exp(-0.5 * ((gx[:, None] - xs[None]) / bw[0]) ** 2),
                       np.exp(-0.5 * ((gy[:, None] - ys[None]) / bw[1]) ** 2))
    with np.errstate(divide="ignore"):
        F = -KT * np.log(P)
    F -= np.nanmin(F[np.isfinite(F)])
    # mask grid points with negligible support (unsampled)
    F[P < P.max() * 1e-8] = np.nan
    return F


def basin_free_energies(F, grid, basins):
    gx, gy = grid
    X, Y = np.meshgrid(gx, gy, indexing="ij")
    out = {}
    for name, (x0, x1, y0, y1) in basins.items():
        m = (X >= x0) & (X <= x1) & (Y >= y0) & (Y <= y1) & np.isfinite(F)
        out[name] = float(-KT * np.log(np.exp(-F[m] / KT).sum())) if m.any() else np.nan
    ref = out.get("inactive", np.nan)
    return {k: (v - ref) * KJ2KCAL for k, v in out.items()}  # kcal/mol relative to inactive


def basin_coverage(F, grid, basins):
    """Fraction of each basin box's grid points with a finite FES value."""
    gx, gy = grid
    X, Y = np.meshgrid(gx, gy, indexing="ij")
    out = {}
    for name, (x0, x1, y0, y1) in basins.items():
        m = (X >= x0) & (X <= x1) & (Y >= y0) & (Y <= y1)
        out[name] = round(float(np.isfinite(F[m]).mean()), 3) if m.any() else 0.0
    return out


def crossings(cv, inactive, active):
    """Committed transitions per walker: a visit to one basin followed by a visit to the
    other, with both CVs inside each basin's box. Returns (n_i2a, n_a2i, event list)."""
    x, y, wk, t = cv["d_tm36"], cv["rmsd_npxxy"], cv["walker"], cv["time"]
    def inside(b):
        return (x >= b[0]) & (x <= b[1]) & (y >= b[2]) & (y <= b[3])
    ina, act = inside(inactive), inside(active)
    i2a = a2i = 0
    events = []
    for k in np.unique(wk):
        sel = np.where(wk == k)[0]
        state = None
        for i in sel:
            s = "I" if ina[i] else "A" if act[i] else None
            if s and state and s != state:
                if s == "A":
                    i2a += 1
                else:
                    a2i += 1
                events.append({"walker": int(k), "time_ps": float(t[i]), "to": s})
            if s:
                state = s
    return i2a, a2i, events


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--active", nargs="+", required=True)
    ap.add_argument("--inactive", nargs="+", required=True)
    ap.add_argument("--out", type=Path, default=Path("results"))
    ap.add_argument("--burn", type=float, default=0.0, help="fraction of each walker to drop")
    ap.add_argument("--grid", default="0.55,1.65,0.0,0.6")
    ap.add_argument("--bins", type=int, default=110)
    ap.add_argument("--bw", default="0.02,0.01")
    ap.add_argument("--inactive-box", default="0.55,0.85,0.0,0.15",
                    help="x0,x1,y0,y1 in nm; set from measured unbiased endpoints")
    ap.add_argument("--active-box", default="1.10,1.65,0.28,0.6")
    a = ap.parse_args()
    g = [float(v) for v in a.grid.split(",")]
    grid = (np.linspace(g[0], g[1], a.bins), np.linspace(g[2], g[3], a.bins))
    bw = [float(v) for v in a.bw.split(",")]
    basins = {"inactive": [float(v) for v in a.inactive_box.split(",")],
              "active": [float(v) for v in a.active_box.split(",")]}
    (a.out / "fes").mkdir(parents=True, exist_ok=True)
    (a.out / "figures").mkdir(parents=True, exist_ok=True)

    summary, fes = {}, {}
    for start, paths in (("active", a.active), ("inactive", a.inactive)):
        cv = read_colvar(paths)
        if a.burn:
            keep = np.zeros(len(cv["time"]), bool)
            for k in np.unique(cv["walker"]):
                idx = np.where(cv["walker"] == k)[0]
                keep[idx[int(len(idx) * a.burn):]] = True
            cv = {c: v[keep] for c, v in cv.items()}
        b = total_bias(cv)
        F = fes2d(cv["d_tm36"], cv["rmsd_npxxy"], b, grid, bw)
        fes[start] = F
        np.savetxt(a.out / "fes" / f"{start}.dat", F * KJ2KCAL,
                   header=f"FES kcal/mol; rows d_tm36 {g[0]}..{g[1]} nm, cols rmsd_npxxy "
                          f"{g[2]}..{g[3]} nm, {a.bins}x{a.bins}")
        i2a, a2i, ev = crossings(cv, basins["inactive"], basins["active"])
        # time quarters (per walker), for the convergence check
        quarters = []
        tmax = {k: cv["time"][cv["walker"] == k].max() for k in np.unique(cv["walker"])}
        tmin = {k: cv["time"][cv["walker"] == k].min() for k in np.unique(cv["walker"])}
        # quarters of the analysed window (after any burn-in), per walker
        frac = np.array([(cv["time"][i] - tmin[cv["walker"][i]])
                         / (tmax[cv["walker"][i]] - tmin[cv["walker"][i]])
                         for i in range(len(cv["time"]))])
        cum = []
        for q in range(4):
            m = (frac >= q / 4) & (frac <= (q + 1) / 4) if q == 0 else \
                (frac > q / 4) & (frac <= (q + 1) / 4)
            Fq = fes2d(cv["d_tm36"][m], cv["rmsd_npxxy"][m], b[m], grid, bw)
            np.savetxt(a.out / "fes" / f"quarters_{start}_{q+1}.dat", Fq * KJ2KCAL)
            quarters.append(basin_free_energies(Fq, grid, basins))
            mc = frac <= (q + 1) / 4
            Fc = fes2d(cv["d_tm36"][mc], cv["rmsd_npxxy"][mc], b[mc], grid, bw)
            cum.append(basin_free_energies(Fc, grid, basins))
        summary[start] = {
            "frames": int(len(cv["time"])), "walkers": int(len(np.unique(cv["walker"]))),
            "aggregate_ns": float(sum(tmax[k] - tmin[k] for k in tmax) / 1000.0),
            "crossings_inactive_to_active": i2a, "crossings_active_to_inactive": a2i,
            "dG_active_minus_inactive_kcal": basin_free_energies(F, grid, basins)["active"],
            # NaN dG means a basin box has no sampled grid point in this start's FES
            "basins_sampled": basin_coverage(F, grid, basins),
            "dG_by_quarter_kcal": [q["active"] for q in quarters],
            "dG_cumulative_kcal": [q["active"] for q in cum],
            "rct_last_kJ": float(cv["opes.rct"][-1]) if "opes.rct" in cv else None,
            "min_lig_cont": float(cv["lig_cont"].min()) if "lig_cont" in cv else None,
            "max_d_salt_nm": float(cv["d_salt"].max()) if "d_salt" in cv else None,
        }
        (a.out / "fes" / f"crossings_{start}.json").write_text(json.dumps(ev, indent=0))

    # start-independence (§10 criterion 2): difference over the region both sampled
    D = (fes["active"] - fes["inactive"]) * KJ2KCAL
    both = np.isfinite(D)
    # restrict to the "sampled range": F < 8 kcal/mol in both
    rel = both & (fes["active"] * KJ2KCAL < 8) & (fes["inactive"] * KJ2KCAL < 8)
    D_aligned = D - np.nanmean(D[rel]) if rel.any() else D
    np.savetxt(a.out / "fes" / "difference_active_minus_inactive.dat", D_aligned)
    shared = both  # every grid point either start sampled, regardless of depth
    D_shared = D - np.nanmean(D[shared]) if shared.any() else D
    summary["start_independence"] = {
        "grid_points_sampled_by_both": int(shared.sum()),
        "rms_diff_where_both_sampled_kcal":
            float(np.sqrt(np.nanmean(D_shared[shared] ** 2))) if shared.any() else None,
        "grid_points_compared": int(rel.sum()),
        "max_abs_diff_kcal": float(np.nanmax(np.abs(D_aligned[rel]))) if rel.any() else None,
        "rms_diff_kcal": float(np.sqrt(np.nanmean(D_aligned[rel] ** 2))) if rel.any() else None,
        "note": "difference after removing the arbitrary per-FES constant; region where both "
                "FES < 8 kcal/mol",
    }
    (a.out / "fes" / "summary.json").write_text(json.dumps(summary, indent=1))
    print(json.dumps(summary, indent=1))

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        ext = [g[2], g[3], g[0], g[1]]
        fig, ax = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
        for k, (start, F) in enumerate(fes.items()):
            im = ax[k].imshow(F * KJ2KCAL, origin="lower", extent=ext, aspect="auto",
                              vmin=0, vmax=12, cmap="viridis")
            ax[k].contour(grid[1], grid[0], F * KJ2KCAL, levels=np.arange(0, 12, 1),
                          colors="w", linewidths=0.4)
            ax[k].set_title(f"FES, {start} start (kcal/mol)")
            ax[k].set_xlabel("CV2 NPxxYA RMSD (nm)"); ax[k].set_ylabel("CV1 R167–T281 Cα (nm)")
            fig.colorbar(im, ax=ax[k])
        im = ax[2].imshow(D_aligned, origin="lower", extent=ext, aspect="auto",
                          vmin=-3, vmax=3, cmap="RdBu_r")
        ax[2].set_title("active − inactive start (kcal/mol)")
        ax[2].set_xlabel("CV2 NPxxYA RMSD (nm)")
        fig.colorbar(im, ax=ax[2])
        fig.savefig(a.out / "figures" / "fes_both_starts.png", dpi=150)
    except ImportError:
        pass


if __name__ == "__main__":
    main()
