"""v0.3 efficacy readout, exactly as pre-registered in PLAN_v03.md.

Per replica (unbiased, COLVAR_measure every 10 ps):
  active        CV1 d_tm36 >= 1.00 nm AND CV2 rmsd_npxxy >= 0.25 nm
  frac_active   fraction of frames active in the analysis window (default 100-500 ns)
  first_exit    start of the first run of >= 10 ns contiguous non-active frames
                (None = never; right-censored)
  gates         ligand bound: d_salt <= 0.6 nm in >= 95% of frames, and lig_cont over
                the last 50 ns >= 50% of its 0-10 ns mean. A failing replica is reported
                and excluded from the comparison, never silently dropped.
  nb_retention  (Nb39 systems) mean nb_cont over the last 100 ns / mean over 0-10 ns

Per pair (DAMGO vs naloxone, same Nb39 condition):
  difference of mean frac_active, SE of the difference = sqrt(SE_dam^2 + SE_nlx^2) with
  SE = SD/sqrt(n) over valid replicas; one-sided Mann-Whitney U (DAMGO > naloxone);
  rank-matched pairings (i-th highest DAMGO vs i-th highest naloxone).
  verdict:  "discrimination"  diff > 2 SE and DAMGO higher in >= 3 of 4 rank pairings
            "reverse"         diff < -2 SE (naloxone holds the active state better)
            "inconclusive"    otherwise (not the same as "no effect")
  The primary readout is the free_* pair; nb39_* is secondary (PLAN_v03.md).

Usage: python analysis/efficacy.py runs/v03 --out results/v03 [--window 100 500]
"""
import argparse
import glob
import json
from pathlib import Path

import numpy as np
from scipy.stats import mannwhitneyu

CV1_MIN, CV2_MIN = 1.00, 0.25
EXIT_NS = 10.0
PAIRS = {"free": ("free_dam", "free_nlx"), "nb39": ("nb39_dam", "nb39_nlx")}


def read(path):
    with open(path) as f:
        fields = f.readline().split()[2:]
    d = np.loadtxt(path, comments="#")
    return {k: d[:, i] for i, k in enumerate(fields)}


def runs_of(mask):
    """(start, end) index pairs of True runs."""
    m = np.concatenate([[False], mask, [False]]).astype(int)
    e = np.flatnonzero(np.diff(m))
    return list(zip(e[::2], e[1::2]))


def replica(path, window):
    cv = read(path)
    t = cv["time"] / 1000.0                                  # ns
    dt = float(np.median(np.diff(t)))
    act = (cv["d_tm36"] >= CV1_MIN) & (cv["rmsd_npxxy"] >= CV2_MIN)
    w = (t >= window[0]) & (t <= window[1])
    out = {"file": path, "ns": round(float(t[-1]), 2),
           "window_complete": bool(t[-1] >= window[1] - dt),
           "frac_active": round(float(act[w].mean()), 4) if w.any() else None,
           "frac_active_all": round(float(act.mean()), 4)}
    first = None
    for a, b in runs_of(~act):
        if (b - a) * dt >= EXIT_NS:
            first = round(float(t[a]), 2)
            break
    out["first_persistent_exit_ns"] = first
    salt_ok = float((cv["d_salt"] <= 0.6).mean())
    lc0 = cv["lig_cont"][t <= 10].mean()
    lc1 = cv["lig_cont"][t >= t[-1] - 50].mean()
    out["gate"] = {"d_salt_le_0.6_frac": round(salt_ok, 4),
                   "lig_cont_last50_over_first10": round(float(lc1 / lc0), 3)}
    out["valid"] = bool(salt_ok >= 0.95 and lc1 >= 0.5 * lc0)
    if "nb_cont" in cv:
        n0 = cv["nb_cont"][t <= 10].mean()
        n1 = cv["nb_cont"][t >= t[-1] - 100].mean()
        out["nb_retention"] = round(float(n1 / n0), 3)
    out["mean_last100"] = {k: round(float(cv[k][t >= t[-1] - 100].mean()), 3)
                           for k in ("d_tm36", "rmsd_npxxy", "d_hbond", "na_site")}
    return out, (t, cv, act)


def compare(dam, nlx):
    a = np.array([r["frac_active"] for r in dam if r["valid"]])
    b = np.array([r["frac_active"] for r in nlx if r["valid"]])
    res = {"n_valid": [int(a.size), int(b.size)]}
    if a.size < 2 or b.size < 2:
        res["verdict"] = "insufficient valid replicas"
        return res
    se = np.sqrt(a.var(ddof=1) / a.size + b.var(ddof=1) / b.size)
    diff = float(a.mean() - b.mean())
    ranked = int(sum(x > y for x, y in zip(sorted(a, reverse=True), sorted(b, reverse=True))))
    U, p = mannwhitneyu(a, b, alternative="greater")
    res.update({"mean_frac_active": [round(float(a.mean()), 4), round(float(b.mean()), 4)],
                "diff_dam_minus_nlx": round(diff, 4), "se_diff": round(float(se), 4),
                "rank_pairings_dam_higher": f"{ranked}/{min(a.size, b.size)}",
                "mannwhitney_U": float(U), "mannwhitney_p_one_sided": round(float(p), 4)})
    if se == 0:
        res["verdict"] = "discrimination" if diff > 0 and ranked >= 3 else "inconclusive"
    elif diff > 2 * se and ranked >= 3:
        res["verdict"] = "discrimination"
    elif diff < -2 * se:
        res["verdict"] = "reverse"
    else:
        res["verdict"] = "inconclusive"
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--window", type=float, nargs=2, default=(100.0, 500.0))
    a = ap.parse_args()
    (a.out / "figures").mkdir(parents=True, exist_ok=True)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rep, series = {"definition": {"active": f"d_tm36>={CV1_MIN} & rmsd_npxxy>={CV2_MIN}",
                                  "window_ns": list(a.window), "exit_ns": EXIT_NS}}, {}
    for sysn in [s for p in PAIRS.values() for s in p]:
        files = sorted(glob.glob(f"{a.run}/{sysn}/r*/COLVAR_measure*"))
        if not files:
            continue
        rows = []
        for f in files:
            r, s = replica(f, a.window)
            rows.append(r)
            series.setdefault(sysn, []).append(s)
        rep[sysn] = rows
    rep["comparison"] = {}
    for name, (d, n) in PAIRS.items():
        if d in rep and n in rep:
            c = compare(rep[d], rep[n])
            c["incomplete_replicas"] = [r["file"] for r in rep[d] + rep[n]
                                        if not r["window_complete"]]
            rep["comparison"][name] = c
            if "nb_retention" in rep[d][0]:
                c["nb_retention_mean"] = [round(float(np.mean([r["nb_retention"] for r in rep[x]])), 3)
                                          for x in (d, n)]
    rep["primary_readout"] = "free"

    for sysn, ss in series.items():
        fig, ax = plt.subplots(2, 1, figsize=(10, 5), sharex=True, constrained_layout=True)
        for k, (t, cv, act) in enumerate(ss):
            ax[0].plot(t[::20], cv["d_tm36"][::20], lw=0.5, label=f"r{k}")
            ax[1].plot(t[::20], cv["rmsd_npxxy"][::20], lw=0.5)
        ax[0].axhline(CV1_MIN, color="k", ls="--", lw=0.8)
        ax[1].axhline(CV2_MIN, color="k", ls="--", lw=0.8)
        ax[0].set_ylabel("CV1 TM3–TM6 (nm)"); ax[1].set_ylabel("CV2 NPxxYA RMSD (nm)")
        ax[1].set_xlabel("time (ns)"); ax[0].legend(ncol=4, fontsize=7)
        ax[0].set_title(f"{sysn}: dashed = active-state thresholds (active = both above)")
        fig.savefig(a.out / "figures" / f"timeseries_{sysn}.png", dpi=130)
        plt.close(fig)
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    names = [s for s in ("nb39_dam", "nb39_nlx", "free_dam", "free_nlx") if s in rep]
    for i, s in enumerate(names):
        v = [r["frac_active"] for r in rep[s] if r["frac_active"] is not None]
        ax.scatter([i] * len(v), v, c=["tab:blue" if r["valid"] else "tab:gray" for r in rep[s]
                                       if r["frac_active"] is not None], zorder=3)
        ax.hlines(np.mean(v), i - 0.25, i + 0.25, color="k")
    ax.set_xticks(range(len(names)), names)
    ax.set_ylabel(f"fraction active, {a.window[0]:.0f}–{a.window[1]:.0f} ns")
    ax.set_title("per replica (grey = failed ligand gate); bar = mean")
    fig.savefig(a.out / "figures" / "frac_active.png", dpi=130)
    (a.out / "efficacy.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep["comparison"], indent=1))


if __name__ == "__main__":
    main()
