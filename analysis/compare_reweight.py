"""Cross-check analysis/fes.py against PLUMED REWEIGHT_BIAS + HISTOGRAM (runbook §7.6).

Both use the total bias and the same bandwidths. PLUMED's grid (GRID_BIN=100 points per
axis + 1) differs from fes.py's, so fes.py is interpolated (linear, on finite points) onto
PLUMED's grid. Reports, per start: position of each global minimum, and the RMS / max
difference (kcal/mol, constant removed) where both are < 8 kcal/mol.

Usage: python analysis/compare_reweight.py results/fes_prod2
"""
import json
import sys
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator

KJ2KCAL = 1 / 4.184


def main():
    out = Path(sys.argv[1])
    gx, gy = np.linspace(0.48, 1.56, 109), np.linspace(0.0, 0.52, 109)
    rep = {}
    for s in ("active", "inactive"):
        p = out / "plumed_reweight" / f"fes_plumed_{s}.dat"
        if not p.exists():
            rep[s] = "missing"
            continue
        d = np.loadtxt(p, comments="#")
        x, y = d[:, 0], d[:, 1]
        Fp = d[:, 2] * KJ2KCAL
        Fp[~np.isfinite(Fp)] = np.nan
        Fp -= np.nanmin(Fp)
        Ff = np.loadtxt(out / "fes" / f"{s}.dat")          # already kcal/mol, min 0
        fi = RegularGridInterpolator((gx, gy), Ff, bounds_error=False, fill_value=np.nan)
        Fq = fi(np.column_stack([x, y]))
        m = np.isfinite(Fp) & np.isfinite(Fq) & (Fp < 8) & (Fq < 8)
        D = Fp[m] - Fq[m]
        D -= D.mean()
        ip, iq = np.nanargmin(Fp), np.nanargmin(Fq)
        rep[s] = {"points_compared": int(m.sum()),
                  "rms_diff_kcal": round(float(np.sqrt(np.mean(D ** 2))), 3),
                  "max_abs_diff_kcal": round(float(np.abs(D).max()), 3),
                  "min_plumed": [round(float(x[ip]), 3), round(float(y[ip]), 3)],
                  "min_fes_py": [round(float(x[iq]), 3), round(float(y[iq]), 3)]}
    (out / "plumed_reweight" / "comparison.json").write_text(json.dumps(rep, indent=1))
    print(json.dumps(rep, indent=1))


if __name__ == "__main__":
    main()
