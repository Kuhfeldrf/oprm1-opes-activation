"""Bilayer embedding QC for one frame or a trajectory (runbook §8 failure signatures:
"the TM bundle loses its fold", "receptor exits the membrane", mis-built membrane).

Per frame:
  thickness_PP_A    mean z separation of upper/lower leaflet phosphorus (POPC "PC" P31)
  apl_A2            (box area - receptor footprint) / lipids per leaflet (POPC + CHL;
                    mixture APL). Footprint = 1 A grid cells within 2 A of any receptor
                    heavy atom in the |z - mid| < 5 A slab. apl_box_A2 is the uncorrected value.
  rec_center_dz_A   receptor TM-core centre of mass minus bilayer midplane (z)
  tm_tilt_deg       angle between the TM-bundle principal axis and the membrane normal
  lipid_in_bundle   lipid heavy atoms within 5 A of the TM-bundle axis and inside the
                    hydrophobic slab (lipids should surround, not enter, the bundle)
  water_in_core     water oxygens inside the hydrophobic core (|z - mid| < 10 A),
                    excluding those within 12 A of the bundle axis (receptor pore water
                    near the sodium pocket is physical)
  tm_helix_frac     fraction of TM-core residues in alpha-helix (MDAnalysis DSSP)

Usage: python analysis/membrane_qc.py TPR TRAJ_OR_GRO [--stride N] [--json out.json]
Reference ranges for POPC:CHOL 7:3 at 310 K: P-P ~ 40-44 A; APL ~ 48-55 A^2.
"""
import argparse
import json
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "cv"))
from scaffold import SCAFFOLD  # noqa: E402

OFFSET = 67   # topology residue index (1-based) + 67 = canonical P35372 number
TM_CORE = sorted(set(SCAFFOLD) | set(range(275, 300)) | set(range(316, 334)))


def footprint(xyz, mid, box, slab=5.0, r=2.0):
    """Receptor cross-sectional area (A^2) at the bilayer midplane on a 1 A grid."""
    s = xyz[np.abs(xyz[:, 2] - mid) < slab][:, :2] % box[:2]
    gx, gy = np.meshgrid(np.arange(0, box[0], 1.0), np.arange(0, box[1], 1.0), indexing="ij")
    G = np.column_stack([gx.ravel(), gy.ravel()])
    occ = np.zeros(len(G), bool)
    for p in s:
        d = np.abs(G - p)
        d = np.minimum(d, box[:2] - d)
        occ |= (d ** 2).sum(1) < r * r
    return float(occ.sum())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("tpr"); ap.add_argument("traj")
    ap.add_argument("--stride", type=int, default=1)
    ap.add_argument("--json", default=None)
    ap.add_argument("--no-dssp", action="store_true")
    a = ap.parse_args()
    u = mda.Universe(a.tpr, a.traj)
    # GROMACS (single rank) wraps atoms into the box at neighbour-search steps, so once
    # the receptor diffuses across a box edge it is split in the trajectory. A split
    # receptor puts the bundle centre/axis inside the lipids and corrupts DSSP. Make the
    # receptor+ligand whole, centre it, and wrap everything else by residue.
    from MDAnalysis import transformations as trans
    solute = u.select_atoms("not resname PA PC OL CHL WAT Na+ Cl-")
    rest = u.atoms - solute
    u.trajectory.add_transformations(trans.unwrap(solute), trans.center_in_box(solute),
                                     trans.wrap(rest, compound="residues"))
    P = u.select_atoms("resname PC and name P31")
    lip_heavy = u.select_atoms("resname PA PC OL CHL and not element H")
    rec_heavy = u.select_atoms("not resname PA PC OL CHL WAT Na+ Cl- DAM and not element H")
    ow = u.select_atoms("resname WAT and name O")
    resids = " ".join(str(n - OFFSET) for n in TM_CORE)
    tm_ca = u.select_atoms(f"(resid {resids}) and name CA and not resname DAM WAT PA PC OL CHL")
    n_lip_leaf = None
    rows = []
    try:
        from MDAnalysis.analysis.dssp import DSSP
        have_dssp = not a.no_dssp
    except ImportError:
        have_dssp = False
    for ts in u.trajectory[::a.stride]:
        z = P.positions[:, 2]
        mid = z.mean()
        up, lo = z[z > mid], z[z <= mid]
        if n_lip_leaf is None:
            heads = u.select_atoms("(resname PC and name P31) or (resname CHL and name O1)")
            hz = heads.positions[:, 2]
            n_lip_leaf = ((hz > mid).sum(), (hz <= mid).sum())
        box = ts.dimensions
        X = tm_ca.positions
        c = X.mean(0)
        w, v = np.linalg.eigh(np.cov((X - c).T))
        axis = v[:, -1] * np.sign(v[2, -1])
        tilt = float(np.degrees(np.arccos(abs(axis[2]))))
        L = lip_heavy.positions
        d_axis = np.linalg.norm(np.cross(L - c, axis), axis=1)
        in_slab = np.abs(L[:, 2] - mid) < 15
        Wp = ow.positions
        d_w = np.linalg.norm(np.cross(Wp - c, axis), axis=1)
        rows.append({
            "time_ns": float(ts.time / 1000),
            "thickness_PP_A": float(up.mean() - lo.mean()),
            "apl_A2": float((box[0] * box[1] - footprint(rec_heavy.positions, mid, box))
                            / np.mean(n_lip_leaf)),
            "apl_box_A2": float(box[0] * box[1] / np.mean(n_lip_leaf)),
            "rec_center_dz_A": float(c[2] - mid),
            "tm_tilt_deg": tilt,
            "lipid_in_bundle": int(((d_axis < 5) & in_slab).sum()),
            "water_in_core": int(((np.abs(Wp[:, 2] - mid) < 10) & (d_w > 12)).sum()),
        })
    if have_dssp:
        prot = u.select_atoms(f"resid {resids} and not resname DAM")
        d = DSSP(prot).run(step=a.stride)
        frac = (np.array(d.results.dssp) == "H").mean(axis=1)
        for r, f in zip(rows, frac):
            r["tm_helix_frac"] = float(f)
    keys = [k for k in rows[0] if k != "time_ns"]
    summary = {"frames": len(rows), "lipids_per_leaflet": [int(x) for x in n_lip_leaf]}
    for k in keys:
        vals = np.array([r[k] for r in rows])
        summary[k] = {"mean": float(vals.mean()), "min": float(vals.min()), "max": float(vals.max()),
                      "last": float(vals[-1])}
    flags = []
    if not 36 <= summary["thickness_PP_A"]["last"] <= 48:
        flags.append("P-P thickness outside 36-48 A")
    if summary["lipid_in_bundle"]["max"] > 0:
        flags.append("lipid heavy atoms inside the TM bundle")
    if abs(summary["rec_center_dz_A"]["last"]) > 8:
        flags.append("receptor TM core off-centre in the bilayer by > 8 A")
    if summary["water_in_core"]["last"] > 20:
        flags.append("water inside the hydrophobic core (bilayer defect)")
    if "tm_helix_frac" in summary and summary["tm_helix_frac"]["last"] < 0.8:
        flags.append("TM helicity < 0.8 (fold loss)")
    summary["flags"] = flags
    print(json.dumps(summary, indent=1))
    if a.json:
        Path(a.json).write_text(json.dumps({"summary": summary, "frames": rows}, indent=1))


if __name__ == "__main__":
    main()
