"""Runbook §6.4: CV1/CV2 and the R167-T281 sidechain distance on deposited
coordinates of every candidate structure. CV2 = RMSD of NPxxYA CA after
superposing the TM1-TM5 scaffold CA onto the inactive reference (no refit on
NPxxYA) -- the same split as the dual-weight PLUMED reference (§7.2).

Usage: python cv/measure_deposited.py structures/raw REF_PDBID REF_CHAIN"""
import sys
from pathlib import Path
import numpy as np
import gemmi
sys.path.insert(0, str(Path(__file__).parent))
from scaffold import SCAFFOLD, NPXXYA, CV1_PAIR, HBOND_PAIR

RECEPTOR = {"9MQH": ["A"], "9MQJ": ["A"], "9MQI": ["A"], "9PXU": ["R"],
            "10TM": ["R", "F"], "9PPQ": ["R"], "8F7Q": ["R", "M"], "8Y72": ["R"]}
STATE = {"9MQH": "inactive", "9MQJ": "inactive", "9MQI": "inactive",
         "9PXU": "inactive", "10TM": "active", "8F7Q": "active", "8Y72": "active",
         "9PPQ": "unknown"}

def atoms(ch, name="CA"):
    out = {}
    for r in ch:
        if r.het_flag == "A":
            a = r.find_atom(name, "*")
            if a:
                out[r.seqid.num] = np.array(a.pos.tolist())
    return out

def kabsch(P, Q):
    """Rotation/translation mapping P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, S, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, pc, qc

raw = Path(sys.argv[1]); ref_id, ref_ch = sys.argv[2], sys.argv[3]
def chain(pid, c):
    return gemmi.read_structure(str(raw / f"{pid}.cif"))[0][c]
ref = atoms(chain(ref_id, ref_ch))
print(f"# reference for CV2: {ref_id}/{ref_ch}; units Angstrom")
print(f"{'pdb':5} {'ch':2} {'state':8} {'CV1_CA':>7} {'CV2_rmsd':>8} {'scaf_rmsd':>9} {'nscaf':>5} {'CZ-OG1':>7}")
for pid, chs in RECEPTOR.items():
    for c in chs:
        ch = chain(pid, c)
        ca = atoms(ch)
        cv1 = np.linalg.norm(ca[CV1_PAIR[0]] - ca[CV1_PAIR[1]])
        scaf = [n for n in SCAFFOLD if n in ca and n in ref]
        P = np.array([ca[n] for n in scaf]); Q = np.array([ref[n] for n in scaf])
        R, pc, qc = kabsch(P, Q)
        fit = lambda X: (X - pc) @ R.T + qc
        srmsd = np.sqrt(((fit(P) - Q) ** 2).sum(1).mean())
        X = fit(np.array([ca[n] for n in NPXXYA])); Y = np.array([ref[n] for n in NPXXYA])
        cv2 = np.sqrt(((X - Y) ** 2).sum(1).mean())
        a = atoms(ch, HBOND_PAIR[0][1]).get(HBOND_PAIR[0][0])
        b = atoms(ch, HBOND_PAIR[1][1]).get(HBOND_PAIR[1][0])
        hb = f"{np.linalg.norm(a - b):7.2f}" if a is not None and b is not None else "    n/a"
        print(f"{pid:5} {c:2} {STATE[pid]:8} {cv1:7.2f} {cv2:8.2f} {srmsd:9.2f} {len(scaf):5d} {hb}")
