"""Pass/fail for the two-CV smoke test (runbook §7.4, four conditions; GPU offload is
not required here because the example runs on CPU)."""
import json
import sys

rows = [list(map(float, l.split())) for l in open(sys.argv[1]) if not l.startswith("#")]
exp = json.load(open(sys.argv[2]))
t, d, r, b1, b2 = zip(*rows)
checks = {
    "COLVAR_smoke written with all steps": len(rows) >= exp["min_rows"],
    "CV1 starts at the expected value": abs(d[0] - exp["cv1_start_nm"]) < exp["tol_start_nm"],
    "CV2 starts at the expected value": abs(r[0] - exp["cv2_start_nm"]) < exp["tol_start_nm"],
    "CV1 changes": max(d) - min(d) > exp["min_change_nm"],
    "CV2 changes": max(r) - min(r) > exp["min_change_nm"],
    "CV1 moves toward its restraint (1.00 nm)": abs(d[-1] - 1.00) < abs(d[0] - 1.00),
    "bias r1 non-zero": max(b1) > 0,
    "bias r2 non-zero": max(b2) > 0,
}
for k, v in checks.items():
    print(f"{'PASS' if v else 'FAIL'}  {k}")
print(f"CV1 {d[0]:.4f} -> {d[-1]:.4f} nm;  CV2 {r[0]:.4f} -> {r[-1]:.4f} nm")
sys.exit(0 if all(checks.values()) else 1)
