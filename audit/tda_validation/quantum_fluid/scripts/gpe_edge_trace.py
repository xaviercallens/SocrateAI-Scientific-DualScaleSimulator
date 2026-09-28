"""POST-HOC trace of the E2a discrepancies (b0 > N_winding): count winding
vortices whose centre lies OUTSIDE the r < 0.7 R_eff analysis disk but within
0.7 R_eff + delta (delta = 0.25, 0.5 a_ho), i.e. vortices whose low-density core
straddles the disk edge and is therefore seen by the lower-star filtration of
the masked disk but not by the all-corners-inside winding count.
Command: prlimit --as=8589934592 -- .venv-tda/bin/python gpe_edge_trace.py
Output: ../results/gpe_edge_trace.json
"""
import glob, json, os, sys
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qf_common as qc  # noqa: E402

out = {}
for fn in sorted(glob.glob(os.path.join(qc.DATA_ROOT, "gpe", "psi_Omega*.npz"))):
    d = np.load(fn); psi = d["psi"].astype(np.complex128); x = d["x"]; g = float(d["g"]); Om = float(d["Omega"])
    dx = x[1] - x[0]
    R = np.sqrt(2 * np.sqrt(g / np.pi) / np.sqrt(1 - Om ** 2))
    q = qc.plaquette_winding(np.angle(psi)); q[-1, :] = 0; q[:, -1] = 0
    jj, ii = np.nonzero(q)
    r = np.hypot(x[ii] + dx / 2, x[jj] + dx / 2)
    out[os.path.basename(fn)[4:-4]] = {"R_eff": float(R),
        "n_vortex_centres_in_(0.7R, 0.7R+0.25]": int(((r > 0.7 * R) & (r <= 0.7 * R + 0.25)).sum()),
        "n_vortex_centres_in_(0.7R, 0.7R+0.5]": int(((r > 0.7 * R) & (r <= 0.7 * R + 0.5)).sum())}
json.dump(out, open(os.path.join(qc.RESULTS, "gpe_edge_trace.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
