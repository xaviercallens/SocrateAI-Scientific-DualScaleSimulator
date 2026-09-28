#!/usr/bin/env python
"""
Part B exploration log (tier X): fixed list of extra K3 / K3xT2 / null configurations
(beyond the pre-registered N-grid studies), every outcome recorded including memory/time failures.
Same worker, same 8 GiB cap, per-run timeout B_TIMEOUT. Resumable (chunked calls).
Run: B_TIMEOUT=560 B_CALL_BUDGET=300 timeout 595 prlimit --as=8589934592 -- <venv-python> B_explore.py
Writes B_explore.json.
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from B_driver import run

HERE = os.path.dirname(os.path.abspath(__file__))
K3 = dict(space="K3", method="rips", expected=[1, 0, 22], save_diagram_dims=[2])
SPECS = [
    dict(K3, N=4000, tau=0.8, seed=0, fields=[3, 2]),
    dict(K3, N=4000, tau=0.8, seed=1, fields=[3, 2]),
    dict(K3, N=4000, tau=0.8, seed=2, fields=[3, 2]),
    dict(K3, N=4000, tau=0.8, seed=0, fields=[3], fps=True),
    dict(K3, N=3000, tau=0.9, seed=0, fields=[3]),
    dict(K3, N=2000, tau=1.0, seed=0, fields=[3]),
    dict(K3, N=2000, tau=1.2, seed=0, fields=[3]),
    dict(K3, N=4000, tau=0.95, seed=0, fields=[3]),
    dict(space="null4", method="rips", expected=[1, 0, 0], null_D=16, null_diam=1.414, save_diagram_dims=[2],
         N=4000, tau=0.8, seed=0, fields=[3, 2]),
    dict(space="K3xT2", method="rips", expected=[1, 2, 23], torus_scale=0.5, save_diagram_dims=[2],
         N=4000, tau=0.8, seed=0, fields=[3]),
    dict(space="K3xT2", method="rips", expected=[1, 2, 23], torus_scale=0.5, save_diagram_dims=[2],
         N=8000, tau=0.8, seed=0, fields=[3]),
    # appended after spec 8 (diameter-matched null) timed out: density-matched null, see worker 'null4m'
    dict(space="null4m", method="rips", expected=[1, 0, 0], save_diagram_dims=[2], N=4000, tau=0.8, seed=0, fields=[3, 2]),
    # appended: K3 x T2 was cheap at N<=8000 (6-dim sample is sparse at tau=0.8), so push N
    dict(space="K3xT2", method="rips", expected=[1, 2, 23], torus_scale=0.5, save_diagram_dims=[2], N=16000, tau=0.8, seed=0, fields=[3]),
    dict(space="K3xT2", method="rips", expected=[1, 2, 23], torus_scale=0.5, save_diagram_dims=[2], N=32000, tau=0.8, seed=0, fields=[3]),
    dict(space="K3xT2", method="rips", expected=[1, 2, 23], torus_scale=0.5, save_diagram_dims=[2], N=64000, tau=0.8, seed=0, fields=[3]),
    # appended: re-runs of the largest completed N of the orbifold / T4 grid studies (same spec and seed as in
    # B_orbifold_rips.json / B_T4_rips.json) to record the beta(eps) curve added to the worker afterwards
    dict(space="orbifold", method="rips", tau=1.2, expected=[1, 0, 6], N=12800, seed=0, fields=[3, 2]),
    dict(space="T4", method="rips", tau=1.2, expected=[1, 4, 6, 4, 1], N=12800, seed=0, fields=[3, 2]),
    # appended after the advisor review: witness complexes at a smaller relaxation tau, where the tau=0.8 witness
    # study (B_K3_witness.json) timed out at its first grid point
    dict(space="K3", method="witness", witness_ratio=20, tau=0.4, expected=[1, 0, 22], N=20000, seed=0, fields=[3], save_diagram_dims=[2]),
    dict(space="K3", method="witness", witness_ratio=20, tau=0.5, expected=[1, 0, 22], N=40000, seed=0, fields=[3], save_diagram_dims=[2]),
    dict(space="T3", method="witness", witness_ratio=20, tau=0.6, expected=[1, 3, 3, 1], N=16000, seed=0, fields=[3]),
]
path = os.path.join(HERE, "B_explore.json")
out = json.load(open(path)) if os.path.exists(path) else {"tier": "X", "runs": {}}
t0 = time.time()
budget = float(os.environ.get("B_CALL_BUDGET", "1e9"))
for i, s in enumerate(SPECS):
    key = str(i)
    if key in out["runs"]:
        continue
    if time.time() - t0 > budget:
        print("paused"); break
    r = run(s)
    out["runs"][key] = r
    json.dump(out, open(path, "w"), indent=1)
    if r.get("failed"):
        print(i, s["space"], s["N"], s["tau"], "FAILED", r["reason"], r["stderr_tail"][-80:], flush=True)
    else:
        f = r["by_field"]["F3"]
        print(i, s["space"], s["N"], s["tau"], r["info"]["n_simplices"], r["wall"], "s", r["maxrss_MB"], "MB",
              "ratio", round(f["best_window_ratio"], 3), "curve", {k: v for k, v in list(f["beta_curve"].items())[-6:]}, flush=True)
out["complete"] = all(str(i) in out["runs"] for i in range(len(SPECS)))
json.dump(out, open(path, "w"), indent=1)
