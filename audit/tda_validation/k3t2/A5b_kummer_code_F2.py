#!/usr/bin/env python
"""
A5b: F2 sensitivity of the chain-level Kummer construction to WHICH points are resolved.
Predictions pre-registered in expectations_A5b.json (Kummer code RM(1,4)), committed before this ran.
Tier B (exact) against tier-L predictions. Writes A5b_results.json.
Run: prlimit --as=8589934592 -- <venv-python> A5b_kummer_code_F2.py
"""
import itertools, json, os, random, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_cells import cubical_torus, product, homology_mod_p_colred, rss_mb
from lib_kummer import build_orbifold_pieces, bockstein_phi2, assemble

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = json.load(open(os.path.join(HERE, "expectations_A5b.json")))["cases"]
N = 6
T0 = time.time()
q, fixed, U, links, checks, _ = build_orbifold_pieces(N)
phi2s = {p: bockstein_phi2(q, links[p])[0] for p in fixed}
bits = [tuple(x // (N // 2) for x in p) for p in fixed]
hyperplanes = []
for a in itertools.product((0, 1), repeat=4):
    if any(a):
        for c in (0, 1):
            hyperplanes.append(frozenset(i for i, b in enumerate(bits) if sum(x * y for x, y in zip(a, b)) % 2 == c))
assert len(set(hyperplanes)) == 30
rng = random.Random(20260919)
while True:
    J8 = frozenset(rng.sample(range(16), 8))
    if J8 not in hyperplanes:
        break
idx = {b: i for i, b in enumerate(bits)}
plane = {idx[(x0, x1, 0, 0)] for x0 in (0, 1) for x1 in (0, 1)}
indep = {idx[b] for b in [(0, 0, 0, 0), (1, 0, 0, 0), (0, 1, 0, 0), (0, 0, 1, 0)]}
cases = {
    "J_affine_hyperplane_x0_eq_0_k8": set(range(16)) - {i for i, b in enumerate(bits) if b[0] == 0},
    "J_non_hyperplane_8set_k8": set(range(16)) - set(J8),
    "resolved_affine_2plane_k4": plane,
    "resolved_4_affinely_independent_k4": indep,
    "resolved_single_point_k1": {idx[(0, 0, 0, 0)]},
}
T2, _ = cubical_torus(2, 2)
out = {"tier": "B vs pre-registered tier-L predictions", "script": "A5b_kummer_code_F2.py", "seed": 20260919, "cases": {}}
for name, resolved in cases.items():
    kinds = ["resolve" if i in resolved else "cone" for i in range(16)]
    X, _ = assemble(q, fixed, U, links, kinds, phi2s)
    P = product(X, T2)
    r = {"resolved_points_bits": sorted(bits[i] for i in resolved), "k": len(resolved),
         "X": {f"F{p}": homology_mod_p_colred(X, p) for p in (2, 3)},
         "X_ranks_d_k": {f"F{p}": homology_mod_p_colred(X, p, True)[1]["ranks_d_k"] for p in (2, 3)},
         "X_x_T2": {f"F{p}": homology_mod_p_colred(P, p) for p in (2, 3)},
         "X_x_T2_ranks_d_k": {f"F{p}": homology_mod_p_colred(P, p, True)[1]["ranks_d_k"] for p in (2, 3)},
         "chi_X": X.chi(), "chi_XxT2": P.chi(), "d2_X": len(X.check_d2()), "expected": EXP[name]}
    def kun(a):
        o = [0] * (len(a) + 2)
        for i, x in enumerate(a):
            for j, y in enumerate((1, 2, 1)):
                o[i + j] += x * y
        return o
    r["PASS"] = (r["X"]["F2"] == EXP[name]["F2"] and r["X"]["F3"] == EXP[name]["F3"] and
                 r["X_x_T2"]["F2"] == kun(EXP[name]["F2"]) and r["X_x_T2"]["F3"] == kun(EXP[name]["F3"]) and r["d2_X"] == 0)
    out["cases"][name] = r
    print(name, r["X"], r["X_x_T2"], "PASS" if r["PASS"] else "FAIL", flush=True)
out["all_PASS"] = all(r["PASS"] for r in out["cases"].values())
out["seconds_total"] = round(time.time() - T0, 2)
out["maxrss_MB"] = rss_mb()
json.dump(out, open(os.path.join(HERE, "A5b_results.json"), "w"), indent=1)
print("all_PASS", out["all_PASS"], out["seconds_total"])
