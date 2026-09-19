#!/usr/bin/env python
"""
A5: negative controls for the chain-level K3 x T^2 route.
  (a) partial resolution X_k: resolve k of the 16 singular points (Cyl(phi)), cone the rest,
      k = 0, 8, 15, 16 (points chosen by a seeded random permutation, seed 20260919),
      homology of X_k and of the DIRECT product X_k x T^2 over F2, F3, F5.
      The value k only changes WHICH complexes are glued in; no formula in k is typed.
  (b) wrong gluing class (phi_2 = 0 at all 16 points): must differ from K3 over F2.
  (c) T^2 replaced by S^2 (cubical boundary of [-1,1]^3): orbifold x S^2, K3 x S^2.
  (d) T^2 replaced by the Klein bottle (T^2(Z_4 x Z_3)/tau): K3 x Klein, orbifold x Klein over F2 vs F3.
Tier B. Writes A5_results.json. Run: prlimit --as=8589934592 -- <venv-python> A5_negative_controls.py
"""
import json, os, random, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_cells import (cubical, cubical_torus, cubical_sphere, quotient_by_involution, klein_action, product,
                       homology_mod_p_colred, rss_mb)
from lib_kummer import build_orbifold_pieces, bockstein_phi2, assemble

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = json.load(open(os.path.join(HERE, "expectations.json")))["partA"]["A5_controls"]
FIELDS = (2, 3, 5)
N = 6
SEED = 20260919
T0 = time.time()
out = {"tier": "B", "script": "A5_negative_controls.py", "N": N, "seed": SEED}

q, fixed, U, links, checks, _ = build_orbifold_pieces(N)
phi2s = {p: bockstein_phi2(q, links[p])[0] for p in fixed}
T2, _ = cubical_torus(2, 2)
S2, _ = cubical_sphere(3)
Tk, ik = cubical([4, 3], [True, True])
KB, _, _, _ = quotient_by_involution(Tk, ik, klein_action(2, 3))


def hom(C):
    r = {"fvector": C.fvector(), "n_cells": len(C), "chi_fvector": C.chi(), "d2_violations": len(C.check_d2())}
    for p in FIELDS:
        t = time.time()
        b, st = homology_mod_p_colred(C, p, True)
        r[f"F{p}"] = b
        r[f"ranks_d_k_F{p}"] = st["ranks_d_k"]
        r[f"sec_F{p}"] = round(time.time() - t, 2)
    r["chi_betti_F3"] = sum((-1) ** k * x for k, x in enumerate(r["F3"]))
    r["maxrss_MB"] = rss_mb()
    return r


rng = random.Random(SEED)
order = list(range(16))
rng.shuffle(order)
# (a) partial resolution
out["partial_resolution"] = {}
for k in (0, 8, 15, 16):
    chosen = set(order[:k])
    kinds = ["resolve" if i in chosen else "cone" for i in range(16)]
    X, _ = assemble(q, fixed, U, links, kinds, phi2s)
    e = EXP["partial_resolution_X_k"][str(k)]
    rX = hom(X)
    rP = hom(product(X, T2))
    rec = {"k": k, "resolved_points": [list(fixed[i]) for i in sorted(chosen)], "X_k": rX, "X_k_x_T2": rP,
           "expected_X_k_F3": e["X_k_F3"], "expected_X_k_x_T2_F3": e["X_k_x_T2_F3"]}
    rec["PASS"] = (rX["F3"] == e["X_k_F3"] and rX["F5"] == e["X_k_F3"] and rP["F3"] == e["X_k_x_T2_F3"]
                   and rP["F5"] == e["X_k_x_T2_F3"] and rX["chi_fvector"] == e["chi_X_k"]
                   and rP["chi_fvector"] == 0 and rX["d2_violations"] == 0 and rP["d2_violations"] == 0)
    out["partial_resolution"][str(k)] = rec
    print("k", k, "X_k", {f: rX[f] for f in ("F2", "F3", "F5")}, "xT2", {f: rP[f] for f in ("F2", "F3")},
          "PASS" if rec["PASS"] else "FAIL", flush=True)
vals = [tuple(out["partial_resolution"][str(k)]["X_k_x_T2"]["F3"]) for k in (0, 8, 15, 16)]
out["k_series_distinct_F3"] = len(set(vals)) == 4
out["k_series_rule"] = EXP["failure_rule"]

# (b) wrong gluing class
Xw, _ = assemble(q, fixed, U, links, ["trivial2"] * 16, phi2s)
rw = hom(Xw)
rwp = hom(product(Xw, T2))
out["wrong_class_phi2_zero"] = {"X": rw, "X_x_T2": rwp,
                                "F3_equals_K3": rw["F3"] == [1, 0, 22, 0, 1],
                                "F2_differs_from_K3": rw["F2"] != [1, 0, 22, 0, 1],
                                "F2_x_T2_differs_from_K3xT2": rwp["F2"] != [1, 2, 23, 44, 23, 2, 1],
                                "note": EXP["wrong_gluing_class_phi2_eq_0_all16"]["meaning"]}
out["wrong_class_phi2_zero"]["PASS"] = (out["wrong_class_phi2_zero"]["F3_equals_K3"] and
                                        out["wrong_class_phi2_zero"]["F2_differs_from_K3"] and
                                        out["wrong_class_phi2_zero"]["F2_x_T2_differs_from_K3xT2"])
print("wrong class", {f: rw[f] for f in ("F2", "F3")}, {f: rwp[f] for f in ("F2", "F3")}, flush=True)

# (c)/(d) replace T^2
Xr, _ = assemble(q, fixed, U, links, ["resolve"] * 16, phi2s)
Xo, _ = assemble(q, fixed, U, links, ["cone"] * 16, phi2s)
exS, exK = EXP["S2_instead_of_T2"], EXP["Klein_instead_of_T2"]
out["factor_S2_Klein"] = {}
out["factor_S2_Klein"]["S2_alone"] = hom(S2)
out["factor_S2_Klein"]["Klein_alone"] = hom(KB)
for name, A, B, ex in (("orbifold_x_S2", Xo, S2, {"F3": exS["orbifold_x_S2_F3"], "F5": exS["orbifold_x_S2_F3"]}),
                       ("K3_x_S2", Xr, S2, {"F2": exS["K3_x_S2_F2"], "F3": exS["K3_x_S2_F3"], "F5": exS["K3_x_S2_F3"]}),
                       ("orbifold_x_Klein", Xo, KB, {"F3": exK["orbifold_x_Klein_F3"], "F5": exK["orbifold_x_Klein_F3"]}),
                       ("K3_x_Klein", Xr, KB, {"F2": exK["K3_x_Klein"]["F2"], "F3": exK["K3_x_Klein"]["F3"],
                                               "F5": exK["K3_x_Klein"]["F3"]})):
    r = hom(product(A, B))
    r["expected"] = ex
    r["PASS"] = all(r[f] == v for f, v in ex.items()) and r["d2_violations"] == 0
    out["factor_S2_Klein"][name] = r
    print(name, r["fvector"], {f: r[f] for f in ("F2", "F3", "F5")}, "PASS" if r["PASS"] else "FAIL", flush=True)
kk = out["factor_S2_Klein"]["K3_x_Klein"]
out["factor_S2_Klein"]["K3xKlein_F2_equals_K3xT2_F2"] = kk["F2"] == [1, 2, 23, 44, 23, 2, 1]
out["factor_S2_Klein"]["K3xKlein_F3_differs_from_K3xT2"] = kk["F3"] != [1, 2, 23, 44, 23, 2, 1]

out["all_PASS"] = (all(r["PASS"] for r in out["partial_resolution"].values()) and out["k_series_distinct_F3"]
                   and out["wrong_class_phi2_zero"]["PASS"] and
                   all(v["PASS"] for k, v in out["factor_S2_Klein"].items() if isinstance(v, dict) and "PASS" in v))
out["seconds_total"] = round(time.time() - T0, 2)
out["maxrss_MB"] = rss_mb()
json.dump(out, open(os.path.join(HERE, "A5_results.json"), "w"), indent=1)
print("all_PASS", out["all_PASS"], out["seconds_total"], "s", out["maxrss_MB"], "MB")
