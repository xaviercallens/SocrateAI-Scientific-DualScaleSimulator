#!/usr/bin/env python
"""
A1: cubical T^2, T^4, T^6 on Z_N^n (periodic), Betti over Z/2, Z/3, Z/5,
plus rank-engine controls (field sensitivity RP^2, RP^3; spheres; Klein bottle;
product routine vs direct torus; dense rank and GUDHI periodic cubical cross-checks).
Tier B (exact arithmetic mod p). Writes A1_results.json.
Run: prlimit --as=8589934592 -- <venv-python> A1_tori_and_engine_controls.py
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_cells import (cubical_torus, cubical_sphere, quotient_by_involution, antipodal_action_box,
                       klein_action, cubical, product, homology_mod_p, homology_mod_p_colred,
                       betti_dense, rss_mb)
import gudhi
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = json.load(open(os.path.join(HERE, "expectations.json")))["partA"]
FIELDS = (2, 3, 5)
out = {"tier": "B", "script": "A1_tori_and_engine_controls.py", "tori": [], "controls": {}}
t_start = time.time()


def full(cc, name, expected_by_field, dense=False):
    rec = {"name": name, "fvector": cc.fvector(), "n_cells": len(cc), "chi_fvector": cc.chi(),
           "d2_violations": len(cc.check_d2()), "by_field": {}}
    for p in FIELDS:
        b1, s1 = homology_mod_p(cc, p, True)
        b2, s2 = homology_mod_p_colred(cc, p, True)
        e = expected_by_field.get(f"F{p}") if isinstance(expected_by_field, dict) else expected_by_field
        r = {"betti_engine1": b1, "betti_engine2": b2, "ranks_d_k": s2["ranks_d_k"],
             "sec_engine1": s1["seconds"], "sec_engine2": s2["seconds"],
             "engines_agree": b1 == b2, "chi_betti": sum((-1) ** k * x for k, x in enumerate(b1)),
             "expected": e}
        if dense:
            r["betti_dense"] = betti_dense(cc, p)
            r["dense_agrees"] = r["betti_dense"] == b1
        r["PASS"] = (b1 == e) and b1 == b2 and r["chi_betti"] == rec["chi_fvector"] and rec["d2_violations"] == 0 \
            and (not dense or r["dense_agrees"])
        rec["by_field"][f"F{p}"] = r
    rec["PASS"] = all(v["PASS"] for v in rec["by_field"].values())
    rec["maxrss_MB"] = rss_mb()
    return rec


# ---- tori ------------------------------------------------------------------
for n, Ns in ((2, (2, 3, 4)), (4, (2, 3, 4)), (6, (2, 3))):
    for N in Ns:
        t0 = time.time()
        cc, _ = cubical_torus(N, n)
        rec = full(cc, f"T^{n} cubical Z_{N}^{n}", EXP["A1_tori"]["expected"][f"T^{n}"], dense=(len(cc) <= 700))
        # GUDHI periodic cubical complex (independent code): all top cells filtration 0
        g = gudhi.PeriodicCubicalComplex(top_dimensional_cells=np.zeros([N] * n),
                                         periodic_dimensions=[True] * n)
        gb = {}
        for p in (2, 3):
            g.compute_persistence(homology_coeff_field=p)
            gb[f"F{p}"] = [int(x) for x in g.betti_numbers()]
        rec["gudhi_periodic_cubical_betti"] = gb
        rec["gudhi_agrees"] = all(gb[f] == rec["by_field"][f]["betti_engine1"] for f in gb)
        rec["PASS"] = rec["PASS"] and rec["gudhi_agrees"]
        rec["seconds_total"] = round(time.time() - t0, 3)
        out["tori"].append(rec)
        print(rec["name"], rec["fvector"], {f: v["betti_engine1"] for f, v in rec["by_field"].items()},
              "gudhi", gb, "PASS" if rec["PASS"] else "FAIL", flush=True)

# ---- product routine vs direct construction --------------------------------
T4, _ = cubical_torus(2, 4)
T2, _ = cubical_torus(3, 2)
P = product(T4, T2)
rec = full(P, "T^4(Z_2) x T^2(Z_3) via product()", EXP["A1_tori"]["expected"]["T^6"])
out["controls"]["product_T4xT2"] = rec
print(rec["name"], rec["fvector"], rec["by_field"]["F3"]["betti_engine1"], rec["PASS"], flush=True)

# ---- field-sensitivity controls --------------------------------------------
ec = EXP["engine_controls"]
S2, i2 = cubical_sphere(3)
out["controls"]["S2"] = full(S2, "S^2 = boundary [-1,1]^3", {f"F{p}": [1, 0, 1] for p in FIELDS}, dense=True)
S3, i3 = cubical_sphere(4)
out["controls"]["S3"] = full(S3, "S^3 = boundary [-1,1]^4", {f"F{p}": ec["S3"] for p in FIELDS}, dense=True)
RP2, info2, _, _ = quotient_by_involution(S2, i2, antipodal_action_box(3))
out["controls"]["RP2"] = full(RP2, "RP^2 = S^2/antipodal (cubical)", ec["RP2"], dense=True)
out["controls"]["RP2"]["n_fixed_cells"] = info2["n_fixed_cells"]
RP3, info3, _, _ = quotient_by_involution(S3, i3, antipodal_action_box(4))
out["controls"]["RP3"] = full(RP3, "RP^3 = S^3/antipodal (cubical)", ec["RP3"], dense=True)
out["controls"]["RP3"]["n_fixed_cells"] = info3["n_fixed_cells"]
for Mx, Ny in ((1, 2), (2, 3)):
    T, it = cubical([2 * Mx, Ny], [True, True])
    K, ik, _, _ = quotient_by_involution(T, it, klein_action(Mx, Ny))
    exp_k = EXP["A5_controls"]["Klein_instead_of_T2"]["Klein"]
    rec = full(K, f"Klein = T^2(Z_{2*Mx} x Z_{Ny})/tau", {"F2": exp_k["F2"], "F3": exp_k["F3"], "F5": exp_k["F3"]}, dense=True)
    out["controls"][f"Klein_{Mx}_{Ny}"] = rec
for k, v in out["controls"].items():
    print(k, v["fvector"], {f: x["betti_engine1"] for f, x in v["by_field"].items()}, "PASS" if v["PASS"] else "FAIL")

out["all_PASS"] = all(r["PASS"] for r in out["tori"]) and all(r["PASS"] for r in out["controls"].values())
out["seconds_total"] = round(time.time() - t_start, 2)
out["maxrss_MB"] = rss_mb()
out["hardware"] = {"nproc": os.cpu_count(), "note": "GCP VM, 29 GB RAM, 8 vCPU; single-threaded pure Python"}
json.dump(out, open(os.path.join(HERE, "A1_results.json"), "w"), indent=1)
print("all_PASS", out["all_PASS"], "seconds", out["seconds_total"], "maxrss_MB", out["maxrss_MB"])
