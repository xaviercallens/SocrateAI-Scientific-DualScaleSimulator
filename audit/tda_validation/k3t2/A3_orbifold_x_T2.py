#!/usr/bin/env python
"""
A3: DIRECT product (T^4/Z2) x T^2 as a 6-dimensional cellular complex
(product cells, d(a x b) = da x b + (-1)^|a| a x db), homology over Z/2, Z/3, Z/5
by sparse rank computation mod p. No Kunneth formula is used in the computation.
T^2 factor is the cubical Z_M^2 complex (M >= 2, NONZERO differential; the minimal
1-2-1 CW with zero differential is deliberately not used).
Writes A3_results.json. Run: prlimit --as=8589934592 -- <venv-python> A3_orbifold_x_T2.py [--big]
  --big adds N=6 x Z_2^2 (~166k cells).
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_cells import (cubical_torus, quotient_by_involution, neg_action, product, homology_mod_p,
                       homology_mod_p_colred, rss_mb)

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = json.load(open(os.path.join(HERE, "expectations.json")))["partA"]["A3_orbifold_x_T2"]
FIELDS = (2, 3, 5)
configs = [(2, 3), (4, 2), (4, 3)]
if "--big" in sys.argv:
    configs.append((6, 2))

out = {"tier": "B", "script": "A3_orbifold_x_T2.py", "argv": sys.argv[1:], "runs": []}
for N, M in configs:
    t0 = time.time()
    cc, index = cubical_torus(N, 4)
    q, info, _, _ = quotient_by_involution(cc, index, neg_action(N))
    T2, _ = cubical_torus(M, 2)
    P = product(q, T2)
    tb = time.time() - t0
    d2 = len(P.check_d2())
    rec = {"orbifold_N": N, "T2_M": M, "fvector": P.fvector(), "n_cells": len(P), "chi_fvector": P.chi(),
           "d2_violations": d2, "build_seconds": round(tb, 2), "by_field": {}}
    for p in FIELDS:
        b2, s2 = homology_mod_p_colred(P, p, True)
        r = {"betti_engine2": b2, "ranks_d_k": s2["ranks_d_k"], "sec_engine2": s2["seconds"]}
        if len(P) <= 60000 or "--both" in sys.argv:
            b1, s1 = homology_mod_p(P, p, True)
            r.update({"betti_engine1": b1, "sec_engine1": s1["seconds"], "engines_agree": b1 == b2})
        e = EXP["expected"].get(f"F{p}")
        chib = sum((-1) ** k * x for k, x in enumerate(b2))
        r["chi_betti"] = chib
        r["expected"] = e if e else "not pre-specified (F2; chi=0 only)"
        r["PASS"] = (b2 == e if e else True) and r.get("engines_agree", True) and chib == 0
        rec["by_field"][f"F{p}"] = r
        print(N, M, p, b2, r.get("engines_agree"), s2["seconds"], "s", flush=True)
    rec["PASS"] = all(v["PASS"] for v in rec["by_field"].values()) and d2 == 0 and rec["chi_fvector"] == 0
    rec["seconds_total"] = round(time.time() - t0, 2)
    rec["maxrss_MB"] = rss_mb()
    out["runs"].append(rec)
    print(N, M, rec["fvector"], "PASS" if rec["PASS"] else "FAIL", rec["seconds_total"], "s", rec["maxrss_MB"], "MB", flush=True)
out["all_PASS"] = all(r["PASS"] for r in out["runs"])
json.dump(out, open(os.path.join(HERE, "A3_results.json"), "w"), indent=1)
print("all_PASS", out["all_PASS"])
