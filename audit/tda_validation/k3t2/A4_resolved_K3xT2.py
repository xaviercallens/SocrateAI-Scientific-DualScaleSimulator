#!/usr/bin/env python
"""
A4: resolved Kummer K3 and K3 x T^2 at chain level.
  X_res = U  union_{16 L_p}  16 x Cyl(phi_p : C(L_p) -> C(S^2))      (see lib_kummer.py)
  then the DIRECT product X_res x T^2 (cubical Z_M^2), homology by rank mod p.
Everything computed; the single tier-L input is the class of phi_2 (generator of
H^2(RP^3;Z)=Z/2, Gysin with Euler number +-2), whose Lefschetz-duality consequence
H_*(Cyl,L) = (0,0,1,0,1) is checked here at chain level.
Controls inside this script: link = RP^3 (all 16, three fields); U + 16 cones must
reproduce the direct A2 orbifold numbers; wrong-class local model (phi_2 = 0).
Writes A4_results.json. Run: prlimit --as=8589934592 -- <venv-python> A4_resolved_K3xT2.py [--big]
  --big also runs X_res x T^2(Z_3^2) and engine 1 on the products.
"""
import json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_cells import (cubical_torus, product, homology_mod_p, homology_mod_p_colred, subcomplex, rss_mb)
from lib_kummer import build_orbifold_pieces, bockstein_phi2, assemble, local_model_relative

HERE = os.path.dirname(os.path.abspath(__file__))
EXPA = json.load(open(os.path.join(HERE, "expectations.json")))["partA"]
EXP = EXPA["A4_resolved_K3_and_K3xT2"]
A2 = json.load(open(os.path.join(HERE, "A2_results.json")))
FIELDS = (2, 3, 5)
N = 6
T0 = time.time()
out = {"tier": "B (with one labelled tier-L input: class of phi_2)", "script": "A4_resolved_K3xT2.py",
       "argv": sys.argv[1:], "N": N}

q, fixed, U, links, checks, sizes = build_orbifold_pieces(N)
out["decomposition_checks"] = checks
out["sizes"] = sizes
print("checks", checks, sizes, flush=True)

# ---- links and local models --------------------------------------------------
link_recs = []
phi2s = {}
for p in fixed:
    L, _ = subcomplex(q, links[p])
    hb = {f"F{pp}": homology_mod_p_colred(L, pp) for pp in FIELDS}
    phi2, diag = bockstein_phi2(q, links[p])
    phi2s[p] = phi2
    lm = {}
    for kind in ("resolve", "trivial2"):
        R, _ = local_model_relative(q, links[p], kind, phi2)
        lm[kind] = {f"F{pp}": homology_mod_p_colred(R, pp) for pp in FIELDS}
    rec = {"p": list(p), "L_fvector": L.fvector(), "L_betti": hb, "phi2_diag": diag, "H_rel_Cyl_L": lm}
    rec["L_is_RP3_homology"] = hb == {k: v for k, v in EXPA["engine_controls"]["RP3"].items()}
    lmx = EXP["local_model_checks"]
    rec["local_model_PASS"] = (lm["resolve"]["F2"] == lmx["H_*(Cyl(pi), L; F2)"] and
                               lm["resolve"]["F3"] == lmx["H_*(Cyl(pi), L; F3)"] and
                               lm["trivial2"]["F2"] == lmx["wrong_class_control_phi2_eq_0"]["H_*(Cyl,L;F2)"] and
                               lm["trivial2"]["F3"] == lmx["wrong_class_control_phi2_eq_0"]["H_*(Cyl,L;F3)"] and
                               diag["integral_cocycle_violations"] == 0 and diag["phi2_mod2_not_coboundary"])
    link_recs.append(rec)
out["links"] = link_recs
out["all_links_RP3"] = all(r["L_is_RP3_homology"] for r in link_recs)
out["all_local_models_PASS"] = all(r["local_model_PASS"] for r in link_recs)
print("links RP3:", out["all_links_RP3"], "local models:", out["all_local_models_PASS"],
      link_recs[0]["H_rel_Cyl_L"], link_recs[0]["phi2_diag"], flush=True)


def hom_all(C, name, expected=None, engines=(2, 1)):
    rec = {"name": name, "fvector": C.fvector(), "n_cells": len(C), "chi_fvector": C.chi(),
           "d2_violations": len(C.check_d2()), "by_field": {}}
    for p in FIELDS:
        r = {}
        if 2 in engines:
            b, s = homology_mod_p_colred(C, p, True)
            r.update({"betti": b, "ranks_d_k": s["ranks_d_k"], "sec_engine2": s["seconds"]})
        if 1 in engines:
            b1, s1 = homology_mod_p(C, p, True)
            r.update({"betti_engine1": b1, "sec_engine1": s1["seconds"]})
            if "betti" not in r:
                r["betti"] = b1
            r["engines_agree"] = b1 == r["betti"]
        r["chi_betti"] = sum((-1) ** k * x for k, x in enumerate(r["betti"]))
        e = expected.get(f"F{p}") if expected else None
        r["expected"] = e
        r["PASS"] = (e is None or r["betti"] == e) and r.get("engines_agree", True) and r["chi_betti"] == rec["chi_fvector"]
        rec["by_field"][f"F{p}"] = r
        print(" ", name, f"F{p}", r["betti"], "PASS" if r["PASS"] else "FAIL", r.get("sec_engine2"), r.get("sec_engine1"), flush=True)
    rec["PASS"] = all(v["PASS"] for v in rec["by_field"].values()) and rec["d2_violations"] == 0
    rec["maxrss_MB"] = rss_mb()
    return rec


# ---- control: U + 16 cones must equal the direct orbifold (A2, N=6) ------------
Xc, _ = assemble(q, fixed, U, links, ["cone"] * 16, phi2s)
a2 = next(r for r in A2["runs"] if r["N"] == N)
exp_cone = {f: a2["by_field"][f]["betti"] for f in ("F2", "F3", "F5")}
out["cone_reassembly"] = hom_all(Xc, "U + 16 cones (must equal A2 direct orbifold)", exp_cone)

# ---- resolved K3 ----------------------------------------------------------------
t = time.time()
Xr, ginfo = assemble(q, fixed, U, links, ["resolve"] * 16, phi2s)
out["K3"] = hom_all(Xr, "X_res = U + 16 Cyl(phi) (Kummer K3)", {f"F{p}": EXP["K3"][f"F{p}"] for p in FIELDS})
out["K3"]["poincare_symmetric_all_fields"] = all(v["betti"] == v["betti"][::-1] for v in out["K3"]["by_field"].values())
out["K3"]["chi_expected"] = EXP["K3"]["chi"]
out["K3"]["build_seconds"] = round(time.time() - t, 2)

# ---- K3 x T^2, direct product --------------------------------------------------
Ms = [2, 3] if "--big" in sys.argv else [2]
out["K3xT2"] = []
for M in Ms:
    t = time.time()
    T2, _ = cubical_torus(M, 2)
    P = product(Xr, T2)
    eng = (2, 1) if "--big" in sys.argv and M == 2 else (2,)
    rec = hom_all(P, f"X_res x T^2(Z_{M}^2) direct product", {f"F{p}": EXP["K3xT2"][f"F{p}"] for p in FIELDS}, eng)
    rec["poincare_symmetric_all_fields"] = all(v["betti"] == v["betti"][::-1] for v in rec["by_field"].values())
    rec["seconds_total"] = round(time.time() - t, 2)
    out["K3xT2"].append(rec)
    del P

out["all_PASS"] = (all(checks.values()) and out["all_links_RP3"] and out["all_local_models_PASS"] and
                   out["cone_reassembly"]["PASS"] and out["K3"]["PASS"] and out["K3"]["poincare_symmetric_all_fields"]
                   and all(r["PASS"] for r in out["K3xT2"]))
out["seconds_total"] = round(time.time() - T0, 2)
out["maxrss_MB"] = rss_mb()
json.dump(out, open(os.path.join(HERE, "A4_results.json"), "w"), indent=1)
print("all_PASS", out["all_PASS"], out["seconds_total"], "s", out["maxrss_MB"], "MB")
