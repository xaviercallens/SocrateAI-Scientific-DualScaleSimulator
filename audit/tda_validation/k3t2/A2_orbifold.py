#!/usr/bin/env python
"""
A2: the orbifold T^4/Z2 (v -> -v on Z_N^4, N even) as a cellular complex
(orbits of the cubical cells, coinvariant boundary), homology over Z/2, Z/3, Z/5.
Cross-checks (tier B, exact):
  * the involution fixes exactly 16 cells, all vertices, pointwise (no cell of dim>=1
    fixed setwise for N even) -> the orbit complex IS the CW complex of the quotient;
  * d^2 = 0 on the quotient; chi(f-vector) = chi(Betti);
  * Lefschetz: sum_k (-1)^k tr(sigma_# | C_k(T^4)) computed at chain level;
  * transfer: the action of sigma on H_k(T^4;F_3) computed on explicit cycles
    (coordinate subtori z_I): (a) z_I are cycles, (b) they are independent modulo
    boundaries, (c) sigma(z_I) - (-1)^k z_I is a boundary (rank tests), so
    dim H_k(T^4;F3)^G is computed, then compared with H_k(T^4/G;F3);
  * N = 2, 4, 6 (N-independence); engine 1 vs engine 2.
Writes A2_results.json. Run: prlimit --as=8589934592 -- <venv-python> A2_orbifold.py
"""
import itertools, json, os, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from lib_cells import (cubical_torus, quotient_by_involution, neg_action, homology_mod_p,
                       homology_mod_p_colred, rss_mb)

HERE = os.path.dirname(os.path.abspath(__file__))
EXP = json.load(open(os.path.join(HERE, "expectations.json")))["partA"]["A2_orbifold_T4_mod_Z2"]
FIELDS = (2, 3, 5)


def sparse_rank(cols, p):
    """rank over F_p of a list of sparse columns (dict row->int); column reduction."""
    piv = {}
    r = 0
    for c in cols:
        col = {k: v % p for k, v in c.items() if v % p}
        while col:
            low = max(col)
            if low in piv:
                o = piv[low]
                f = (col[low] * pow(o[low], p - 2, p)) % p
                for x, v in o.items():
                    nv = (col.get(x, 0) - f * v) % p
                    if nv:
                        col[x] = nv
                    else:
                        col.pop(x, None)
            else:
                piv[low] = col
                r += 1
                break
    return r


def transfer_check(N, p=3):
    cc, index = cubical_torus(N, 4)
    act = neg_action(N)
    n = len(cc)
    # chain-level Lefschetz number
    trace = [0] * 5
    for i, (v, S) in enumerate(cc.label):
        (w, S2), s = act(v, S)
        if index[(w, S2)] == i:
            trace[len(S)] += s
    L = sum((-1) ** k * t for k, t in enumerate(trace))

    def sigma(chain):
        out = {}
        for c, x in chain.items():
            (w, S2), s = act(*cc.label[c])
            j = index[(w, S2)]
            out[j] = out.get(j, 0) + s * x
        return out

    def bd_chain(chain):
        out = {}
        for c, x in chain.items():
            for f, y in cc.bd[c].items():
                out[f] = out.get(f, 0) + x * y
        return {k: v for k, v in out.items() if v % p}

    res = {"N": N, "p": p, "chain_traces_sigma_by_dim": trace, "lefschetz_number_chain_level": L, "degrees": {}}
    for k in range(5):
        B = [cc.bd[c] for c in range(n) if cc.dims[c] == k + 1]
        rB = sparse_rank(B, p)
        Z = []
        for I in itertools.combinations(range(4), k):
            z = {}
            for vals in itertools.product(range(N), repeat=k):
                v = [0] * 4
                for i, x in zip(I, vals):
                    v[i] = x
                z[index[(tuple(v), I)]] = 1
            Z.append(z)
        cycles_ok = all(not bd_chain(z) for z in Z)
        indep = sparse_rank(B + Z, p) - rB
        sign = (-1) ** k
        diffs_are_boundaries = []
        for z in Z:
            sz = sigma(z)
            w = dict(sz)
            for c, x in z.items():
                w[c] = w.get(c, 0) - sign * x
            w = {a: b for a, b in w.items() if b % p}
            diffs_are_boundaries.append(sparse_rank(B + [w], p) == rB if w else True)
        # invariant dimension: action is sign * identity on the span of the z_I
        inv_dim = len(Z) if sign == 1 else 0
        res["degrees"][k] = {"n_cycles_z_I": len(Z), "all_cycles": cycles_ok,
                             "independent_mod_boundaries": indep,
                             "sigma_acts_as": sign, "sigma_z_minus_sign_z_is_boundary_all": all(diffs_are_boundaries),
                             "dim_invariants": inv_dim}
    res["invariant_betti_F3"] = [res["degrees"][k]["dim_invariants"] for k in range(5)]
    return res


out = {"tier": "B", "script": "A2_orbifold.py", "runs": [], "transfer": None}
for N in (2, 4, 6):
    t0 = time.time()
    cc, index = cubical_torus(N, 4)
    q, info, _, _ = quotient_by_involution(cc, index, neg_action(N))
    rec = {"N": N, "cover_fvector": cc.fvector(), "quotient_fvector": q.fvector(), "n_cells": len(q),
           "n_fixed_cells": info["n_fixed_cells"], "fixed_cell_dims": info["fixed_cell_dims"],
           "fixed_cells": [list(v) for v, S in info["fixed_cells"]],
           "d2_violations": len(q.check_d2()), "chi_fvector": q.chi(), "by_field": {}}
    for p in FIELDS:
        b1, s1 = homology_mod_p(q, p, True)
        b2, s2 = homology_mod_p_colred(q, p, True)
        e = EXP["expected_odd_p"].get(f"F{p}")
        chib = sum((-1) ** k * x for k, x in enumerate(b1))
        rec["by_field"][f"F{p}"] = {"betti": b1, "betti_engine2": b2, "ranks_d_k": s2["ranks_d_k"],
                                    "engines_agree": b1 == b2, "chi_betti": chib,
                                    "expected": e if e else "not pre-specified (chi=8 only)",
                                    "sec": [s1["seconds"], s2["seconds"]],
                                    "PASS": (b1 == e if e else True) and b1 == b2 and chib == EXP["chi"]}
    rec["seconds"] = round(time.time() - t0, 2)
    rec["maxrss_MB"] = rss_mb()
    rec["PASS"] = all(v["PASS"] for v in rec["by_field"].values()) and rec["d2_violations"] == 0 \
        and rec["n_fixed_cells"] == 16 and rec["fixed_cell_dims"] == [0] and rec["chi_fvector"] == EXP["chi"]
    # 2-primary torsion count from UCT: b_k(F2) = r_k + t_k + t_{k-1}, r_k := b_k(F3) (F3 = F5 checked)
    r = rec["by_field"]["F3"]["betti"]
    f2 = rec["by_field"]["F2"]["betti"]
    t = []
    prev = 0
    for k in range(5):
        tk = f2[k] - r[k] - prev
        t.append(tk)
        prev = tk
    rec["derived_2torsion_summands_per_degree_H_k(Z)"] = t
    rec["derived_note"] = ("UCT: b_k(F2) = rank_k + t_k + t_{k-1}; rank_k taken from F3 (F3 = F5 found). "
                           "t_k = number of cyclic 2-primary summands in H_k(T^4/Z2; Z). Derived, not a separate computation.")
    out["runs"].append(rec)
    print(N, rec["quotient_fvector"], {f: v["betti"] for f, v in rec["by_field"].items()},
          "fixed", rec["n_fixed_cells"], "t2", t, "PASS" if rec["PASS"] else "FAIL", rec["seconds"], "s", flush=True)

tr = transfer_check(4, 3)
tr["matches_quotient_F3"] = tr["invariant_betti_F3"] == out["runs"][1]["by_field"]["F3"]["betti"]
tr["lefschetz_expected"] = EXP["lefschetz_number_sigma_on_T4"]
tr["PASS"] = tr["matches_quotient_F3"] and tr["lefschetz_number_chain_level"] == 16 and all(
    d["all_cycles"] and d["independent_mod_boundaries"] == d["n_cycles_z_I"] and d["sigma_z_minus_sign_z_is_boundary_all"]
    for d in tr["degrees"].values())
out["transfer"] = tr
print("transfer", tr["invariant_betti_F3"], "L=", tr["lefschetz_number_chain_level"], "PASS" if tr["PASS"] else "FAIL")

# external comparison: simplicial Freudenthal track (read-only, other worktree, committed file)
ext = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity_v2/D-tda/02_invariance_and_quotient_results.json"
try:
    d = json.load(open(ext))
    out["external_freudenthal_track_file"] = ext
    out["external_freudenthal_track_excerpt"] = {k: v for k, v in d.items() if "quot" in k.lower() or "N=" in k}
except Exception as e:  # report, do not repair
    out["external_freudenthal_track_error"] = repr(e)
out["all_PASS"] = all(r["PASS"] for r in out["runs"]) and tr["PASS"]
out["maxrss_MB"] = rss_mb()
json.dump(out, open(os.path.join(HERE, "A2_results.json"), "w"), indent=1, default=str)
print("all_PASS", out["all_PASS"])
