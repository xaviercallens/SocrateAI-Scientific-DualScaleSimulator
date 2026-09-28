#!/usr/bin/env python
"""Track D v3, step 0: premise checks for the hybrid Mayer-Vietoris route, N in {4,6,8}.

Run:  cd audit/k3t2_rigidity_v3/D-tda && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python 00_premise.py 3 4 5 6 7 8
(arguments = the grid sizes N; the committed run used 3 4 5 6 7 8; odd N are negative controls, not valid quotients)

Checks per N (all counts computed):
 (1) fixed points of v -> -v on Z_N^4; orbit collisions inside a simplex (must be 0).
 (2) REGULARITY of the quotient (f-vector test): the vertex-canonicalised image complex Q equals the
     orbit space T^4/Z2 iff f_j(Q) = f_j(T)/2 for j>=1 and f_0(Q) = (f_0(T)+#fixed)/2.
 (3) OPEN-STAR disjointness (v3 fix; v2 used closed stars, stricter than Mayer-Vietoris needs):
     open stars of p != q are disjoint iff no simplex of Q contains both p and q (p,q non-adjacent).
     Closed-star disjointness is reported too, for comparison with v2.
 (4) link of EVERY singular point: Betti numbers over Z/2, Z/3, Z/5 (RP^3 expected: (1,1,1,1) mod 2, (1,0,0,1) odd).
 (5) link of an ordinary vertex: S^3 Betti (1,0,0,1) over Z/2, Z/3, Z/5.
 (6) chi(Q) from the f-vector vs (chi(T^4)+#fixed)/2 (orbifold formula, both computed).
"""
import json, sys, time
from lib_common import *

Ns = [int(a) for a in sys.argv[1:]] or [4, 6, 8]
FIELDS = (2, 3, 5)
OUT = {"args": Ns}


def link_of(q_tops, v):
    return [tuple(sorted(s - {v})) for s in q_tops if v in s]


for N in Ns:
    t0 = time.time()
    tops, q, fixed, collide = build_quotient(N)
    fT = fvector(all_faces_by_dim(tops))
    faces_q = all_faces_by_dim([tuple(s) for s in q])
    fQ = fvector(faces_q)
    nfix = len(fixed)
    reg = (fQ[0] == (fT[0] + nfix) // 2 and (fT[0] + nfix) % 2 == 0 and
           all(fQ[j] * 2 == fT[j] for j in range(1, len(fT))))
    chiT, chiQ = chi_from_f(fT), chi_from_f(fQ)
    # adjacency / stars
    fset = set(fixed)
    edges = faces_q[1]
    adjacent_pairs = [sorted(e) for e in edges if e <= fset]
    stars = {p: set().union(*[s for s in q if p in s]) for p in fixed}
    closed_viol = sum(1 for i in range(nfix) for j in range(i + 1, nfix) if stars[fixed[i]] & stars[fixed[j]])
    links = {}
    for p in fixed:
        lt = link_of(q, p)
        links[p] = {f"Z{f}": pad(gudhi_betti(lt, f), 4) for f in FIELDS}
    agree = {f"Z{f}": len({tuple(l[f"Z{f}"]) for l in links.values()}) == 1 for f in FIELDS}
    rep = links[fixed[0]]
    ordinary = next(v for s in q for v in sorted(s) if v not in fset)
    lo = {f"Z{f}": pad(gudhi_betti(link_of(q, ordinary), f), 4) for f in FIELDS}
    OUT[f"N={N}"] = {
        "N": N, "n_fixed_points": nfix, "orbit_collisions_in_simplices": collide,
        "fvector_T4": fT, "fvector_Q": fQ, "chi_T4_from_f": chiT, "chi_Q_from_f": chiQ,
        "chi_Q_orbifold_formula_(chiT+n_fixed)/2": (chiT + nfix) // 2,
        "chi_Q_matches_orbifold_formula": chiQ == (chiT + nfix) // 2,
        "quotient_is_regular_f_halving": reg,
        "n_adjacent_singular_pairs(open_stars_meet)": len(adjacent_pairs),
        "open_stars_pairwise_disjoint": len(adjacent_pairs) == 0,
        "n_closed_star_violating_pairs": closed_viol,
        "closed_stars_pairwise_disjoint(v2_criterion)": closed_viol == 0,
        "singular_link_betti_by_field_all_points": {str(p): l for p, l in links.items()},
        "all_singular_links_agree_per_field": agree,
        "singular_link_representative": rep,
        "singular_link_is_RP3_homology": rep["Z2"] == [1, 1, 1, 1] and rep["Z3"] == [1, 0, 0, 1] and rep["Z5"] == [1, 0, 0, 1] and all(agree.values()),
        "ordinary_vertex_link": lo,
        "ordinary_link_is_S3": all(lo[k] == [1, 0, 0, 1] for k in lo),
        "runtime_sec": round(time.time() - t0, 1),
    }
    o = OUT[f"N={N}"]
    o["PREMISE_VALID_FOR_MV"] = bool(o["quotient_is_regular_f_halving"] and o["open_stars_pairwise_disjoint"]
                                     and o["singular_link_is_RP3_homology"] and o["ordinary_link_is_S3"]
                                     and o["orbit_collisions_in_simplices"] == 0)
    print(N, {k: o[k] for k in ("n_fixed_points", "quotient_is_regular_f_halving", "open_stars_pairwise_disjoint",
                                "closed_stars_pairwise_disjoint(v2_criterion)", "singular_link_is_RP3_homology",
                                "ordinary_link_is_S3", "chi_Q_from_f", "PREMISE_VALID_FOR_MV", "runtime_sec")}, flush=True)
    (HERE / "00_premise_results.json").write_text(json.dumps(OUT, indent=1))
