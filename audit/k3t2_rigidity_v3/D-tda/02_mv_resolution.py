#!/usr/bin/env python
"""Track D v3, step 2: Mayer-Vietoris resolution of the singular points of T^4/Z2, b0..b4 all COMPUTED.

Run:  cd audit/k3t2_rigidity_v3/D-tda && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python 02_mv_resolution.py 4 6 8 --chain-max-N 8
(positional = grid sizes N, only those with PREMISE_VALID_FOR_MV in 00_premise_results.json are used;
 --chain-max-N M = largest N for which the chain-level rank of the H_3 map is computed; committed run: 4 6 8 --chain-max-N 8)

Decomposition (open-star criterion, v3):
  A = Q minus open stars of the singular set P  (deformation retract: FULL subcomplex U of simplices avoiding P;
      v2 used only the tops avoiding P, which drops the link simplices)
  B = disjoint union of k pieces (cone, or disc bundle over S^2), A cap B ~ k copies of the link ~ RP^3.
Field F = Z/3 and Z/5 (odd; H_1 = H_2 = 0 for RP^3, so phi_1, phi_2, phi_4 have zero domain).
Ranks needed: phi_0 (component bookkeeping, exact) and phi_3 : H_3(k links) -> H_3(U;F).
phi_3 is computed at CHAIN LEVEL: fundamental cycle of each link (nullspace of the link boundary map mod p),
rank of {[c_p]} in H_3(U) = rank([d_4 | c_1..c_k]) - rank(d_4) by sparse elimination mod p.
Cross-check (identity, not an input): rank phi_3 = b3(U) - b3(Q) from GUDHI.
Validation: MV with all B = cone reproduces the GUDHI Betti numbers of Q for every field.
Local model of a resolved point (disc bundle over S^2 with boundary the link) is a DECLARED INPUT (tier L).
"""
import json, sys, time
import itertools
from lib_common import *

args = sys.argv[1:]
chain_max = 6
if "--chain-max-N" in args:
    i = args.index("--chain-max-N"); chain_max = int(args[i + 1]); del args[i:i + 2]
Ns = [int(a) for a in args] or [4, 6, 8]
prem = json.loads((HERE / "00_premise_results.json").read_text())
FIELDS = (3, 5)
OUT = {"args": sys.argv[1:], "chain_max_N": chain_max}
S2b = pad([1, 0, 1], 5)
firstN_exported = False

for N in Ns:
    if not prem.get(f"N={N}", {}).get("PREMISE_VALID_FOR_MV"):
        OUT[f"N={N}"] = {"skipped": "premise not valid / not computed"}; continue
    t0 = time.time()
    tops, qset, fixed, _ = build_quotient(N)
    q = [tuple(sorted(s)) for s in qset]
    P = set(fixed); k = len(P)
    order = {v: i for i, v in enumerate(sorted({v for s in q for v in s}))}
    fq = fvector(all_faces_by_dim(q))
    U_facets = [tuple(v for v in s if v not in P) for s in q]
    U_facets = [s for s in U_facets if s]
    fU = fvector(all_faces_by_dim(U_facets))
    entry = {"N": N, "n_singular": k, "chi_Q_from_f": chi_from_f(fq), "chi_U_from_f": chi_from_f(fU), "fvector_Q": fq, "fvector_U": fU}
    # cone (closed star of one singular vertex) and link
    p0 = fixed[0]
    star = [s for s in q if p0 in s]
    fcone = fvector(all_faces_by_dim(star))
    link_tops = {p: [tuple(sorted(set(s) - {p})) for s in q if p in s] for p in fixed}
    fL = fvector(all_faces_by_dim(link_tops[p0]))
    entry.update(chi_cone_from_f=chi_from_f(fcone), chi_link_from_f=chi_from_f(fL))
    entry["betti_Q"] = {f"Z{f}": pad(gudhi_betti(q, f)) for f in (2, 3, 5)}
    entry["betti_U"] = {f"Z{f}": pad(gudhi_betti(U_facets, f)) for f in (2, 3, 5)}
    entry["betti_cone_Z3"] = pad(gudhi_betti(star, 3))
    # components of U (H_0)
    parent = {}
    def find(x):
        while parent.setdefault(x, x) != x:
            parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for s in U_facets:
        for v in s[1:]:
            parent[find(v)] = find(s[0])
    comps_U = len({find(v) for s in U_facets for v in s})
    link_comp = {p: {find(v) for s in link_tops[p] for v in s} for p in fixed}
    assert all(len(c) == 1 for c in link_comp.values())
    entry["components_U"] = comps_U
    # chain-level H_3 map
    chain_ok = N <= chain_max
    rank3 = {}
    if chain_ok:
        # tetrahedra of U: 4-subsets of tops avoiding P
        tetra_idx = {}
        for s in q:
            so = tuple(sorted(s, key=order.get))
            for c in itertools.combinations(so, 4):
                if not (set(c) & P):
                    tetra_idx.setdefault(c, len(tetra_idx))
        cols4 = []
        for s in q:
            if set(s) & P: continue
            so = tuple(sorted(s, key=order.get))
            cols4.append({tetra_idx[so[:i] + so[i + 1:]]: (-1) ** i for i in range(5)})
        # link fundamental cycles per field
        cycles = {}
        for p in fixed:
            lt = [tuple(sorted(t, key=order.get)) for t in link_tops[p]]
            tri_idx = {}
            bcols = []
            for t in lt:
                col = {}
                for i in range(4):
                    f = t[:i] + t[i + 1:]
                    col[tri_idx.setdefault(f, len(tri_idx))] = (-1) ** i
                bcols.append(col)
            cycles[p] = (lt, bcols)
        for F in FIELDS:
            sr = SparseRank(F)
            for c in cols4: sr.add(c)
            r4 = sr.rank
            dims_null = {}
            added = 0
            for p in fixed:
                lt, bcols = cycles[p]
                ns = nullspace_gfp(bcols, F)
                dims_null[p] = len(ns)
                assert len(ns) == 1, "link fundamental class not 1-dim"
                cyc = {tetra_idx[lt[i]]: v for i, v in ns[0].items()}
                added += sr.add(cyc)
            rank3[F] = added
            entry.setdefault("chain_level", {})[f"Z{F}"] = {"rank_d4_U": r4, "link_cycle_space_dims_all_1": all(v == 1 for v in dims_null.values()),
                                                           "rank_phi3": added}
    # MV
    for F in FIELDS:
        fn = f"Z{F}"
        bU, bQ = entry["betti_U"][fn], entry["betti_Q"][fn]
        bL = pad(gudhi_betti(link_tops[p0], F), 4)
        assert all(pad(gudhi_betti(link_tops[p], F), 4) == bL for p in fixed)
        r3_id = bU[3] - bQ[3]           # identity from the singular MV sequence (cone pieces have H_3 = 0)
        r3 = rank3.get(F, r3_id)
        # phi_0: k link components -> U component + own B component
        sr0 = SparseRank(F)   # rows: link components; columns: U components (by id) and one B component per link
        ucomp_ids = {c: i for i, c in enumerate(sorted({find(v) for s in U_facets for v in s}))}
        for i, p in enumerate(fixed):
            sr0.add({ucomp_ids[next(iter(link_comp[p]))]: 1, comps_U + i: 1})
        r0 = sr0.rank
        ranks = {0: r0, 3: r3}
        bAB = [k * x for x in bL]
        pieces = {}
        for j in range(k + 1):
            bB = [(k - j) * c + j * s for c, s in zip(pad(entry["betti_cone_Z3"]), S2b)]
            m = mv_betti(bU, bB, bAB, ranks)
            b = [x["b_n"] for x in m]
            pieces[j] = {"betti": b, "chi_from_betti": sum((-1) ** i * x for i, x in enumerate(b))}
        m0 = [x["b_n"] for x in mv_betti(bU, [k * c for c in pad(entry["betti_cone_Z3"])], bAB, ranks)]
        m_full = pieces[k]
        b = m_full["betti"]
        entry[fn] = {"betti_link": bL, "rank_phi0": r0, "rank_phi3_used": r3, "rank_phi3_source": "chain level" if F in rank3 else "GUDHI identity b3(U)-b3(Q) (chain level skipped)",
                     "rank_phi3_identity_b3U_minus_b3Q": r3_id, "chain_equals_identity": (F not in rank3) or rank3[F] == r3_id,
                     "MV_singular_reproduces_GUDHI_Q": m0 == bQ, "MV_singular_betti": m0, "GUDHI_betti_Q": bQ,
                     "resolved_betti_b0_b4": b, "chi_resolved_from_betti": m_full["chi_from_betti"],
                                          "scan_over_j_resolved": {str(j): v for j, v in pieces.items()}}
    S2f = chi_from_f(fvector(all_faces_by_dim(list(itertools.combinations(range(4), 3)))))
    for F in FIELDS:
        e = entry[f"Z{F}"]
        e["chi_S2_used_in_fvector_formula"] = S2f
        e["chi_resolved_from_fvectors"] = entry["chi_Q_from_f"] + k * (S2f - entry["chi_cone_from_f"])
        e["chi_betti_equals_fvector"] = e["chi_resolved_from_betti"] == e["chi_resolved_from_fvectors"]
    entry["runtime_sec"] = round(time.time() - t0, 1)
    OUT[f"N={N}"] = entry
    e3 = entry["Z3"]
    print(N, "k", k, "b", e3["resolved_betti_b0_b4"], "chi", e3["chi_resolved_from_betti"], e3["chi_resolved_from_fvectors"],
          "chain==id", e3["chain_equals_identity"], "MVsing", e3["MV_singular_reproduces_GUDHI_Q"], entry["runtime_sec"], flush=True)
    (HERE / "02_mv_results.json").write_text(json.dumps(OUT, indent=1, default=str))
    if not firstN_exported:
        ex = {"source": f"02_mv_resolution.py, N={N}, Z/3", "chi": e3["chi_resolved_from_betti"], "b2": e3["resolved_betti_b0_b4"][2],
              "n_singular": k, "betti_b0_b4": e3["resolved_betti_b0_b4"], "status": "computed by exact Mayer-Vietoris on GUDHI Betti numbers; depends on declared inputs (local model of a resolved point)"}
        (HERE / "exports.json").write_text(json.dumps(ex, indent=1)); firstN_exported = True
