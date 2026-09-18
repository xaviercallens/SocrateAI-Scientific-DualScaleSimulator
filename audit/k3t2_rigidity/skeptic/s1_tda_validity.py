#!/usr/bin/env python
"""
SKEPTIC audit, Track D (TDA): independent re-implementation (does NOT import
lib_freudenthal or any Track D code) of the Freudenthal T^4 triangulation and
the T^4/Z_2 canonical-collapse quotient, followed by validity checks:

  V1  every simplex of the GUDHI SimplexTree has all its codim-1 faces present
      (walk st.get_simplices(), not the construction list)
  V2  SimplexTree simplex count == independently computed face-set size
  V3  quotient vertex count == orbit count 16 + (N^4 - 16)/2 (orbit count from
      the group action, not from the complex)
  V4  chi from f-vector == alternating Betti sum (GUDHI, Z/3 and Z/5)
  V5  2-to-1 check: in every dim >= 1 the quotient has exactly half the T^4
      simplices (the check that actually certifies K/G = |K|/G here)
  V6  Mayer-Vietoris premise used in 03_resolution_hybrid.py: are the closed
      stars of the 16 singular vertices pairwise DISJOINT in the quotient, so
      that A n B is 16 disjoint RP^3 links?  Measured at N=4 and N=6, plus the
      Betti numbers of the union of the 16 links (16 disjoint RP^3 over Z/3
      would give (16,0,0,16)).
  V7  negative control of Track D (free translation) at N=4: run V4 and V5 on
      it -- shows that V4 has no teeth (it passes on an invalid complex) and V5
      is the discriminating check.

Only python ints/tuples for combinatorics; GUDHI for homology.
Output: skeptic/s1_tda_validity_results.json
"""
import itertools
import json
from collections import Counter

import gudhi

OUT_PATH = ("/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/"
            "k3t2_rigidity/skeptic/s1_tda_validity_results.json")
D = 4


def kuhn_top(N):
    tops = []
    for v in itertools.product(range(N), repeat=D):
        for perm in itertools.permutations(range(D)):
            cur = list(v)
            vs = [tuple(cur)]
            for k in perm:
                cur[k] = (cur[k] + 1) % N
                vs.append(tuple(cur))
            assert len(set(vs)) == D + 1
            tops.append(tuple(vs))
    return tops


def faces_of(tops):
    F = set()
    for s in tops:
        for r in range(1, len(s) + 1):
            for c in itertools.combinations(s, r):
                F.add(frozenset(c))
    return F


def fvec(F):
    c = Counter(len(f) - 1 for f in F)
    return [c[d] for d in range(max(c) + 1)]


def gudhi_complex(tops):
    verts = sorted({v for s in tops for v in s})
    idx = {v: i for i, v in enumerate(verts)}
    st = gudhi.SimplexTree()
    for s in tops:
        st.insert([idx[v] for v in s], filtration=0.0)
    return st, idx


def betti(st, field):
    st.compute_persistence(homology_coeff_field=field, min_persistence=0,
                           persistence_dim_max=True)
    return st.betti_numbers()


def faces_closed(st):
    S = {tuple(sorted(s)) for s, _ in st.get_simplices()}
    missing = 0
    for s in S:
        if len(s) > 1:
            for f in itertools.combinations(s, len(s) - 1):
                if f not in S:
                    missing += 1
    return missing, len(S)


def quotient(tops, act):
    canon = lambda v: min(v, act(v))
    q = [tuple(canon(v) for v in s) for s in tops]
    q = [s for s in q if len(set(s)) == len(s)]
    return q


def validity_block(tops_T4, qtops, n_orbits_expected):
    FT = faces_of(tops_T4)
    FQ = faces_of(qtops)
    fT, fQ = fvec(FT), fvec(FQ)
    st, _ = gudhi_complex(qtops)
    missing, nS = faces_closed(st)
    chi_f = sum((-1) ** d * c for d, c in enumerate(fQ))
    out = {
        "fvector_T4": fT, "fvector_quotient": fQ,
        "V1_missing_codim1_faces_in_simplextree": missing,
        "V2_simplextree_num_simplices": st.num_simplices(),
        "V2_faceset_size": len(FQ),
        "V2_match": st.num_simplices() == len(FQ) == nS,
        "V3_num_vertices": st.num_vertices(),
        "V3_orbit_count_expected": n_orbits_expected,
        "V3_match": st.num_vertices() == n_orbits_expected,
        "chi_from_fvector": chi_f,
    }
    for p in (3, 5):
        st_p, _ = gudhi_complex(qtops)
        b = betti(st_p, p)
        out[f"betti_Z{p}"] = b
        out[f"V4_chi_betti_eq_fvector_Z{p}"] = sum((-1) ** k * x for k, x in enumerate(b)) == chi_f
    out["V5_dims1to4_exactly_half"] = all(2 * fQ[d] == fT[d] for d in range(1, D + 1))
    return out


RES = {}
for N in (4, 6):
    tops = kuhn_top(N)
    neg = lambda v, N=N: tuple((-x) % N for x in v)
    fixed = sorted({v for s in tops for v in s if neg(v) == v})
    nv = N ** D
    orbits = len(fixed) + (nv - len(fixed)) // 2
    q = quotient(tops, neg)
    blk = validity_block(tops, q, orbits)
    blk["num_fixed_vertices"] = len(fixed)

    # V6: closed stars of the singular vertices in the quotient
    canon = lambda v, N=N: min(v, neg(v))
    fixed_q = [canon(p) for p in fixed]
    star_verts = {}
    for p in fixed_q:
        star_verts[p] = {v for s in q if p in s for v in s}
    overlaps = 0
    shared_examples = []
    for a, b in itertools.combinations(fixed_q, 2):
        inter = star_verts[a] & star_verts[b]
        if inter:
            overlaps += 1
            if len(shared_examples) < 3:
                shared_examples.append([list(a), list(b), [list(x) for x in sorted(inter)][:4]])
    link_union = [tuple(v for v in s if v != p) for p in fixed_q for s in q if p in s]
    st_l, _ = gudhi_complex(link_union)
    blk["V6_pairs_of_singular_closed_stars_sharing_a_vertex"] = overlaps
    blk["V6_examples"] = shared_examples
    blk["V6_union_of_16_links_betti_Z3"] = betti(st_l, 3)
    blk["V6_expected_if_disjoint_RP3s_Z3"] = [len(fixed), 0, 0, len(fixed)]
    # complement U as in 03_resolution_hybrid (top simplices avoiding all
    # singular vertices), recomputed independently
    fs = set(fixed_q)
    U = [s for s in q if not (set(s) & fs)]
    st_u, _ = gudhi_complex(U)
    blk["U_betti_Z3_independent"] = betti(st_u, 3)
    # U n (union of closed stars) -- the actual MV intersection
    Uf = faces_of(U)
    Bf = faces_of([s for s in q if set(s) & fs])
    inter = [tuple(f) for f in (Uf & Bf)]
    st_i = gudhi.SimplexTree()
    allv = sorted({v for f in inter for v in f})
    ix = {v: i for i, v in enumerate(allv)}
    for f in inter:
        st_i.insert([ix[v] for v in f], filtration=0.0)
    blk["V6_actual_MV_intersection_UcapB_betti_Z3"] = betti(st_i, 3)
    RES[f"N={N}"] = blk

    if N == 4:
        # V7: Track D's own negative control (translation by (2,0,0,0))
        tr = lambda v, N=N: ((v[0] + 2) % N,) + v[1:]
        qt = quotient(tops, tr)
        RES["N=4_translation_control"] = validity_block(tops, qt, nv // 2)

with open(OUT_PATH, "w") as f:
    json.dump(RES, f, indent=1, default=str)
print(json.dumps(RES, indent=1, default=str))
