"""
lib_kummer.py -- chain-level Kummer construction on the cubical orbifold T^4/Z2.

  X   = T^4/Z2 as the orbit CW complex of the cubical Z_N^4 complex (N even, N >= 6
        so that the closed stars of the 16 singular vertices are pairwise disjoint).
  U   = cells whose closure contains no singular vertex (subcomplex; computed).
  L_p = cells of the closed star box [p-1,p+1]^4 not containing p (subcomplex of U;
        computed; checked to be RP^3 over F2/F3/F5).
  Replacement pieces glued along L_p (algebraic mapping cylinder of phi: C(L_p) -> C(B)):
    'cone'     : B = point,                  phi_0 = augmentation            -> closed star (cone on RP^3)
    'resolve'  : B = S^2 = e0 u e2,           phi_0 = augmentation,
                                              phi_2 = delta(x~)/2, x = generator of H^1(L_p;F2)
                                              (Bockstein: the generator of H^2(RP^3;Z) = Z/2)
                                              -> disc bundle D(O(-2)) = Cyl(RP^3 -> S^2)   [tier-L input:
                                                 the circle-bundle projection pulls [S^2] back to that generator]
    'trivial2' : B = S^2, phi_2 = 0            -> WRONG gluing class (negative control)

Everything is computed from the cells; the only non-computed input is the choice of the
cohomology class of phi_2 (tier L, Gysin with Euler number +-2), and its consequence
H_*(Cyl, L) = (0,0,1,0,1) (Lefschetz duality) is CHECKED at chain level.
"""
from __future__ import annotations

import itertools

from lib_cells import (CC, cubical_torus, quotient_by_involution, neg_action, subcomplex,
                       glue_cylinders, relative_complex, homology_mod_p_colred, nullspace_mod_p,
                       dense_rank_mod_p)


def closure_vertices(v, S, N):
    out = []
    for T in itertools.product(*[(0, 1) if i in S else (0,) for i in range(len(v))]):
        out.append(tuple((x + t) % N for x, t in zip(v, T)))
    return out


def cyc_dist(a, b, N):
    d = (a - b) % N
    return min(d, N - d)


def build_orbifold_pieces(N):
    assert N % 2 == 0 and N >= 6, "need N even >= 6 for disjoint closed stars"
    cc, index = cubical_torus(N, 4)
    q, info, _, _ = quotient_by_involution(cc, index, neg_action(N))
    fixed = [v for v, S in info["fixed_cells"]]
    assert len(fixed) == 16
    contains = []
    for o, (v, S) in enumerate(q.label):
        cv = set(closure_vertices(v, S, N))
        contains.append([p for p in fixed if p in cv])
    U = [o for o in range(len(q)) if not contains[o]]
    Uset = set(U)
    links = {}
    for p in fixed:
        Lp = []
        for o in U:
            v, S = q.label[o]
            if all(all(cyc_dist(w[i], p[i], N) <= 1 for i in range(4)) for w in closure_vertices(v, S, N)):
                Lp.append(o)
        links[p] = Lp
    # checks (reported, not repaired)
    checks = {}
    allL = [o for p in fixed for o in links[p]]
    checks["links_pairwise_disjoint"] = len(allL) == len(set(allL))
    star_cells = [o for o in range(len(q)) if contains[o]]
    checks["each_star_cell_contains_exactly_one_singular_vertex"] = all(len(contains[o]) == 1 for o in star_cells)
    # every face of a star cell is either a star cell of the same p or in L_p
    ok = True
    for o in star_cells:
        p = contains[o][0]
        Lset = set(links[p])
        for f in q.bd[o]:
            if not (contains[f] == [p] or f in Lset):
                ok = False
    checks["star_faces_in_star_or_link"] = ok
    # U closed under faces
    checks["U_closed_under_faces"] = all(f in Uset for o in U for f in q.bd[o])
    checks["link_closed_under_faces"] = all(f in set(links[p]) for p in fixed for o in links[p] for f in q.bd[o])
    sizes = {"X_fvector": q.fvector(), "U_n": len(U), "star_n_total": len(star_cells),
             "L_fvector": subcomplex(q, links[fixed[0]])[0].fvector()}
    return q, fixed, U, links, checks, sizes


def bockstein_phi2(q, Lp):
    """
    Integral 2-cocycle on L_p representing the generator of H^2(L_p;Z)=Z/2:
    x = an F2 1-cocycle not an F2 coboundary; phi2 = delta(x~)/2 with x~ in {0,1}.
    Returns (phi2: dict cell_id -> int, diagnostics).
    """
    L, new = subcomplex(q, Lp)
    old = {v: k for k, v in new.items()}
    c0 = [i for i in range(len(L)) if L.dims[i] == 0]
    c1 = [i for i in range(len(L)) if L.dims[i] == 1]
    c2 = [i for i in range(len(L)) if L.dims[i] == 2]
    c3 = [i for i in range(len(L)) if L.dims[i] == 3]
    i1 = {c: j for j, c in enumerate(c1)}
    i2 = {c: j for j, c in enumerate(c2)}
    # delta_1 : C^1 -> C^2, (delta x)(c) = sum_e [dc:e] x(e); rows = 2-cells
    D1 = [[0] * len(c1) for _ in c2]
    for c in c2:
        for e, a in L.bd[c].items():
            D1[i2[c]][i1[e]] = a
    # delta_0 image: columns = coboundaries of vertices, as vectors on C^1
    D0cols = []
    for v in c0:
        col = [0] * len(c1)
        for e in c1:
            if v in L.bd[e]:
                col[i1[e]] = L.bd[e][v]
        D0cols.append(col)
    r0 = dense_rank_mod_p(D0cols, 2)
    ker = nullspace_mod_p(D1, len(c1), 2)
    x = None
    for z in ker:
        if dense_rank_mod_p(D0cols + [z], 2) > r0:
            x = z
            break
    diag = {"dim_ker_delta1_F2": len(ker), "rank_delta0_F2": r0, "dim_H1_F2": len(ker) - r0}
    assert x is not None, "no nontrivial F2 1-cocycle"
    phi2 = {}
    for c in c2:
        s = sum(a * x[i1[e]] for e, a in L.bd[c].items())
        assert s % 2 == 0, "delta(x~) not even: x not an F2 cocycle"
        if s // 2:
            phi2[old[c]] = s // 2
    # integral cocycle check: (delta phi2)(d) = sum_c [dd:c] phi2(c) = 0 for 3-cells d
    viol = 0
    for d in c3:
        s = sum(a * phi2.get(old[c], 0) for c, a in L.bd[d].items())
        viol += (s != 0)
    diag["integral_cocycle_violations"] = viol
    # nontriviality mod 2: phi2 mod 2 not in image of delta_1 (mod 2)
    D1cols = [[D1[r][j] for r in range(len(c2))] for j in range(len(c1))]
    vec = [phi2.get(old[c], 0) for c in c2]
    r1 = dense_rank_mod_p(D1cols, 2)
    diag["phi2_mod2_not_coboundary"] = dense_rank_mod_p(D1cols + [vec], 2) > r1
    diag["phi2_support"] = len(phi2)
    diag["phi2_values"] = sorted(set(phi2.values()))
    return phi2, diag


def piece(kind, Lp, phi2=None):
    if kind == "cone":
        return {"L": Lp, "base": [(0, {})], "phi": lambda c, k: {0: 1} if k == 0 else {}}
    if kind == "resolve":
        return {"L": Lp, "base": [(0, {}), (2, {})],
                "phi": lambda c, k: ({0: 1} if k == 0 else ({1: phi2[c]} if (k == 2 and c in phi2) else {}))}
    if kind == "trivial2":
        return {"L": Lp, "base": [(0, {}), (2, {})], "phi": lambda c, k: {0: 1} if k == 0 else {}}
    raise ValueError(kind)


def assemble(q, fixed, U, links, kinds, phi2s):
    """kinds: list of 16 strings in the order of `fixed`."""
    pieces = [piece(k, links[p], phi2s.get(p)) for k, p in zip(kinds, fixed)]
    return glue_cylinders(q, U, pieces)


def local_model_relative(q, Lp, kind, phi2=None):
    """Chain complex of (Cyl(phi), L_p) relative: Cyl built on L_p alone, then L_p cells dropped."""
    Z, _ = glue_cylinders(q, Lp, [piece(kind, Lp, phi2)])
    drop = [i for i, lab in enumerate(Z.label) if lab[0] == "Y"]
    return relative_complex(Z, drop), Z
