#!/usr/bin/env python
"""Track D v3, step 3: chain-level cup pairing on H^2(-;Z/2) (attempt at intersection-form parity).

Run:  cd audit/k3t2_rigidity_v3/D-tda && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python 03_cup_parity.py 4
(argument = grid size N; committed run used 4)

Ordered simplicial cup product (Alexander-Whitney) on the Kuhn complex T^4 and on the quotient Q = T^4/Z2 (singular,
NOT resolved). Pairing b(x,y) = sum over 4-simplices (v0..v4) x(v0v1v2) y(v2v3v4) mod 2, evaluated on the mod-2
fundamental class [sum of all 4-simplices] (checked to be a mod-2 cycle). q(x)=b(x,x) is additive on H^2(;Z/2);
'even' means q == 0 on all of H^2.
SCOPE: this does NOT decide the parity of the K3 form. The resolved K3 has no triangulation here, and H^2(Q;Z/2)
is not the K3 lattice mod 2. No odd-form control is available (a control complex with odd form, e.g. CP^2, was not built),
so a result 'q == 0' is only shown on cases expected to be even; nonzero q on Q shows the code can output nonzero.
"""
import json, sys, itertools, random
from lib_common import *

N = int(sys.argv[1]) if len(sys.argv) > 1 else 4


def gf2_rows_rank_basis(vecs):
    """vecs: list of python-int bitsets. returns (basis dict pivot->vec)"""
    piv = {}
    for v in vecs:
        while v:
            h = v.bit_length() - 1
            if h in piv: v ^= piv[h]
            else: piv[h] = v; break
    return piv


def reduce(v, piv):
    while v:
        h = v.bit_length() - 1
        if h in piv: v ^= piv[h]
        else: return v
    return 0


def analyse(tops_sets, label):
    tops = [tuple(sorted(s)) for s in tops_sets]
    allv = sorted({v for s in tops for v in s}); order = {v: i for i, v in enumerate(allv)}
    tops = [tuple(sorted(s, key=order.get)) for s in tops]
    faces = {d: {} for d in range(5)}
    for s in tops:
        for d in range(5):
            for c in itertools.combinations(s, d + 1):
                faces[d].setdefault(c, len(faces[d]))
    n = {d: len(faces[d]) for d in faces}
    # mod-2 cycle check of fundamental class
    cnt = {}
    for s in tops:
        for c in itertools.combinations(s, 4): cnt[c] = cnt.get(c, 0) + 1
    fund_is_cycle = all(v % 2 == 0 for v in cnt.values())
    # coboundaries as bitsets over C^2 (triangle indexing)
    # delta1: C^1 -> C^2, e -> sum of triangles containing e ; delta2: C^2 -> C^3 dual
    # Cocycles Z^2 = {x in C^2 : for every tetra t, sum_{tri in t} x(tri) = 0}
    tet_rows = []
    for t in faces[3]:
        r = 0
        for c in itertools.combinations(t, 3): r |= 1 << faces[2][c]
        tet_rows.append(r)
    # nullspace of tet_rows system (unknowns: triangles)
    piv = {}  # pivot bit -> (row, ) for elimination; solve by reduced echelon
    rows = []
    for r in tet_rows:
        while r:
            h = r.bit_length() - 1
            if h in piv: r ^= piv[h]
            else: piv[h] = r; break
    # reduce to RREF
    hs = sorted(piv)
    for h in hs:
        for h2 in hs:
            if h2 != h and (piv[h2] >> h) & 1: piv[h2] ^= piv[h]
    free = [i for i in range(n[2]) if i not in piv]
    Z = []
    for f in free:
        x = 1 << f
        for h, r in piv.items():
            if (r >> f) & 1: x |= 1 << h
        Z.append(x)
    # coboundaries B^2: image of delta1 (each edge -> triangles containing it)
    B = []
    edge_tris = {e: 0 for e in faces[1]}
    for t, i in faces[2].items():
        for e in itertools.combinations(t, 2): edge_tris[e] |= 1 << i
    B = list(edge_tris.values())
    pB = gf2_rows_rank_basis(B)
    reps = []
    pfull = dict(pB)
    for z in Z:
        r = reduce(z, pfull)
        if r:
            reps.append(z); pfull[r.bit_length() - 1] = r
    h2 = len(reps)
    # cup evaluation
    def bil(x, y):
        s = 0
        for t in tops:
            a = faces[2][(t[0], t[1], t[2])]; b = faces[2][(t[2], t[3], t[4])]
            s ^= ((x >> a) & 1) & ((y >> b) & 1)
        return s
    G = [[bil(x, y) for y in reps] for x in reps]
    sym = all(G[i][j] == G[j][i] for i in range(h2) for j in range(h2))
    # rank of G over GF(2)
    rk = len(gf2_rows_rank_basis([sum(G[i][j] << j for j in range(h2)) for i in range(h2)]))
    qv = [G[i][i] for i in range(h2)]
    # sanity: q(coboundary + cocycle) = q(cocycle); additivity on random combos
    rnd = random.Random(0)
    def comb():
        x = 0
        for z in reps:
            if rnd.random() < .5: x ^= z
        return x
    inv_ok = True
    for _ in range(20):
        x = comb(); c = rnd.choice(B)
        inv_ok &= (bil(x ^ c, x ^ c) == bil(x, x))
    return {"label": label, "n_simplices_by_dim": [n[d] for d in range(5)], "mod2_fundamental_class_is_cycle": fund_is_cycle,
            "dim_H2_Z2": h2, "gram_symmetric": sym, "gram_rank_Z2": rk, "q_on_basis": qv, "q_identically_zero": not any(qv),
            "q_invariant_under_adding_coboundary_(20 random checks)": inv_ok}


OUT = {"N": N}
tops, qset, fixed, _ = build_quotient(N)
OUT["T4"] = analyse([frozenset(s) for s in tops], "Kuhn T^4, expected even (spin) and nondegenerate rank 6 (Poincare duality)")
OUT["Q"] = analyse(list(qset), "quotient T^4/Z2, singular orbifold, NOT the K3")
OUT["scope"] = "does not decide the K3 parity; see docstring and could_not_do"
(HERE / "03_cup_parity_results.json").write_text(json.dumps(OUT, indent=1))
print(json.dumps(OUT, indent=1))
