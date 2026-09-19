"""Track C v3 part 1: E8 Gram from root system, candidate lattices m U + n(-E8), controls.
cd audit/k3t2_rigidity_v3/C-lattices && PY 01_lattices.py
Writes 01_lattices.json. All arithmetic exact (Fraction/sympy)."""
import json, itertools
from fractions import Fraction
from pathlib import Path
import sympy as sp
from decl import HERE

def e8_roots():
    R = []
    for i, j in itertools.combinations(range(8), 2):
        for si in (1, -1):
            for sj in (1, -1):
                v = [Fraction(0)] * 8; v[i] = Fraction(si); v[j] = Fraction(sj); R.append(tuple(v))
    for s in itertools.product((1, -1), repeat=8):
        if s.count(-1) % 2 == 0:
            R.append(tuple(Fraction(x, 2) for x in s))
    return R

def dot(a, b): return sum(x * y for x, y in zip(a, b))

def simple_roots(R):
    f = [Fraction(10**(8 - k)) + Fraction(1, 3 + k) for k in range(8)]  # generic functional
    pos = [r for r in R if dot(f, r) > 0]
    ps = set(pos)
    simple = [r for r in pos if not any(tuple(a - b for a, b in zip(r, q)) in ps for q in pos)]
    return pos, simple

def gram(vs): return sp.Matrix([[sp.Rational(dot(a, b).numerator, dot(a, b).denominator) for b in vs] for a in vs])

def blockdiag(*Ms):
    n = sum(M.rows for M in Ms); G = sp.zeros(n, n); k = 0
    for M in Ms: G[k:k+M.rows, k:k+M.rows] = M; k += M.rows
    return G

def signature(G):
    """exact: number of positive / negative roots of charpoly via sympy exact real-root counting"""
    x = sp.symbols('x'); p = sp.Poly(G.charpoly(x).as_expr(), x)
    npos = p.count_roots(sp.Rational(1, 10**9), None) if False else None
    # count_roots counts with multiplicity? use sqf-part-free approach: factor
    pos = neg = zero = 0
    for fac, mult in sp.factor_list(p.as_expr())[1]:
        q = sp.Poly(fac, x)
        n_ge0 = q.count_roots(0, None); n_le0 = q.count_roots(None, 0); n_all = q.degree()
        z = n_ge0 + n_le0 - n_all  # roots exactly 0 (simple in sqfree factor)
        pos += mult * (n_ge0 - z); neg += mult * (n_le0 - z); zero += mult * z
    return pos, neg, zero

def inertia_ldl(G):
    """independent route: exact symmetric congruence diagonalisation over Q"""
    n = G.rows; A = sp.Matrix(G); pos = neg = zero = 0
    active = list(range(n))
    while active:
        piv = next((i for i in active if A[i, i] != 0), None)
        if piv is None:
            pr = next(((i, j) for i in active for j in active if i < j and A[i, j] != 0), None)
            if pr is None: zero += len(active); break
            i, j = pr; A[i, :] = A[i, :] + A[j, :]; A[:, i] = A[:, i] + A[:, j]; continue
        p = A[piv, piv]
        for i in active:
            if i != piv:
                c = A[i, piv] / p; A[i, :] = A[i, :] - c * A[piv, :]; A[:, i] = A[:, i] - c * A[:, piv]
        if p > 0: pos += 1
        else: neg += 1
        active.remove(piv)
    return pos, neg, zero

res = {}
R = e8_roots(); pos, simple = simple_roots(R)
E8 = gram(simple)
res["e8"] = {"n_roots": len(R), "n_positive": len(pos), "n_simple": len(simple),
  "gram": E8.tolist().__repr__(), "det": int(E8.det()), "diag": [int(E8[i, i]) for i in range(8)],
  "even": all(E8[i, i] % 2 == 0 for i in range(8)), "signature_ldl": inertia_ldl(E8),
  "expected_from_literature": {"det": 1, "n_roots": 240, "source": "FROM MEMORY"}}
# negative controls on E8-ness: Z^8 (odd), D8 (det 4)
Z8 = sp.eye(8); res["controls_rank8"] = {"Z8": {"det": 1, "even": False},
  "D8_from_roots_+-ei+-ej": {}}
D8r = [tuple(Fraction(x) for x in [1 if k == i else (-1 if k == i + 1 else 0) for k in range(8)]) for i in range(7)] + \
      [tuple(Fraction(1 if k in (6, 7) else 0) for k in range(8))]
D8 = gram(D8r); res["controls_rank8"]["D8_from_roots_+-ei+-ej"] = {"det": int(D8.det()), "even": True}

U = sp.Matrix(((0, 1), (1, 0))); mE8 = -E8
def build(m, n): return blockdiag(*([U] * m + [mE8] * n)) if (m + n) else sp.zeros(0, 0)

# chain: chi_top -> signature (see 02); here candidate family scan, EVERY (m,n) with 2m+8n=rank, rank from 02
# rank is computed in 02; to avoid circular import we recompute the chain deterministically here
from k3chain import chain
ch = chain()
rank = ch["b2"]; sig_target = (ch["b2plus"], ch["b2minus"])
cands = []
for m in range(0, rank // 2 + 1):
    for n in range(0, rank // 8 + 1):
        if 2 * m + 8 * n != rank: continue
        G = build(m, n)
        s = inertia_ldl(G); s2 = signature(G)
        cands.append({"m": m, "n": n, "rank": G.rows, "det": int(G.det()), "even": all(G[i, i] % 2 == 0 for i in range(G.rows)),
                      "signature_ldl": s[:2], "signature_charpoly": s2[:2], "routes_agree": s[:2] == s2[:2],
                      "matches_hodge_signature": s[:2] == sig_target})
res["candidates_mU_nE8_rank_b2"] = cands
res["selected"] = [c for c in cands if c["matches_hodge_signature"]]
# also all (m,n) 0<=m,n<=30 with any rank whose signature == target, unique?
allmn = []
for m in range(0, 31):
    for n in range(0, 31):
        if (m, m + 8 * n) == sig_target: allmn.append((m, n))
res["all_mn_0_30_with_signature_equal_hodge"] = allmn
# note: mU+n(-E8) has signature (m, m+8n); realizable objects: all are even unimodular lattices that exist
G322 = build(3, 2)
res["typed_gram_3U_2mE8"] = {"rank": G322.rows, "det": int(G322.det()), "signature_ldl": inertia_ldl(G322)[:2],
   "signature_charpoly": signature(G322)[:2],
   "note": "typed input (decl.py typed_gram_3U_2mE8): rank/signature are consequences of the typed matrix, NOT independent of it"}
# Mukai lattice: H^0+H^2+H^4 = U (+) H^2 -> 4U+2(-E8): rank = sum of Betti numbers
G24 = build(4, 2)
res["mukai"] = {"rank": G24.rows, "signature": inertia_ldl(G24)[:2], "rank_equals_chi_top": G24.rows == ch["chi_top"],
   "note": "rank 24 = sum b_i is chi_top by chi_top = sum(-1)^i b_i with odd Betti zero; the Gram is typed"}
(HERE / "01_lattices.json").write_text(json.dumps(res, indent=2, default=str))
print(json.dumps({k: res[k] for k in ("selected", "all_mn_0_30_with_signature_equal_hodge", "typed_gram_3U_2mE8", "mukai")}, indent=1, default=str))
print(res["e8"]["det"], res["e8"]["n_roots"], res["e8"]["n_simple"], res["e8"]["even"], res["e8"]["signature_ldl"])
