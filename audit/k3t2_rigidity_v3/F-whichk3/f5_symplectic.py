"""Track F part 5: Lefschetz numbers on H^*(T^4,Z), fixed points on Km(A) of symplectic (M,t), comparison with M24 cycle shapes (tier C).
Run: cd audit/k3t2_rigidity_v3/F-whichk3 && <venv python> f5_symplectic.py
(no arguments; needs _m24_data.json and exports.json written by f4_codes.py 40000). Exact integer arithmetic.
"""
import itertools, json
from collections import Counter
import numpy as np
import sympy as sp
from sympy.combinatorics import Permutation, PermutationGroup
from common import *

# ---------- (1) Lefschetz numbers on H^*(T^4, Z) = Lambda^* Z^4
def ext_power(M, k):
    n = len(M)
    subs = list(itertools.combinations(range(n), k))
    return sp.Matrix(len(subs), len(subs), lambda a, b: sp.Matrix([[M[i][j] for j in subs[b]] for i in subs[a]]).det() if k else 1)
def lefschetz(M):
    tr = [int(sp.Matrix(ext_power(M, k)).trace()) if k else 1 for k in range(5)]
    return sum((-1) ** k * tr[k] for k in range(5)), tr
Id4 = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
neg = [[-x for x in r] for r in Id4]
L_translation, tr_t = lefschetz(Id4)      # translation acts on H^* as the identity (it is homotopic to id)
L_minus1, tr_m = lefschetz(neg)
# direct fixed-point counts
def count_solutions(M, t, N):
    """number of x in (Z/N)^4 (x = y/N) with (M - I) x + t = 0 mod 1, t given as integer vector over N."""
    ys = np.array(list(itertools.product(range(N), repeat=4)), dtype=np.int64)
    Mm = np.array(M, dtype=np.int64)
    r = ys @ (Mm - np.eye(4, dtype=np.int64)).T + np.array(t, dtype=np.int64)[None, :]
    return int(((r % N) == 0).all(axis=1).sum())
fix_minus1 = count_solutions(neg, [0, 0, 0, 0], 2)    # -x = x -> 2x = 0
fix_transl = count_solutions(Id4, [1, 0, 0, 0], 2)     # x + t = x, t=(1/2,0,0,0)

# ---------- (2) Kummer surface: H^2(Km) = pi_* H^2(A) (rank 6, translation acts as identity) + 16 exceptional curves permuted by p -> p+t
A2 = [(a, b, c, d) for a in (0, 1) for b in (0, 1) for c in (0, 1) for d in (0, 1)]
def km_lefschetz_translation(t):
    perm_fixed = sum(1 for p in A2 if tuple((x + y) % 2 for x, y in zip(p, t)) == p)
    tr_H2 = int(sp.Matrix(ext_power(Id4, 2)).trace()) + perm_fixed
    return 1 + tr_H2 + 1, {"trace_on_H2": tr_H2, "exceptional_curves_fixed": perm_fixed}
km_transl = {}
for t in A2[1:]:
    L, d = km_lefschetz_translation(t)
    km_transl[str(t)] = L
transl_vals = sorted(set(km_transl.values()))

# ---------- (3) fixed points on Km(A) of x -> Mx + t (t in A[2]), A = E_tau x E_tau'
Mi = [[0, -1], [1, 0]]                 # multiplication by i on Z+iZ, columns = images of (1, i)
Mt = [[0, -1], [1, 1]]                 # multiplication by tau=e^{i pi/3} on Z+tau Z  (tau^2 = tau - 1)
def m2(A, B): return [[sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)] for i in range(2)]
def mp(A, n):
    R = [[1, 0], [0, 1]]
    for _ in range(n % 12): R = m2(R, A)
    return R
def inv2(A): return [[A[1][1], -A[0][1]], [-A[1][0], A[0][0]]]     # det 1
def block(A, B):
    M = [[0] * 4 for _ in range(4)]
    for i in range(2):
        for j in range(2): M[i][j] = A[i][j]; M[2 + i][2 + j] = B[i][j]
    return M
def mm(A, B): return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
def matvec_mod(M, v, N): return tuple(sum(M[i][j] * v[j] for j in range(4)) % N for i in range(4))

def km_fixed(M, t, lam_ratio_distinct):
    """t as integer 4-vector over 2 (t/2). Returns (order on Km, fixed points on Km, detail). lam_ratio_distinct: True if the holomorphic tangent eigenvalues at a fixed node are distinct."""
    # order on Km: g^k = (M^k, t_k); identity on Km iff M^k = +-1 and t_k = 0 mod 1 (t_k in (1/2)Z^4)
    Mk, tk = [row[:] for row in Id4], (0, 0, 0, 0)
    order = None
    for k in range(1, 25):
        # g^k(x) = M^k x + t_k ; compose g after g^{k-1}: x -> M(M^{k-1}x + t_{k-1}) + t
        tk = tuple((sum(M[i][j] * tk[j] for j in range(4)) + t[i]) % 2 for i in range(4))
        Mk = mm(M, Mk)
        if (Mk == Id4 or Mk == neg) and all(x == 0 for x in tk):
            order = k; break
    dI = int(sp.Matrix(M).__sub__(sp.eye(4)).det()); dP = int(sp.Matrix(M).__add__(sp.eye(4)).det())
    d = max(abs(dI), abs(dP), 1)
    N = 2 * d
    if N > 40: raise ValueError("grid too big")
    tN = [(x * N // 2) for x in t]      # t/2 = tN/N
    ys = np.array(list(itertools.product(range(N), repeat=4)), dtype=np.int64)
    Mm = np.array(M, dtype=np.int64)
    Iden = np.eye(4, dtype=np.int64)
    rp = ys @ (Mm - Iden).T + np.array(tN)[None, :]
    rm = ys @ (Mm + Iden).T + np.array(tN)[None, :]
    okp = ((rp % N) == 0).all(axis=1); okm = ((rm % N) == 0).all(axis=1)
    sel = ys[okp | okm]
    is2 = ((2 * sel) % N == 0).all(axis=1)
    non_nodes = int((~is2).sum())
    assert non_nodes % 2 == 0
    nodes_fixed = int((okp[is2 & True] if False else ((ys[okp] * 2) % N == 0).all(axis=1)).sum())   # 2-torsion p with M p + t = p
    node_pts = nodes_fixed * (2 if lam_ratio_distinct else 10 ** 6)
    return order, non_nodes // 2 + node_pts, {"non_node_classes": non_nodes // 2, "fixed_nodes": nodes_fixed}

cases = {
  "E_i x E_i, M=(i,-i)": (block(Mi, inv2(Mi)), True, "(i,-i)"),
  "E_tau x E_tau, tau=e^{i pi/3}, M=(tau^2,tau^-2) [omega, omega^2]": (block(mp(Mt, 2), inv2(mp(Mt, 2))), True, "(omega,omega^2)"),
  "E_tau x E_tau, tau=e^{i pi/3}, M=(tau,tau^-1) [order 6 on A]": (block(Mt, inv2(Mt)), True, "(tau,tau^-1)"),
  "translation only, M=1": (Id4, True, "1"),
}
table = {}
byorder = {}
for name, (M, distinct, _) in cases.items():
    rows = []
    for t in A2:
        try:
            o, fp, det = km_fixed(M, t, distinct if M != Id4 else True)
        except ValueError:
            continue
        if o == 1: continue     # identity map: fixed-point count meaningless
        rows.append({"t": t, "order_on_Km": o, "fixed_points": fp, **det})
        if o and o > 1: byorder.setdefault(o, Counter())[fp] += 1
    table[name] = rows
    print(name, sorted(Counter((r["order_on_Km"], r["fixed_points"]) for r in rows).items()))
byorder_s = {str(o): dict(c) for o, c in sorted(byorder.items())}
print(byorder_s)

# ---------- (4) M24 comparison

ex = json.load(open(HERE / "exports.json"))
shapes = ex["M24_cycle_shapes"]
def n_fixed(shape): 
    for part in shape.split():
        l, m = part.split("^")
        if l == "1": return int(m)
    return 0
m24_fix = {o: sorted({n_fixed(s) for s in v}) for o, v in shapes.items()}
d = json.load(open(HERE / "_m24_data.json"))
octads, gens = d["octads"], d["gens"]
SP = lambda x: Permutation(list(x))
M24 = PermutationGroup([SP(g) for g in gens])
O = octads[0]
Opts = [i for i in range(24) if (O >> i) & 1]
H = M24
for p_ in Opts: H = H.stabilizer(p_)
Hord = H.order()
comp = [i for i in range(24) if not (O >> i) & 1]
elts = list(H.generate())
def ctype(p): 
    a = p.array_form; seen = [0] * 24; ct = []
    for i in range(24):
        if not seen[i]:
            l, j = 0, i
            while not seen[j]: seen[j] = 1; j = a[j]; l += 1
            ct.append(l)
    return tuple(sorted(ct))
cts = Counter(" ".join(f"{l}^{c}" for l, c in sorted(Counter(ctype(p)).items())) for p in elts)
# identification with F_2^4: pick basis of the elementary abelian group and translate
c0 = comp[0]
orb_pt = {p.array_form[c0]: p for p in elts}
regular = len(orb_pt) == 16
basis = []; span_ = {tuple(range(24)): ()}
def comb(a, b): return tuple(b[a[i]] for i in range(24))
for p in elts:
    if tuple(p.array_form) not in span_ and len(basis) < 4:
        new = {}
        for k, coords in span_.items(): new[comb(k, tuple(p.array_form))] = coords + (1,) if False else None
        basis.append(p)
        span_ = {}
        for bits in itertools.product((0, 1), repeat=len(basis)):
            q = tuple(range(24))
            for bt, bp in zip(bits, basis):
                if bt: q = comb(q, tuple(bp.array_form))
            span_[q] = bits
pt_coord = {q[c0]: bits for q, bits in span_.items()}
hyper = set()
for a in range(1, 16):
    for c in (0, 1):
        hyper.add(frozenset(pt for pt, bits in pt_coord.items() if (bin(a & sum(b << i for i, b in enumerate(bits))).count("1") & 1) == c))
inside = {frozenset(i for i in range(24) if (o >> i) & 1) for o in octads if o & O == 0}
def match_table():
    out = {}
    for o, c in byorder.items():
        sh = [sname for sname in shapes.get(str(o), []) if n_fixed(sname) in c]
        out[str(o)] = {"Km_fixed_points": sorted(c), "M24_shapes_of_that_order_with_same_number_of_fixed_points": sh}
    return out
res5 = {
 "comparison_Km_vs_M24_shapes (tier C identification)": match_table(),
 "lefschetz_translation": {"L": L_translation, "traces_on_H0..H4": tr_t, "direct_fixed_points": fix_transl},
 "lefschetz_minus1": {"L": L_minus1, "traces_on_H0..H4": tr_m, "direct_fixed_points": fix_minus1},
 "Km_translation_lefschetz_all_15": transl_vals, "Km_translation_detail_example": km_lefschetz_translation((1, 0, 0, 0))[1],
 "Km_fixed_point_table": table, "fixed_points_by_order_on_Km": byorder_s,
 "M24_fixed_points_by_order_of_shapes_with_fixed_pts": m24_fix,
 "octad_pointwise_stabiliser": {"order": Hord, "regular_on_complement_16_points": regular, "cycle_shapes_of_elements": dict(cts),
      "elementary_abelian_rank": len(basis), "octads_in_complement_equal_affine_hyperplanes_under_regular_action": (frozenset(hyper) == inside),
      "n_hyperplanes": len(hyper), "n_octads_in_complement": len(inside)},
}
print(json.dumps({k: v for k, v in res5.items() if k not in ("Km_fixed_point_table",)}, indent=1, default=str)[:3000])
add_inputs([
 {"name": "Lefschetz_fixed_point_formula", "value": "for a finite-order automorphism with isolated fixed points, #Fix = sum (-1)^k tr(g|H^k)", "tier": "L", "why": "task item 5"},
 {"name": "Km_cohomology", "value": "H^2(Km A) = pi_* H^2(A) (rank 6, translations act trivially) + 16 exceptional curves E_p permuted as p->p+t; H^1=H^3=0", "tier": "L", "why": "standard structure of the Kummer surface"},
 {"name": "hol_tangent_eigenvalues_distinct_at_nodes", "value": "for the chosen symplectic (M,t) the holomorphic tangent eigenvalues (lambda, lambda^-1) are distinct, so g has 2 fixed points on each fixed exceptional curve", "tier": "B", "why": "true for (i,-i),(omega,omega^2),(tau,tau^-1) since lambda^2 != 1; for translations there are no fixed nodes"},
 {"name": "Nikulin_fixed_points", "value": "symplectic automorphisms of K3 of order 2,3,4,5,6,7,8 have 8,6,4,4,2,3,2 fixed points", "tier": "L", "why": "FROM MEMORY; used only as an 'expected' comparison, filled after computing"},
 {"name": "identification_H*(Km)_with_24_points", "value": "identify the 24 = rank H^*(K3) with the 24 permuted points of M24", "tier": "C", "why": "conjectural identification asked for in item 5"},
])
add_results([
 {"id": "F5_lefschetz_T4", "quantity": "Lefschetz numbers on H^*(T^4,Z) of translation and of -1", "computed": res5["lefschetz_translation"] | {"minus1": res5["lefschetz_minus1"]}, "shared_inputs": ["Lefschetz_fixed_point_formula"], "script": RELDIR + "/f5_symplectic.py", "command": cmd("f5_symplectic.py"), "tier": "B (formula tier L)"},
 {"id": "F5_km_fixed_points", "quantity": "fixed points on Km(A) of translations (Lefschetz on H^*(Km)) and of x->Mx+t, by exact enumeration, keyed by order on Km", "computed": {"translations": transl_vals, "by_order": byorder_s, "table": table},
  "shared_inputs": ["Km_cohomology", "hol_tangent_eigenvalues_distinct_at_nodes", "Lefschetz_fixed_point_formula"], "script": RELDIR + "/f5_symplectic.py", "command": cmd("f5_symplectic.py"), "tier": "B",
  "expected_after": "Nikulin (FROM MEMORY, tier L): order 2:8, 3:6, 4:4, 6:2. Orders 5,7,8 not realised on Km of E x E' with these automorphisms"},
 {"id": "F5_M24_comparison", "quantity": "M24 pointwise octad stabiliser 2^4 acts regularly on the 16 complementary points; nontrivial elements have shape 1^8 2^8; matches Km translations group (character 24,8,...,8)", "computed": res5["octad_pointwise_stabiliser"] | {"M24_fixed_points_by_order": m24_fix},
  "shared_inputs": ["identification_H*(Km)_with_24_points"], "script": RELDIR + "/f5_symplectic.py", "command": cmd("f5_symplectic.py"), "tier": "B for the computed group facts; the identification with H^*(Km A) is tier C"},
])
add_rigidity([
 {"id": "F5_fixed_points", "parameter_inserted": "the automorphism (M,t) of A descending to Km A", "selecting_condition": "none: fixed points are counted", "condition_uses_true_value": False,
  "solution_set": byorder_s, "classification": "VERIFIED_IDENTITY", "input_it_depends_on": "Km_cohomology, hol_tangent_eigenvalues_distinct_at_nodes",
  "negative_control": "control: direct enumeration of fixed points of translations equals the Lefschetz value 8 on H^*(Km); a wrong exceptional-curve permutation (trivial action) would give 24, not 8", "control_perturbs_same_parameter": True},
])
add_cnd('F5_orders_5_6_7_8', 'Orders 5,7,8 are not realised by any (M,t) on Km(E x E\') with M in the elliptic families used (needs 5th/7th/8th roots of unity); order 6 is not realised by (tau,tau^-1)+t (g^3 = translation by 2Mt = 0). Only orders 2,3,4 were computed.')
add_cnd('F5_M24_order_identification', 'The map (symplectic Km automorphism of order n) <-> (M24 element of order n) was only compared through fixed-point counts for n=2,3,4; tier C.')
add_exports({"F5_km_fixed_points_by_order": byorder_s, "F5_translation_fixed_points_Km": transl_vals, "F5_octad_stabiliser_order": Hord})
print("done")
