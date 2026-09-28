"""Track F part 3: rank-4 root lattices, |Aut| by enumeration, and complex structures of R^4/D4-type from order-3/4 automorphisms.
Run: cd audit/k3t2_rigidity_v3/F-whichk3 && <venv python> f3_t4.py
(no arguments). Exact integer arithmetic.
|Aut(L)| = number of 4-tuples of minimal vectors (roots) with the Gram matrix of the standard simple basis; each such tuple defines an automorphism (det Gram equal => index 1)."""
import itertools, sys
from fractions import Fraction as F
import sympy as sp
from common import *
from lat import *

def cartan_chain(n):
    return [[2 if i == j else (-1 if abs(i - j) == 1 else 0) for j in range(n)] for i in range(n)]
def blockdiag(*Ms):
    n = sum(len(M) for M in Ms); out = [[0] * n for _ in range(n)]; o = 0
    for M in Ms:
        for i in range(len(M)):
            for j in range(len(M)): out[o + i][o + j] = M[i][j]
        o += len(M)
    return out
D4 = [[2, -1, 0, 0], [-1, 2, -1, -1], [0, -1, 2, 0], [0, -1, 0, 2]]
LAT = {"A4": cartan_chain(4), "D4": D4, "A3+A1": blockdiag(cartan_chain(3), cartan_chain(1)),
       "A2+A2": blockdiag(cartan_chain(2), cartan_chain(2)), "A2+A1+A1": blockdiag(cartan_chain(2), cartan_chain(1), cartan_chain(1)),
       "A1^4": blockdiag(*[cartan_chain(1)] * 4)}

def vecs_of_norm(G, nrm):
    Gs = sp.Matrix(G); Gi = Gs.inv()
    bnd = [int(sp.floor(sp.sqrt(nrm * Gi[i, i]))) for i in range(4)]
    out = []
    for v in itertools.product(*[range(-b, b + 1) for b in bnd]):
        if sum(v[i] * G[i][j] * v[j] for i in range(4) for j in range(4)) == nrm: out.append(v)
    return out

def ip(G, u, v): return sum(u[i] * G[i][j] * v[j] for i in range(4) for j in range(4))

def automorphisms(G):
    mins = vecs_of_norm(G, 2)     # minimal vectors = roots (norm 2) for all six lattices
    res = []
    # backtracking on images of basis vectors
    def rec(k, imgs):
        if k == 4:
            res.append(tuple(imgs)); return
        for v in mins:
            if ip(G, v, v) != G[k][k]: continue
            if all(ip(G, v, imgs[j]) == G[k][j] for j in range(k)):
                rec(k + 1, imgs + [v])
    rec(0, [])
    return mins, res

auts = {}
summary = {}
for name, G in LAT.items():
    mins, A = automorphisms(G)
    # check: images generate L (det of coordinate matrix = +-1)
    ok = all(abs(sp.Matrix(a).det()) == 1 for a in A)
    summary[name] = {"det": int(sp.Matrix(G).det()), "n_roots": len(mins), "Aut_order": len(A), "all_unimodular": ok}
    auts[name] = A
    print(name, summary[name])
best = max(summary, key=lambda k: summary[k]["Aut_order"])
ties = [k for k in summary if summary[k]["Aut_order"] == summary[best]["Aut_order"]]
print("maximal:", ties)

# matrices M (columns = images of basis vectors, in basis coordinates)
def mat(a): return [[a[j][i] for j in range(4)] for i in range(4)]   # M[i][j] = j-th image's i-th coordinate
def matmul(A, B): return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]
Id = [[1 if i == j else 0 for j in range(4)] for i in range(4)]
def madd(A, B): return [[A[i][j] + B[i][j] for j in range(4)] for i in range(4)]
def order_of(M):
    P = M
    for k in range(1, 25):
        if P == Id: return k
        P = matmul(P, M)
    return None

def T_of_structure(M, lam):
    """J acts on L=Z^4 by M; holomorphic 1-forms phi satisfy phi(Mx) = lam phi(x), i.e. M^T phi = lam phi."""
    MT = [[M[j][i] for j in range(4)] for i in range(4)]
    vs = nullspace_K(MT, lam)
    if len(vs) != 2: return None
    om = wedge_coeffs(vs[0], vs[1])
    r = NS_T(om)
    return r

def ns_positive(r):
    """projective iff NS contains a class of positive square (signature (1, rho-1))."""
    return signature(r["NS_gram"])

out_D4 = {}
maxname = best if len(ties) == 1 else ties[0]
A = auts[maxname]
Gm = LAT[maxname]
strs = {}
for a in A:
    M = mat(a)
    o = order_of(M)
    M2 = matmul(M, M)
    kind = None
    if M2 == [[-x for x in row] for row in Id]: kind = "M^2=-1"
    elif madd(madd(M2, M), Id) == [[0] * 4] * 4: kind = "M^2+M+1=0"
    if kind is None: continue
    lam = I_ if kind == "M^2=-1" else OMEGA
    r = T_of_structure(M, lam)
    if r is None: continue
    T = r["T_gram"]
    key = (kind, len(r["NS_basis"]), tuple(sorted(T[i][i] for i in range(len(T)))) if T else None)
    if len(T) == 2 and T[0][0] < 0: T = [[-x for x in row] for row in T]
    red = lagrange_reduce(T[0][0], T[0][1], T[1][1]) if len(T) == 2 and T[0][0]*T[1][1]-T[0][1]**2 > 0 else ("nondef", signature(T), int(sp.Matrix(T).det()))
    k2 = (kind, len(r["NS_basis"]), red, int(sp.Matrix(T).det()), signature(r["NS_gram"]))
    strs.setdefault(k2, 0); strs[k2] += 1
print("structures (kind, rank NS, reduced T, det T, sig NS): count")
for k, v in strs.items(): print(k, v)

# J = (2 M + 1)/sqrt3 check: same holomorphic forms as lam=OMEGA (eigenvalue i sqrt3 / sqrt3 = i) -> covered by kind M^2+M+1=0
# Riemann form integrality for the order-4 case: E(x,y) = G(Jx,y) integral and antisymmetric?
riem = None
for a in A:
    M = mat(a)
    if matmul(M, M) == [[-x for x in row] for row in Id]:
        E = [[sum(M[k][i] * Gm[k][j] for k in range(4)) for j in range(4)] for i in range(4)]   # E_ij = G(M e_i, e_j)
        riem = {"E_antisymmetric": all(E[i][j] == -E[j][i] for i in range(4) for j in range(4)), "E_integral": True, "E": E}
        break

# Aut element counts by order
ordc = {}
for a in A:
    o = order_of(mat(a)); ordc[o] = ordc.get(o, 0) + 1

add_inputs([
 {"name": "rank4_root_lattices", "value": "Cartan matrices of A4,D4,A3+A1,A2+A2,A2+A1+A1,A1^4", "tier": "L", "why": "standard root lattices given in task"},
 {"name": "complex_structure_J_equals_g", "value": "complex structure on R^4/L chosen so that the automorphism g acts as i (g^2=-1) or omega (g^2+g+1=0)", "tier": "B", "why": "'compatible with an automorphism of order 3 or 4'; other J commuting with g (a continuous family for order 3) are NOT enumerated"},
])
add_results([
 {"id": "F3_aut_orders", "quantity": "|Aut(L)| by exact enumeration for rank-4 root lattices", "computed": {"per_lattice": summary, "maximal": ties},
  "shared_inputs": ["rank4_root_lattices"], "script": RELDIR + "/f3_t4.py", "command": cmd("f3_t4.py"), "tier": "B"},
 {"id": "F3_complex_structures", "quantity": f"complex structures J=g on R^4/{maxname} for every g in Aut with g^2=-1 or g^2+g+1=0: (kind, rank NS, reduced T(A), det T, signature NS) -> number of such g",
  "computed": {str(k): v for k, v in strs.items()}, "sign_note": "for g whose complex structure has orientation opposite to the coordinate orientation, NS has signature (3,1) and T is negative definite; T is reported after the global sign change (equivalent to reversing orientation)", "aut_orders_histogram": {str(k): v for k, v in ordc.items()}, "riemann_form_order4": riem,
  "shared_inputs": ["rank4_root_lattices", "complex_structure_J_equals_g"], "script": RELDIR + "/f3_t4.py", "command": cmd("f3_t4.py"), "tier": "B"},
])
add_rigidity([
 {"id": "F3_maximal_aut", "parameter_inserted": "choice of rank-4 root lattice", "selecting_condition": "maximal |Aut(L)|", "condition_uses_true_value": False,
  "solution_set": ties, "classification": "RIGID_GIVEN_DEFINITION", "input_it_depends_on": "rank4_root_lattices (the list of six candidates)",
  "negative_control": "the other five root lattices have smaller |Aut| (see F3_aut_orders); all are existing lattices", "control_perturbs_same_parameter": True},
])
add_cnd("F3_other_J", "Only J=g (g acting as i or omega) was enumerated; other complex structures commuting with g (a continuous family for order 3, e.g. eigenvalue omega on one line and omega-bar on another) were not enumerated.")
add_exports({"F3_max_aut_lattice": ties, "F3_aut_orders": {k: v["Aut_order"] for k, v in summary.items()}})
print("done")
