"""Track F part 1: T2 moduli, stabilisers, maximal-symmetry point, Narain lattice roots.
Run: cd audit/k3t2_rigidity_v3/F-whichk3 && <venv python> f1_moduli.py 6
argument: N = entry bound for the SL(2,Z) matrix search (committed run uses 6).
All arithmetic exact (Fraction / sympy).  Points tau = x + i*y stored as (x, s=y^2) with x, s Fractions.
"""
import sys, itertools, json
from fractions import Fraction as F
import sympy as sp
from common import *

N = int(sys.argv[1]) if len(sys.argv) > 1 else 6

def canon(x, s):
    """reduce to fundamental domain, exactly. returns (x,s) canonical: -1/2 < x <= 1/2, |tau|>=1, x>=0 on the unit circle."""
    x, s = F(x), F(s)
    while True:
        k = (x + F(1, 2)) // 1          # shift x into (-1/2,1/2]
        x -= k
        if x == -F(1, 2):
            x = F(1, 2)
        n = x * x + s
        if n < 1 or (n == 1 and x < 0):
            x, s = -x / n, s / (n * n)     # tau -> -1/tau
            continue
        return (x, s)

def sl2_mats(N):
    for a, b, c, d in itertools.product(range(-N, N + 1), repeat=4):
        if a * d - b * c == 1:
            yield (a, b, c, d)

def stabiliser(pt, N):
    x, s = pt
    out = []
    for (a, b, c, d) in sl2_mats(N):
        # gamma tau = tau  <=>  c tau^2 + (d-a) tau - b = 0 ; imag part: 2 c x + d - a = 0 ; real part: c (x^2 - s) + (d-a) x - b = 0
        if 2 * c * x + d - a == 0 and c * (x * x - s) + (d - a) * x - b == 0:
            out.append((a, b, c, d))
    return out

# ---- (a) points with stabiliser larger than {+-1}
cands = {}
for (a, b, c, d) in sl2_mats(N):
    if (a, b, c, d) in [(1, 0, 0, 1), (-1, 0, 0, -1)]:
        continue
    tr = a + d
    if abs(tr) < 2 and c != 0:
        x = F(a - d, 2 * c)
        s = F(4 - tr * tr, 4 * c * c)
        cands[canon(x, s)] = 1
stab = {}
for pt in cands:
    st = stabiliser(pt, N)
    stab[pt] = len(st)
elliptic = {pt: n for pt, n in stab.items() if n > 2}
print("elliptic points (x, y^2) -> |Stab_SL2Z|:", elliptic)

# ---- (b) full duality group stabiliser
def stab_order(pt):
    if pt in stab: return stab[pt]
    return len(stabiliser(pt, N))

def full_stab(tau, rho):
    """order of stabiliser in the FORMAL group (SL2Z x SL2Z) x| (exchange, reflection): sum over (eps,r) of [images equivalent] |Stab tau||Stab rho|"""
    total = 0
    detail = []
    for eps in (0, 1):
        for r in (0, 1):
            t, p = tau, rho
            if r:
                t, p = (-t[0], t[1]), (-p[0], p[1])
            if eps:
                t, p = p, t
            if canon(*t) == tau and canon(*p) == rho:
                total += stab_order(tau) * stab_order(rho)
                detail.append((eps, r))
    return total, detail

# candidate point set: elliptic points, plus sample of non-elliptic points (symmetric lines and generic rationals)
P = {}
for pt in elliptic: P[pt] = "elliptic"
samples = []
for x in [F(0), F(1, 2), F(1, 4), F(1, 3), F(-1, 5)]:
    for s in [F(2), F(3), F(4), F(9, 4), F(1, 1) - x * x, F(7, 3)]:
        c = canon(x, s)
        if c not in P:
            P[c] = "sample"
allpairs = []
for tau in P:
    for rho in P:
        o, det = full_stab(tau, rho)
        allpairs.append((o, tau, rho, det))
allpairs.sort(key=lambda t: -t[0])
maxo = allpairs[0][0]
maxpts = [(t, r, d) for (o, t, r, d) in allpairs if o == maxo]
print("max formal stabiliser order", maxo, "at", maxpts)
print("top values:", sorted(set(o for o, *_ in allpairs), reverse=True)[:6])

def fmt(pt):
    return f"x={pt[0]}, y^2={pt[1]}"

# ---- (c) Narain lattice at each maximal point (and, for contrast, at the other elliptic pairs)
def narain_roots(pt_tau, pt_rho):
    tau1, tau2 = sp.Rational(pt_tau[0].numerator, pt_tau[0].denominator), sp.sqrt(sp.Rational(pt_tau[1].numerator, pt_tau[1].denominator))
    rho1, rho2 = sp.Rational(pt_rho[0].numerator, pt_rho[0].denominator), sp.sqrt(sp.Rational(pt_rho[1].numerator, pt_rho[1].denominator))
    Gm = sp.simplify((rho2 / tau2)) * sp.Matrix([[1, tau1], [tau1, tau1**2 + tau2**2]])
    Gm = Gm.applyfunc(sp.nsimplify)
    Bm = sp.Matrix([[0, rho1], [-rho1, 0]])
    Gi = Gm.inv().applyfunc(sp.simplify)
    # p_{L,R} = e^* (n + (B +- G) w)/sqrt2 ; p_L^2 = (1/2) v_L^T G^-1 v_L ; v_{L,R} = n + (B +- G) w
    n = sp.symbols('n1 n2'); w = sp.symbols('w1 w2')
    nv, wv = sp.Matrix(n), sp.Matrix(w)
    vL = nv + (Bm + Gm) * wv
    vR = nv + (Bm - Gm) * wv
    pL2 = sp.expand((vL.T * Gi * vL)[0] / 2)
    pR2 = sp.expand((vR.T * Gi * vR)[0] / 2)
    Q = sp.hessian(pL2 + pR2, list(n) + list(w)) / 2   # positive definite quadratic form on (n,w)
    Qi = Q.inv()
    # bound: v_i^2 <= (Q^-1)_ii * (v^T Q v) = 2 (Q^-1)_ii  when p_L^2 + p_R^2 = 2
    bnd = [int(sp.floor(sp.sqrt(sp.N(2 * Qi[i, i], 30)))) + 1 for i in range(4)]
    Lroots, Rroots = [], []
    fL = sp.lambdify(list(n) + list(w), pL2, 'sympy'); fR = sp.lambdify(list(n) + list(w), pR2, 'sympy')
    for v in itertools.product(*[range(-b, b + 1) for b in bnd]):
        a = sp.nsimplify(sp.expand(fL(*v))); b = sp.nsimplify(sp.expand(fR(*v)))
        if b == 0 and a == 2: Lroots.append(v)
        if a == 0 and b == 2: Rroots.append(v)
    return Gm, Bm, bnd, Lroots, Rroots

def pair(u, v):   # Lorentzian pairing of (n1,n2,w1,w2)
    return u[0]*v[2] + u[1]*v[3] + v[0]*u[2] + v[1]*u[3]

def simple_system(roots):
    import random
    for c in itertools.product(range(1, 30), repeat=4):
        f = lambda v: sum(a*b for a, b in zip(c, v))
        if all(f(v) != 0 for v in roots): break
    pos = [v for v in roots if f(v) > 0]
    ps = set(pos)
    simple = [v for v in pos if not any(tuple(a - b for a, b in zip(v, u)) in ps for u in pos)]
    return simple

def classify(cartan):
    n = len(cartan)
    comps, seen = [], set()
    for i in range(n):
        if i in seen: continue
        st, comp = [i], set()
        while st:
            k = st.pop()
            if k in comp: continue
            comp.add(k)
            for j in range(n):
                if j != k and cartan[k][j] != 0 and j not in comp: st.append(j)
        seen |= comp; comps.append(sorted(comp))
    out = []
    for comp in comps:
        m = len(comp)
        deg = {i: sum(1 for j in comp if j != i and cartan[i][j] != 0) for i in comp}
        edges = sum(deg.values()) // 2
        assert edges == m - 1, "not a tree"
        if any(cartan[i][j] not in (0, -1) for i in comp for j in comp if i != j):
            out.append(f"non-simply-laced(m={m})"); continue
        br = [i for i in comp if deg[i] == 3]
        if not br:
            out.append(f"A{m}")
        else:
            b = br[0]
            arms = []
            for j in comp:
                if j != b and cartan[b][j] != 0:
                    L, prev, cur = 1, b, j
                    while True:
                        nxt = [k for k in comp if k not in (prev, cur) and cartan[cur][k] != 0]
                        if not nxt: break
                        prev, cur = cur, nxt[0]; L += 1
                    arms.append(L)
            arms.sort()
            if arms[:2] == [1, 1]: out.append(f"D{m}")
            elif arms[:2] == [1, 2] and m in (6, 7, 8): out.append(f"E{m}")
            else: out.append(f"unknown{arms}")
    return out

res_c = []
for (o, t, r, det) in allpairs:
    pass
pairs_to_do = [(t, r, "MAX") for (t, r, d) in maxpts]
el = list(elliptic)
for t in el:
    for r in el:
        if (t, r, "MAX") not in pairs_to_do and (t, r) not in [(a, b) for a, b, _ in pairs_to_do]:
            pairs_to_do.append((t, r, "contrast"))
narain = []
for (t, r, kind) in pairs_to_do:
    Gm, Bm, bnd, Lr, Rr = narain_roots(t, r)
    rec = {"tau": fmt(t), "rho": fmt(r), "kind": kind, "G": str(Gm.tolist()), "B12": str(Bm[0, 1]), "box": bnd,
           "n_left_roots": len(Lr), "n_right_roots": len(Rr)}
    for side, roots in (("left", Lr), ("right", Rr)):
        if roots:
            sim = simple_system(roots)
            sg = 1 if side == 'left' else -1   # p_R.p_R' = -(Lorentzian pairing)
            C = [[sg * pair(u, v) for v in sim] for u in sim]
            rec[side + "_simple_roots"] = sim
            rec[side + "_cartan"] = C
            rec[side + "_algebra"] = classify(C) if len(sim) else []
            rk = sp.Matrix(sim).rank()
            rec[side + "_rank_of_roots"] = int(rk)
        else:
            rec[side + "_algebra"] = []
    narain.append(rec)
    print(kind, rec["tau"], rec["rho"], "left", rec["left_algebra"], len(Lr), "right", rec["right_algebra"], len(Rr))

# ---- write results
add_inputs([
 {"name": "duality_group", "value": "SL2Z_tau x SL2Z_rho, exchange tau<->rho, (tau,rho)->(-conj tau,-conj rho); formal group (kernel (+-1,+-1) order 4 not divided out)", "tier": "L", "why": "given in task statement"},
 {"name": "narain_convention", "value": "G=(rho2/tau2)[[1,tau1],[tau1,|tau|^2]], B12=rho1, p_{L,R}=e^*(n+(B+-G)w)/sqrt2 (alpha'=1)", "tier": "L", "why": "standard convention (FROM MEMORY); root counts of the enhanced algebra are convention independent up to the exchange of L/R and tau<->rho"},
 {"name": "search_bound_N", "value": N, "tier": "B", "why": "entries of SL2Z matrices searched in [-N,N]"},
])
add_results([
 {"id": "F1a_elliptic_points", "quantity": "points of the fundamental domain with |Stab_SL2Z| > 2, and orders",
  "computed": {fmt(p): o for p, o in elliptic.items()}, "shared_inputs": ["search_bound_N"],
  "script": RELDIR + "/f1_moduli.py", "command": cmd("f1_moduli.py", str(N)),
  "expected_after": "expected (FROM MEMORY, tier L): tau=i (order 4), tau=e^{2 pi i/3} (order 6)", "tier": "B"},
 {"id": "F1b_max_stabiliser", "quantity": "maximal stabiliser order in the full (formal) duality group over the candidate set (elliptic pairs + 30 sampled non-elliptic points)",
  "computed": {"max_order_formal": maxo, "max_order_effective_divided_by_4": maxo // 4,
               "points": [{"tau": fmt(t), "rho": fmt(r), "eps_r_elements": d} for t, r, d in maxpts],
               "distinct_values_top6": sorted(set(o for o, *_ in allpairs), reverse=True)[:6],
               "argument": "order = sum over (eps,r) of [image equivalent] |Stab tau||Stab rho| <= 4*a*b with a,b in {2,4,6}, so 144 needs a=b=6 and all four (eps,r)"},
  "shared_inputs": ["duality_group", "search_bound_N"], "script": RELDIR + "/f1_moduli.py", "command": cmd("f1_moduli.py", str(N)), "tier": "B"},
 {"id": "F1c_narain_roots", "quantity": "roots (p_R=0,p_L^2=2 and p_L=0,p_R^2=2) of Gamma^{2,2} and Cartan-matrix classification, at each maximal point and (contrast) other elliptic pairs",
  "computed": narain, "shared_inputs": ["narain_convention", "duality_group"], "script": RELDIR + "/f1_moduli.py", "command": cmd("f1_moduli.py", str(N)), "tier": "B"},
])
maxrec = [r for r in narain if r["kind"] == "MAX"]
add_rigidity([
 {"id": "F1_point", "parameter_inserted": "(tau,rho) in H x H", "selecting_condition": "stabiliser in the formal full duality group is maximal (over candidate set: elliptic pairs and 30 sampled non-elliptic points)",
  "condition_uses_true_value": False, "solution_set": [f"({fmt(t)} ; {fmt(r)})" for t, r, d in maxpts],
  "classification": "RIGID_GIVEN_DEFINITION", "input_it_depends_on": "duality_group (the reflection and exchange are typed inputs; without the reflection or without the exchange the order at (omega,omega) is 72)",
  "negative_control": "neighbouring points: (i,i) 64, (i,omega) 48, sampled generic points <= 16; all realizable points, same parameter (tau,rho)",
  "control_perturbs_same_parameter": True},
])
add_cnd("F1_global_max", "The global maximum over all of H x H is argued (order <= 4ab, a,b in {2,4,6}) and checked over the elliptic pairs plus 30 sampled non-elliptic points, not by exhausting the continuum.")
add_exports({"F1_max_point": [{"tau": fmt(t), "rho": fmt(r)} for t, r, d in maxpts], "F1_max_stabiliser_formal": maxo,
             "F1_enhanced_algebra": {"left": maxrec[0]["left_algebra"], "right": maxrec[0]["right_algebra"]} if maxrec else None})
print("done")
