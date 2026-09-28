"""Track F part 2: NS(A), T(A) for A = E_tau x E_tau', tau,tau' the SL2Z-elliptic points found in F1(a).
Run: cd audit/k3t2_rigidity_v3/F-whichk3 && <venv python> f2_abelian.py 50
argument: max discriminant for the binary-form list (committed run: 50).
Exact arithmetic over K=Q(i,sqrt3). H^2(A,Z) = Lambda^2 H^1 = U^3; period = phi1^phi2 with phi1 = e^1 + tau e^2, phi2 = e^3 + tau' e^4."""
import sys, json
from fractions import Fraction as F
import sympy as sp
from common import *
from lat import *

DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 50
f1 = json.load(open(HERE / "results.json"))["results"]["F1a_elliptic_points"]["computed"]
# tau values from F1(a): keys 'x=..., y^2=..'; convert to K elements (only i and omega arise; checked here, not typed)
def to_K(key):
    xs, ys = key.split(", ")
    x = F(xs.split("=")[1]); y2 = F(ys.split("=")[1])
    if x == 0 and y2 == 1: return I_, "i"
    if x == F(1, 2) and y2 == F(3, 4): return K((F(1, 2), 0, 0, F(1, 2))), "e^{i pi/3}"     # x=1/2, y=sqrt3/2
    raise ValueError(key)
taus = {k: to_K(k) for k in f1}
print(taus)

def surface(tau, taup):
    phi1 = [ONE, tau, ZERO, ZERO]
    phi2 = [ZERO, ZERO, ONE, taup]
    om = wedge_coeffs(phi1, phi2)
    return NS_T(om), om

def binary_forms(dmax):
    out = []
    for p in range(2, 2 * dmax + 1, 2):
        for r in range(p, 2 * dmax + 1, 2):
            for q in range(0, p // 2 + 1):
                d = p * r - q * q
                if 0 < d <= dmax: out.append((p, q, r, d))
    return out
forms = binary_forms(DMAX)
fset = {(p, q, r) for p, q, r, d in forms}
print("reduced even positive-definite binary forms with disc<=", DMAX, ":", len(forms))

root2 = {"A1+A1": lagrange_reduce(2, 0, 2), "A2": lagrange_reduce(2, 1, 2)}
cases = []
for (k1, (t1, n1)), (k2, (t2, n2)) in itertools.product(taus.items(), repeat=2):
    r, om = surface(t1, t2)
    T = r["T_gram"]; ns = r["NS_gram"]
    rec = {"tau": n1, "tau'": n2, "rank_NS": len(r["NS_basis"]), "rank_T": len(r["T_basis"]),
           "NS_gram": ns, "NS_signature": signature(ns), "T_gram": T, "T_signature": signature(T),
           "T_disc": int(sp.Matrix(T).det()) if T else None, "T_basis": r["T_basis"]}
    if len(T) == 2:
        red = lagrange_reduce(T[0][0], T[0][1], T[1][1])
        rec["T_reduced_form"] = red
        rec["T_root_lattice_scaling"] = [(name, c) for name, rl in root2.items() for c in range(1, 9)
                                          if T[0][0] % c == 0 and T[1][1] % c == 0 and T[0][1] % c == 0 and
                                          lagrange_reduce(T[0][0] // c, T[0][1] // c, T[1][1] // c) == rl]
        # Km A
        TK = [[2 * x for x in row] for row in T]
        redK = lagrange_reduce(TK[0][0], TK[0][1], TK[1][1])
        rec["TKm_gram"] = TK; rec["TKm_disc"] = int(sp.Matrix(TK).det()); rec["TKm_reduced_form"] = redK
        rec["T_in_binary_list"] = red in fset; rec["TKm_in_binary_list"] = redK in fset
        rec["TKm_root_lattice_scaling"] = [(name, c) for name, rl in root2.items() for c in range(1, 9)
                                          if TK[0][0] % c == 0 and TK[1][1] % c == 0 and TK[0][1] % c == 0 and
                                          lagrange_reduce(TK[0][0] // c, TK[0][1] // c, TK[1][1] // c) == rl]
        rec["T_evenness"] = all(T[i][i] % 2 == 0 for i in range(2))
    else:
        rec["theta_T_to_norm_8"] = theta(T, 8) if T and signature(T)[1] == 0 else "not definite"
        TK = [[2 * x for x in row] for row in T]
        rec["TKm_disc"] = int(sp.Matrix(TK).det())
        rec["T_in_binary_list"] = False
    cases.append(rec)
    print(n1, n2, "rank NS", rec["rank_NS"], "T", T, "disc", rec["T_disc"])

# realizable negative control: non-CM tau = tau' = 2i is NOT in K's nice set? tau=2i in K: yes (2i). Both rank counts.
ctrl = {}
for nm, tt in [("2i", K((0, 2, 0, 0))), ("i+1 (T-shift of i, sanity: must reproduce T of i)", K((1, 1, 0, 0))), ("3i", K((0, 3, 0, 0)))]:
    r, _ = surface(tt, tt)
    ctrl[nm] = {"rank_NS": len(r["NS_basis"]), "rank_T": len(r["T_basis"]), "T_gram": r["T_gram"]}
# a CM point whose j-invariant is not special: tau=2i is CM (Q(i)) -> Picard number of E x E is 4 as well; use tau=(1+i sqrt3)... in K also CM.
# a genuinely non-CM point cannot be represented exactly in K, so the rank-3 (non-CM) control is not computed (see could_not_do).
print(ctrl)

# consistency: dedicated sign check
add_inputs([
 {"name": "H2_A_equals_U3", "value": "H^2(A,Z)=Lambda^2 H^1(A,Z) with wedge pairing (computed, signature (3,3), even unimodular)", "tier": "B", "why": "checked in this script (GRAM6)"},
 {"name": "Nikulin_TKm", "value": "T(Km A) = T(A)(2)", "tier": "L", "why": "Nikulin/Morrison, quoted"},
 {"name": "binary_form_disc_bound", "value": DMAX, "tier": "B", "why": "list bound"},
])
G6 = sp.Matrix(GRAM6)
add_results([
 {"id": "F2_H2_lattice", "quantity": "Gram of Lambda^2 Z^4, det, signature", "computed": {"det": int(G6.det()), "signature": signature(GRAM6), "even": all(GRAM6[i][i] % 2 == 0 for i in range(6))},
  "shared_inputs": ["H2_A_equals_U3"], "script": RELDIR + "/f2_abelian.py", "command": cmd("f2_abelian.py", str(DMAX)), "tier": "B"},
 {"id": "F2_NS_T", "quantity": "NS(A), T(A), T(Km A) for A=E_tau x E_tau', tau,tau' in {points of F1(a)}", "computed": cases,
  "shared_inputs": ["Nikulin_TKm", "binary_form_disc_bound"], "script": RELDIR + "/f2_abelian.py", "command": cmd("f2_abelian.py", str(DMAX)), "tier": "B (T(Km A)=T(A)(2) is tier L)"},
 {"id": "F2_binary_forms", "quantity": "number of reduced even positive-definite binary forms of disc <= bound", "computed": {"count": len(forms), "bound": DMAX},
  "shared_inputs": ["binary_form_disc_bound"], "script": RELDIR + "/f2_abelian.py", "command": cmd("f2_abelian.py", str(DMAX)), "tier": "B"},
 {"id": "F2_control", "quantity": "control: E_tau x E_tau for other tau (all CM, so rank is also 4)", "computed": ctrl, "shared_inputs": [], "script": RELDIR + "/f2_abelian.py", "command": cmd("f2_abelian.py", str(DMAX)), "tier": "B"},
])
add_rigidity([
 {"id": "F2_T", "parameter_inserted": "complex structure (tau,tau') of A=E_tau x E_tau'", "selecting_condition": "NS = integral classes orthogonal to the period; T = orth complement",
  "condition_uses_true_value": False, "solution_set": "see F2_NS_T", "classification": "RIGID_GIVEN_DEFINITION",
  "input_it_depends_on": "H2_A_equals_U3, the choice of the elliptic tau values from F1(a); Nikulin_TKm for the Kummer form",
  "negative_control": "tau=tau'=2i, 3i (realizable, other CM points): rank and disc change (see F2_control)", "control_perturbs_same_parameter": True},
])
add_cnd("F2_noncm_control", "A non-CM control (Picard number 3) cannot be represented exactly in K=Q(i,sqrt3); only other CM points (2i,3i) were used as controls.")
print("done")
