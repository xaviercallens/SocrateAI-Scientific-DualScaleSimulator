"""Skeptic check for Track E: is the (A.6)-(A.7) diagonal parametrisation of Gamma_{3,3} faithful
to TT's integral lattice (A.2)-(A.3)?  Exact integer arithmetic only.
 - det of U^3 (TT A.3) vs det of 2*eta_{3,3} (TT A.6/A.7 Gram): index of the diagonal lattice.
 - rank-4 (U+U): enumerate flux pairs (alpha,beta) in 2*(U+U) (TT eq 2.5: even coefficients in a
   Z-basis of H^2(K3,Z)) with alpha^2=beta^2>0, alpha.beta=0, alpha^2<=24, hyperbolic coefficients
   |c|<=BOUND; count how many are ALSO in Track E's family 2*L_diag (even diagonal coordinates).
 - 8 | alpha^2 for all alpha in 2*Gamma (Gamma even): checked on every enumerated vector.
Command: cd audit/k3t2_rigidity_v2/skeptic_math && <venv python> skE_lattice.py
"""
import itertools, json, os
import sympy
U = sympy.Matrix([[0, 1], [1, 0]])
U3 = sympy.diag(U, U, U)
D = sympy.diag(2, 2, 2, -2, -2, -2)
out = {"det_U3": int(U3.det()), "det_diag_2eta": int(D.det()),
       "index_of_diag_sublattice": int(sympy.sqrt(abs(D.det() / U3.det())))}
# rank 4: hyperbolic coords (a1,b1,a2,b2) -> vector a1 e1 + b1 f1 + a2 e2 + b2 f2
def gram(x, y):
    return x[0] * y[1] + x[1] * y[0] + x[2] * y[3] + x[3] * y[2]
# diagonal basis: u_i = e_i + f_i (norm 2), w_i = e_i - f_i (norm -2);
# a e + b f = ((a+b)/2) u + ((a-b)/2) w  -> diagonal coordinates
def diag_coords(x):
    return [sympy.Rational(x[0] + x[1], 2), sympy.Rational(x[2] + x[3], 2),
            sympy.Rational(x[0] - x[1], 2), sympy.Rational(x[2] - x[3], 2)]
def in_E_family(x):  # E: coefficients in diagonal basis are even integers
    return all(c.q == 1 and c.p % 2 == 0 for c in diag_coords(x))
BOUND = 3
vecs = [tuple(2 * c for c in v) for v in itertools.product(range(-BOUND, BOUND + 1), repeat=4)]
pos = {}
for v in vecs:
    n = gram(v, v)
    assert n % 8 == 0
    if 0 < n <= 24:
        pos.setdefault(n, []).append(v)
tot = inE = 0
example_missing = None
by_norm = {}
for n, vs in pos.items():
    t = e = 0
    for a in vs:
        for b in vs:
            if gram(a, b) == 0:
                t += 1
                if in_E_family(a) and in_E_family(b):
                    e += 1
                elif example_missing is None:
                    example_missing = {"alpha_hyperbolic_coords": a, "beta_hyperbolic_coords": b,
                                       "alpha_diag_coords": [str(c) for c in diag_coords(a)],
                                       "beta_diag_coords": [str(c) for c in diag_coords(b)],
                                       "alpha_sq": n}
    by_norm[n] = {"TT_admissible_pairs": t, "of_which_in_trackE_family": e}
    tot += t; inE += e
out.update({"rank4_hyperbolic_coefficient_bound_(even_entries_2c,|c|<=)": BOUND,
            "all_alpha_sq_divisible_by_8": True,
            "alpha_sq_values_realised": sorted(pos),
            "pairs_by_alpha_sq": by_norm, "total_pairs": tot, "in_trackE_family": inE,
            "example_TT_admissible_pair_missed_by_trackE": example_missing,
            "TT_4.15_example_2e1diag_in_hyperbolic_coords": [2, 2, 0, 0],
            "reading": "E's diagonal-basis family is a proper subfamily (index-8 sublattice) of TT's integral "
                       "flux lattice; 8|alpha^2 and the admissible set {8,16,24} survive (Gamma even), but E's "
                       "counts undercount TT-admissible flux pairs at any stated bound."})
# Tadpole bookkeeping re-derived: O7 at fixed points of z->-z on T^2 = (1/2)Z^2/Z^2
from fractions import Fraction as Fr
fp_T2 = [(a, b) for a in (Fr(0), Fr(1, 4), Fr(1, 2), Fr(3, 4)) for b in (Fr(0), Fr(1, 4), Fr(1, 2), Fr(3, 4))
         if (2 * a).denominator == 1 and (2 * b).denominator == 1]
nO7 = len(fp_T2); qO7_D7units = -4
nD7 = -nO7 * qO7_D7units          # 7-brane tadpole: nD7 * 1 + nO7 * (-4) = 0
D3_total = nO7 * 2 + nD7 * 1      # TT lines 164-169: 2 per O7, 1 per D7 (tier L)
betti = [1, 0, 22, 0, 1]
bK3K3 = [sum(betti[i] * betti[k - i] for i in range(5) if 0 <= k - i < 5) for k in range(9)]
chiK3K3 = sum((-1) ** k * b for k, b in enumerate(bK3K3))
# TT (4.15) in diag basis: alpha=2e1, beta=2e2; alpha_y=-beta, beta_y=alpha (4.8); N_flux (2.8)
dg = [2, 2, -2, -2]
dotd = lambda x, y: sum(g * a * b for g, a, b in zip(dg, x, y))
ax, bx = [2, 0, 0, 0], [0, 2, 0, 0]; ay = [-c for c in bx]; by = ax
Nflux = -dotd(bx, ay) + dotd(by, ax)
out["tadpole"] = {"n_O7_T2_fixed_points": nO7, "n_D7": nD7, "D3_total": D3_total,
                  "chi_K3xK3_over_24": chiK3K3 // 24, "chi_K3xK3": chiK3K3,
                  "TT_4.15_N_flux": Nflux, "TT_4.16_N_D3": D3_total - Nflux // 2}
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "skE_lattice_results.json"), "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
