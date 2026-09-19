"""Skeptic (math lens): independent minimal adjudication of rows orientifold#6,#7,#8,#9,#11 and D-tda#56.
Route independent of comparison/third_computation_orientifold.py (which used TT's O7/D7 counts and
Q(O7^-) = -4 as tier-L inputs). Here the 7-brane content is derived from the F-theory / Sen-limit
geometry, with these tier-L inputs only:
  L1  Kodaira fibre Euler numbers: e(I_0*) = 6, e(I_1) = 1, and e(smooth elliptic fibre) = 0;
      for an elliptic K3 over P^1, sum of e(singular fibres) = e(K3).
  L2  Sen limit: one I_0* fibre = one O7^- plus 4 D7 (monodromy -1); D7 count read off as
      e(I_0*) - e(O7^-) with e(O7^-) = 2 (an O7^- resolves into two mutually non-local 7-branes).
  L3  each 7-brane wrapping K3 carries chi(K3)/24 units of D3 charge (per (p,q) 7-brane).
Computed: fixed points of z -> -z on T^2 and T^4 by enumeration; real dimension of the fixed locus of
each involution acting on R^{3,1} x K3 x T^2 (resp. T^4 x T^2), hence the O-plane type Op with p+1 = dim;
e(K3) from Track D's blind Betti numbers (read, not typed); chi(K3xK3) by Kunneth.
Negative control: 16 I_0* fibres (the sealed v3.20.0 TadpoleCancellation O7 count) would need
sum e = 96 > e(K3), impossible.
D-tda#56: sigma(K3) from Track C's derived signature (read), Dirac index -sigma/8; Todd genus chi(O)
(Noether, c1 = 0) and A-hat = (2c2 - c1^2)/24 equal iff c1^2 = 0.
Command: cd audit/k3t2_rigidity_v2/skeptic_math && \
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python skX_rows_orientifold_dirac.py
"""
import itertools, json, os
from fractions import Fraction as Fr

here = os.path.dirname(os.path.abspath(__file__))
V2 = os.path.dirname(here)


def fixed_pts(n):  # x in (1/2^k)Z^n/Z^n, k=3 grid, with -x == x mod 1
    g = [Fr(i, 8) for i in range(8)]
    return [x for x in itertools.product(g, repeat=n) if all((2 * xi) % 1 == 0 for xi in x)]


nT2, nT4 = len(fixed_pts(2)), len(fixed_pts(4))
# fixed-locus dimension: spacetime R^{3,1} (4) + K3 (4) + T^2 (2); involution acts as -1 on T^2 only
dim_O_T2 = 4 + 4 + 0          # K3 fixed pointwise, T^2 -> points
dim_O_T4 = 4 + 0 + 2          # T^4 -> points (Kummer involution), T^2 untouched
D = json.load(open(os.path.join(V2, "D-tda", "results.json")))
b = [D["resolved_K3_betti_numbers"]["N=6"][f"b{k}"] for k in range(5)]
eK3 = sum((-1) ** k * x for k, x in enumerate(b))
# Kunneth
bb = [0] * 9
for i, x in enumerate(b):
    for j, y in enumerate(b): bb[i + j] += x * y
eK3K3 = sum((-1) ** k * x for k, x in enumerate(bb))
e_I0star, e_O7 = 6, 2                       # L1, L2
n_I0star = nT2                              # one I_0* over each branch point of T^2 -> P^1 (Sen limit)
sum_e = n_I0star * e_I0star
n_D7 = n_I0star * (e_I0star - e_O7)         # D7 per I_0* = 4
# D7-charge: monodromy of I_0* is -1 (total 7-brane charge of the stack zero in the Sen limit)
Q_O7_in_D7_units = -(e_I0star - e_O7)
n_7branes = sum_e                           # number of (p,q) 7-branes = sum of Euler numbers (I_1 count)
D3_per_7brane = Fr(eK3, 24)                 # L3
D3_total = n_7branes * D3_per_7brane
D3_split = {"O7": n_I0star * e_O7 * D3_per_7brane, "D7": n_D7 * D3_per_7brane}
neg = {"n_I0star_if_16_O7": nT4, "sum_e": nT4 * e_I0star, "exceeds_e(K3)": nT4 * e_I0star > eK3}
C = json.load(open(os.path.join(V2, "C-lattices", "exports.json")))
p, q = C["derived_signature_p_q"]
sigma = p - q
dirac = Fr(-sigma, 8)
c1sq, c2 = 0, eK3
todd = Fr(c1sq + c2, 12); ahat = Fr(2 * c2 - c1sq, 24); hirz_sigma = Fr(c1sq - 2 * c2, 3)
out = {
    "fixed_points": {"T2": nT2, "T4": nT4},
    "O_plane_from_fixed_locus_dim": {"involution_on_T2": f"O{dim_O_T2 - 1}", "involution_on_T4": f"O{dim_O_T4 - 1}"},
    "e_K3_from_trackD_betti": eK3, "betti_K3_trackD": b,
    "F_theory": {"n_I0star": n_I0star, "sum_e_singular_fibres": sum_e, "equals_e_K3": sum_e == eK3,
                 "n_D7": n_D7, "n_O7": n_I0star, "Q_O7_in_D7_units": Q_O7_in_D7_units,
                 "total_O7_charge": n_I0star * Q_O7_in_D7_units, "total_D7_charge": n_D7,
                 "n_7branes": n_7branes, "D3_total": str(D3_total),
                 "D3_split": {k: str(v) for k, v in D3_split.items()},
                 "chi_K3xK3": eK3K3, "chi_K3xK3_over_24": str(Fr(eK3K3, 24)),
                 "chi_K3_over_24_(per_7brane)": str(D3_per_7brane)},
    "negative_control_16_O7": neg,
    "D-tda#56": {"sigma_from_trackC_signature": sigma, "dirac_index_minus_sigma_over_8": str(dirac),
                 "hirzebruch_sigma_(c1^2-2c2)/3": str(hirz_sigma), "todd_chiO": str(todd), "ahat": str(ahat),
                 "todd_equals_ahat_here": todd == ahat,
                 "todd_minus_ahat_symbolic": str(__import__("sympy").simplify(
                     (lambda a, c: (a + c) / 12 - (2 * c - a) / 24)(*__import__("sympy").symbols("c1sq c2"))))},
}
out["expected_after_computing"] = {
    "source": "comparison.json lean_value (sealed v3.20.0)",
    "orientifold#6_total_d7_charge": 64, "orientifold#7_total_o7_charge": -64,
    "orientifold#8_totalO7Charge": -64, "orientifold#9_totalD7Charge": 64,
    "orientifold#11_d3TadpoleTarget": 1, "D-tda#56_DiracIndex(-16)": 2,
}
json.dump(out, open(os.path.join(here, "skX_rows_orientifold_dirac_results.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
