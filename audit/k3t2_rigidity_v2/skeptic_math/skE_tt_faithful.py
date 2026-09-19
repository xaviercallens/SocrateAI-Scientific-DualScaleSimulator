"""Skeptic (math lens): line-by-line faithfulness of Track E to Tripathy-Trivedi (pinned text
sources/hep-th_0301139_TripathyTrivedi.txt), with exact arithmetic.

 (1) TT App. A (A.6) forms WITHOUT the sqrt2 prefactor (TT line 1881-1882 claims they 'form a basis of
     H^2(T^4,Z)'), expanded in the six elementary forms dx^a^dx^b; wedge pairing computed from
     permutation signs with TT's orientation (A.7) (integral over T^4 = 2 x (-1/2) = -1).
     Output: Gram of elementary basis (expect +-U^3, det -1), Gram of the (A.6) combos, index of
     their span in H^2(T^4,Z) = |det(coefficient matrix)|.  TT lines 1821-1822 'change of basis such
     that H33 = 2 eta33' is tested as an integral statement via det.
 (2) Which Gram TT's own worked examples use: (4.15)-(4.16), (4.31)-(4.32) and the §5 example
     (lines 1185-1198) evaluated with the diagonal Gram 2*eta and with the hyperbolic Gram (A.3).
 (3) 8 | alpha^2 in the literal reading of TT eq (2.5)/line 223 (even coefficients in a Z-basis of the
     unimodular Gamma_{3,19} of (A.2)-(A.4)): Gamma even  =>  alpha = 2v, alpha^2 = 4 v^2 in 8Z.
     Checked symbolically on a general v; Gamma_{3,19} Gram typed from (A.2)-(A.4) (det, evenness).
     Negative control: an ODD unimodular lattice of the same signature allows alpha^2 = 4.
 (4) Tadpole arithmetic: (2.3) with (4.13); TT line 941 prose vs (2.3); TT line 1197 and (4.32).
 (5) Track E's orbit group (coordinate permutations within sign blocks x sign flips, order 64 for
     rank 4): fraction of elements that reverse the orientation of the positive-definite block.
Command: cd audit/k3t2_rigidity_v2/skeptic_math && \
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python skE_tt_faithful.py
"""
import itertools, json, os
import sympy
from fractions import Fraction as Fr

out = {}
# ---------- (1) wedge pairing on T^4, coords (x1,x2,y1,y2) = 0,1,2,3 ----------
pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
def perm_sign(p):
    s = 1
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]: s = -s
    return s
VOL_INTEGRAL = -1  # TT (A.7): -1/2 on T^4/Z2, i.e. -1 on the cover T^4
def wedge(a, b):
    if len(set(a + b)) < 4: return 0
    return perm_sign(a + b) * VOL_INTEGRAL
G_elem = sympy.Matrix(6, 6, lambda i, j: wedge(pairs[i], pairs[j]))
w = {p: i for i, p in enumerate(pairs)}
def form(*terms):
    v = [0] * 6
    for c, p in terms: v[w[p]] += c
    return v
# TT (A.6) without sqrt2 (x1=0,x2=1,y1=2,y2=3)
A6 = [form((1, (0, 1)), (-1, (2, 3))),   # e1 = dx1dx2 - dy1dy2
      form((1, (0, 3)), (-1, (1, 2))),   # e2 = dx1dy2 - dx2dy1
      form((1, (0, 2)), (1, (1, 3))),    # e3 = dx1dy1 + dx2dy2
      form((1, (0, 1)), (1, (2, 3))),    # e4 = dx1dx2 + dy1dy2
      form((1, (0, 3)), (1, (1, 2))),    # e5 = dx1dy2 + dx2dy1
      form((1, (0, 2)), (-1, (1, 3)))]   # e6 = dx1dy1 - dx2dy2
C = sympy.Matrix(A6)
G_A6 = C * G_elem * C.T
out["1_T4_lattice"] = {
    "Gram_elementary_forms": [list(map(int, G_elem.row(i))) for i in range(6)],
    "det_elementary": int(G_elem.det()),
    "Gram_A6_combos_no_sqrt2": [int(G_A6[i, i]) for i in range(6)],
    "A6_Gram_offdiagonal_all_zero": all(G_A6[i, j] == 0 for i in range(6) for j in range(6) if i != j),
    "det_A6_Gram": int(G_A6.det()),
    "index_of_A6_span_in_H2(T4,Z)": abs(int(C.det())),
    "reading": "TT's (A.6) combos have Gram 2*eta33 but span an index-|det C| sublattice of H^2(T^4,Z) "
               "(unimodular, det -1); 'form a basis of H^2(T^4,Z)' and 'change of basis such that H33 = 2 eta33' hold over Q only.",
}
# ---------- (2) TT worked examples in both Grams ----------
eta_diag = sympy.diag(2, 2, 2, -2, -2, -2)
U = sympy.Matrix([[0, 1], [1, 0]]); H33 = sympy.diag(U, U, U)
def vec(**kw):
    v = sympy.zeros(6, 1)
    for k, c in kw.items(): v[int(k[1:]) - 1] = c
    return v
def dot(G, a, b): return (a.T * G * b)[0, 0]
ex = {}
for name, G in (("diag_2eta_(A.6-A.7)", eta_diag), ("hyperbolic_H33_(A.3)", H33)):
    ax, bx = 2 * vec(e1=1), 2 * vec(e2=1)                                   # (4.15)
    r415 = {"alpha_x^2": int(dot(G, ax, ax)), "beta_x^2": int(dot(G, bx, bx)), "alpha.beta": int(dot(G, ax, bx)),
            "N_D3_from_(4.14)": int(24 - dot(G, ax, ax))}
    ax, ay = 2 * vec(e1=1, e2=-1), 2 * vec(e1=1, e2=1, e4=1)                 # (4.31)
    bx, by = -4 * vec(e2=1), 2 * vec(e1=2, e4=1, e5=1)
    r431 = {"alpha_xx": int(dot(G, ax, ax)), "alpha_yy": int(dot(G, ay, ay)),
            "N_flux_(4.32)=ax.by-bx.ay": int(dot(G, ax, by) - dot(G, bx, ay))}
    bx5, ay5 = 2 * vec(e4=1), 2 * vec(e4=1)                                   # §5, lines 1185-1198
    r5 = {"N_flux_(2.8)=-bx.ay+by.ax": int(-dot(G, bx5, ay5))}
    ex[name] = {"(4.15)-(4.16)": r415, "(4.31)-(4.32)": r431, "sec5_lines1185-1198": r5}
out["2_TT_examples_by_Gram"] = ex
# ---------- (3) 8 | alpha^2 in the literal (unimodular) reading ----------
E8 = sympy.Matrix([[2, -1, 0, 0, 0, 0, 0, 0], [-1, 2, -1, 0, 0, 0, 0, 0], [0, -1, 2, -1, 0, 0, 0, -1],
                   [0, 0, -1, 2, -1, 0, 0, 0], [0, 0, 0, -1, 2, -1, 0, 0], [0, 0, 0, 0, -1, 2, -1, 0],
                   [0, 0, 0, 0, 0, -1, 2, 0], [0, 0, -1, 0, 0, 0, 0, 2]])       # typed from TT (A.4)
G319 = sympy.diag(H33, -E8, -E8)
v = sympy.Matrix(sympy.symbols("v1:23", integer=True))
q = sympy.expand((v.T * G319 * v)[0, 0])
half = sympy.expand(q / 2)
out["3_eight_divides_alpha_sq"] = {
    "det_E8_(A.4)": int(E8.det()), "E8_(A.4)_positive_definite": all(E8[:k, :k].det() > 0 for k in range(1, 9)),
    "det_Gamma319": int(G319.det()), "Gamma319_diagonal_all_even": all(G319[i, i] % 2 == 0 for i in range(22)),
    "v^2/2_has_integer_coefficients": all(c.is_integer for c in sympy.Poly(half, *v).coeffs()),
    "alpha=2v_=>_alpha^2=8*(integer_form)": sympy.expand(4 * q - 8 * half) == 0 and all(c.is_integer for c in sympy.Poly(half, *v).coeffs()),
}
vc = 2 * sympy.Matrix([1, 0, 0, 0, 0, 0])
nc = int((vc.T * sympy.diag(1, 1, 1, -1, -1, -1) * vc)[0, 0])
out["3_eight_divides_alpha_sq"]["negative_control_odd_lattice_I(3,3)"] = {
    "Gram": "diag(1,1,1,-1,-1,-1)", "v": [1, 0, 0, 0, 0, 0], "alpha=2v_alpha_sq": nc, "divisible_by_8": nc % 8 == 0}
# admissible alpha^2 <= 24 in 2*(U+U) realised with an orthogonal equal-norm partner (small search)
real = set()
for a in itertools.product(range(-2, 3), repeat=4):
    va = sympy.Matrix(a) * 2
    Gu = sympy.diag(U, U)
    n = (va.T * Gu * va)[0, 0]
    if 0 < n <= 24:
        for b in itertools.product(range(-2, 3), repeat=4):
            vb = sympy.Matrix(b) * 2
            if (vb.T * Gu * vb)[0, 0] == n and (va.T * Gu * vb)[0, 0] == 0:
                real.add(int(n)); break
out["3_eight_divides_alpha_sq"]["alpha_sq_realised_in_2(U+U)_coeff_|c|<=4"] = sorted(real)
# ---------- (4) tadpole arithmetic ----------
def ND3(Nflux): return Fr(24) - Fr(Nflux, 2)   # (2.3)
# (2.8) N_flux = -bx.ay + by.ax with (4.8) ay=-bx, by=ax; dot products as symbols (bilinear, symmetric)
aa, bb, ab = sympy.symbols("ax.ax bx.bx ax.bx")
Nflux_28 = -(-bb) + aa            # -bx.(-bx) + ax.ax
Nflux_413 = sympy.simplify(Nflux_28.subs(bb, aa))   # (4.9) bx^2 = ax^2
out["4_tadpole"] = {
    "(4.13)_N_flux_symbolic": str(Nflux_413),
    "(2.3)_N_D3_for_N_flux=16,32,48": [str(ND3(x)) for x in (16, 32, 48)],
    "line941_claim_N_flux=24_and_N_D3=0_consistent_with_(2.3)": ND3(24) == 0,
    "N_flux_required_for_N_D3=0": str(sympy.solve(sympy.Rational(24) - sympy.Symbol("Nf") / 2, sympy.Symbol("Nf"))[0]),
    "line1197_N_flux=8(sec5_example,diag Gram)_gives_N_D3": str(ND3(ex["diag_2eta_(A.6-A.7)"]["sec5_lines1185-1198"]["N_flux_(2.8)=-bx.ay+by.ax"])),
    "line1197_text_says": "N_flux/2 = 4, so that 20 D3-branes",
    "(4.32)_N_flux=32_gives_N_D3": str(ND3(32)),
    "(4.32)_text_says": "8 D3 branes",
}
# ---------- (5) E's orbit group vs orientation of the positive block ----------
cnt = rev = 0
for perm_pos in itertools.permutations(range(2)):
    for perm_neg in itertools.permutations(range(2)):
        for signs in itertools.product((1, -1), repeat=4):
            M = sympy.zeros(4, 4)
            for i, j in enumerate(perm_pos): M[j, i] = signs[i]
            for i, j in enumerate(perm_neg): M[2 + j, 2 + i] = signs[2 + i]
            G4 = sympy.diag(2, 2, -2, -2)
            assert M.T * G4 * M == G4
            cnt += 1
            if M[:2, :2].det() < 0: rev += 1
out["5_orbit_group"] = {"order": cnt, "elements_reversing_positive_block_orientation": rev,
                        "reading": "these are products with an odd number of reflections in norm +2 vectors; tier-L (from memory, "
                                   "Borcea/Donaldson): diffeomorphisms of K3 act through O^+(H^2), which preserves the orientation of "
                                   "positive-definite 3-planes. So E's group is not shown to lie in the physical identification group "
                                   "and its orbit counts are not certified upper bounds."}
here = os.path.dirname(os.path.abspath(__file__))
json.dump(out, open(os.path.join(here, "skE_tt_faithful_results.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
