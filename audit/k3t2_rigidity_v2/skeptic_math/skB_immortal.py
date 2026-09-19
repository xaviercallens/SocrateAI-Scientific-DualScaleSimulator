"""Skeptic (math lens) independent re-derivation for Track B.
Checks, with own product-formula series (no Track B code or cache):
  (i)  solve B^2 = a*A*A21 + b*E4*A^2 + c*A*Hhat exactly (overdetermined, unknowns free) and
       translate to Track B's G_2/A = 9B^2/(4A)+3E4A/4 normalisation => implied (N, E4A coeff, M).  A=phi_{-2,1}, B=phi_{0,1}, Hhat = sum H(4n-l^2) q^n y^l with H(0)=-1/12,
       A_{2,1} = sum_s q^{s^2+s} y^{2s+1}/(1-q^s y)^2.  The s=0 term y/(1-y)^2 is multiplied by A
       EXACTLY (each q^n coefficient of A is divisible by (1-y)^2), so no region/strip convention
       enters; s!=0 terms are expanded in q^{|s|} (valid for |q|<|y|<|q|^-1).
  (ii) negative control: H(3) -> H(3)+1/6 makes the system inconsistent.
  (iii) Goettsche prod(1-q^n)^-24 coefficients 0..8.
Command: cd audit/k3t2_rigidity_v2/skeptic_math && <venv python> skB_immortal.py
"""
import json, os, math
from fractions import Fraction as Fr
from collections import defaultdict
QMAX = 5

def mul(a, b):
    c = defaultdict(Fr)
    for (n1, l1), x in a.items():
        for (n2, l2), y in b.items():
            if n1 + n2 <= QMAX:
                c[(n1 + n2, l1 + l2)] += x * y
    return {k: v for k, v in c.items() if v != 0}

def add(*ss, coeffs=None):
    c = defaultdict(Fr)
    for i, s in enumerate(ss):
        w = coeffs[i] if coeffs else 1
        for k, v in s.items():
            c[k] += w * v
    return {k: v for k, v in c.items() if v != 0}

one = {(0, 0): Fr(1)}
def factor(n, lpow, sign):  # (1 + sign q^n y^lpow)
    return {(0, 0): Fr(1), (n, lpow): Fr(sign)}

def inv_factor(n, sign):  # 1/(1 + sign q^n), q-only
    return {(n * k, 0): Fr((-sign) ** k) for k in range(QMAX // n + 1)}

# A = phi_{-2,1} = (y - 2 + 1/y) prod (1-q^n y)^2 (1-q^n/y)^2 / (1-q^n)^4
Aprod = dict(one)
for n in range(1, QMAX + 1):
    for f in (factor(n, 1, -1), factor(n, 1, -1), factor(n, -1, -1), factor(n, -1, -1)):
        Aprod = mul(Aprod, f)
    for _ in range(4):
        Aprod = mul(Aprod, inv_factor(n, -1))
A = mul({(0, 1): Fr(1), (0, 0): Fr(-2), (0, -1): Fr(1)}, Aprod)

# B = phi_{0,1} = 4 sum_i (theta_i(z)/theta_i(0))^2, product forms in q^{1/2}: use n in half units
# To keep integer q-exponents, work with variable Q2 = q^{1/2} then keep integer powers.
QM2 = 2 * QMAX
def mul2(a, b):
    c = defaultdict(Fr)
    for (n1, l1), x in a.items():
        for (n2, l2), y in b.items():
            if n1 + n2 <= QM2:
                c[(n1 + n2, l1 + l2)] += x * y
    return {k: v for k, v in c.items() if v != 0}
def inv2(n, sign):
    return {(n * k, 0): Fr((-sign) ** k) for k in range(QM2 // n + 1)}
def ratio_sq(half_shift, sign):
    # prod_{n>=1} (1 + sign q^{m} y)(1 + sign q^{m}/y)/(1 + sign q^{m})^2, m = n - half_shift/2, squared
    r = {(0, 0): Fr(1)}
    for n in range(1, QMAX + 2):
        e = 2 * n - half_shift  # exponent in q^{1/2} units
        if e > QM2: break
        for _ in range(2):
            r = mul2(r, {(0, 0): Fr(1), (e, 1): Fr(sign)})
            r = mul2(r, {(0, 0): Fr(1), (e, -1): Fr(sign)})
            r = mul2(r, inv2(e, sign)); r = mul2(r, inv2(e, sign))
    return r
t2 = mul2({(0, 1): Fr(1, 4), (0, 0): Fr(2, 4), (0, -1): Fr(1, 4)}, ratio_sq(0, +1))
t3 = ratio_sq(1, +1)
t4 = ratio_sq(1, -1)
B2 = add(t2, t3, t4, coeffs=[4, 4, 4])
assert all(n % 2 == 0 for (n, l) in B2), "half-integer q powers must cancel"
B = {(n // 2, l): v for (n, l), v in B2.items()}

# Hurwitz class numbers by reduced forms (own implementation)
def hurwitz(D):
    if D == 0: return Fr(-1, 12)
    if D < 0 or D % 4 in (1, 2): return Fr(0)
    tot = Fr(0)
    b = D % 2
    while b * b <= D // 3 + 1:
        if (b * b + D) % 4 == 0:
            ac = (b * b + D) // 4
            a = max(b, 1)
            while a * a <= ac:
                if ac % a == 0:
                    c = ac // a
                    if a <= c and b <= a:
                        w = Fr(1)
                        if a == b == c: w = Fr(1, 3)
                        elif a == c and b == 0: w = Fr(1, 2)
                        # count (a,b,c) and (a,-b,c) unless b==0 or b==a or a==c
                        mult = 1 if (b == 0 or b == a or a == c) else 2
                        tot += w * mult
                a += 1
        b += 2
    return tot
Hhat = {}
for n in range(0, QMAX + 1):
    for l in range(-2 * int(math.isqrt(4 * n)) - 1, 2 * int(math.isqrt(4 * n)) + 2):
        v = hurwitz(4 * n - l * l)
        if v: Hhat[(n, l)] = v

# A * A_{2,1}
def poly_div_1my_sq(coeffs_by_l):
    # divide Laurent poly sum c_l y^l by (1-y)^2, exact; return dict
    if not coeffs_by_l: return {}
    lo, hi = min(coeffs_by_l), max(coeffs_by_l)
    p = [coeffs_by_l.get(l, Fr(0)) for l in range(lo, hi + 1)]
    # divide by 1 - 2y + y^2 (ascending)
    qt = []
    rem = p[:]
    for i in range(len(p) - 2):
        c = rem[i]; qt.append(c)
        rem[i] -= c; rem[i + 1] += 2 * c; rem[i + 2] -= c
    assert all(x == 0 for x in rem), "not divisible by (1-y)^2"
    return {lo + i: c for i, c in enumerate(qt) if c != 0}
AA21 = {}
# s = 0: A * y/(1-y)^2
for n in range(0, QMAX + 1):
    row = {l: v for (m, l), v in A.items() if m == n}
    for l, v in poly_div_1my_sq(row).items():
        AA21[(n, l + 1)] = AA21.get((n, l + 1), Fr(0)) + v
A21rest = defaultdict(Fr)
for s in list(range(1, QMAX + 2)) + list(range(-QMAX - 2, 0)):
    for k in range(0, 4 * QMAX + 5):
        if s >= 1:
            n, l = s * s + s + s * k, 2 * s + 1 + k
        else:
            n, l = s * s - s - s * k, 2 * s - 1 - k
        if n <= QMAX: A21rest[(n, l)] += k + 1
AA21 = add(AA21, mul(A, dict(A21rest)))

BB = mul(B, B)
AH = mul(A, Hhat)
E4 = {(m, 0): Fr(240 * sum(d ** 3 for d in range(1, m + 1) if m % d == 0)) for m in range(1, QMAX + 1)}
E4[(0, 0)] = Fr(1)
E4AA = mul(E4, mul(A, A))
# Solve B^2 = a*(A*A21) + b*(E4*A^2) + c*(A*Hhat) exactly (overdetermined linear system, no target typed).
import sympy
a_, b_, c_ = sympy.symbols("a b c")
keys = set(BB) | set(AA21) | set(E4AA) | set(AH)
eqs = [sympy.Rational(BB.get(k, 0)) - a_ * sympy.Rational(AA21.get(k, 0)) - b_ * sympy.Rational(E4AA.get(k, 0))
       - c_ * sympy.Rational(AH.get(k, 0)) for k in keys]
sol = sympy.solve(eqs, [a_, b_, c_], dict=True)
sol = {str(k): str(v) for k, v in sol[0].items()} if sol else None
# translate to Track B's normalisation: G2/A = 9B^2/(4A) + 3E4A/4 ; G2/A - N A21 = 3E4A - M Hhat
# => 9/4 * (a A21 + b E4 A + c Hhat) + 3/4 E4 A - N A21 = 3 E4 A - M Hhat
if sol:
    a, b, c = (sympy.Rational(sol[x]) for x in "abc")
    implied = {"N": str(sympy.Rational(9, 4) * a), "E4A_coefficient": str(sympy.Rational(9, 4) * b + sympy.Rational(3, 4)),
               "M": str(-sympy.Rational(9, 4) * c)}
else:
    implied = None
# negative control: perturb Hhat's H(3) by +1/6 -> system must become inconsistent
Hbad = dict(Hhat)
for k in list(Hbad):
    if 4 * k[0] - k[1] ** 2 == 3: Hbad[k] += Fr(1, 6)
AHb = mul(A, Hbad)
eqsb = [sympy.Rational(BB.get(k, 0)) - a_ * sympy.Rational(AA21.get(k, 0)) - b_ * sympy.Rational(E4AA.get(k, 0))
        - c_ * sympy.Rational(AHb.get(k, 0)) for k in set(keys) | set(AHb)]
solb = sympy.solve(eqsb, [a_, b_, c_], dict=True)
target = None
scan = None
# Goettsche
g = [Fr(1)] + [Fr(0)] * 8
for n in range(1, 9):
    for _ in range(24):
        for i in range(n, 9):
            g[i] += g[i - n]
out = {"QMAX": QMAX,
       "B_q0": {str(l): str(v) for (n, l), v in B.items() if n == 0},
       "A_q0": {str(l): str(v) for (n, l), v in A.items() if n == 0},
       "H_first": {D: str(hurwitz(D)) for D in (0, 3, 4, 7, 8, 11, 12, 15, 16, 19, 20, 23, 24)},
       "exact_solution_B2_eq_a_AA21_plus_b_E4A2_plus_c_AHhat": sol,
       "n_equations": len(keys),
       "implied_trackB_constants": implied,
       "negative_control_H3_plus_1_6_solutions": [str(x) for x in solb],
       "goettsche_0_8": [str(x) for x in g]}
json.dump(out, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "skB_immortal_results.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
