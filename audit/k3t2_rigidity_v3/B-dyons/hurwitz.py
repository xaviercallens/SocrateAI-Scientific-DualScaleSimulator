"""
Hurwitz class numbers H(D), D <= DMAX, by (1) exact reduced-form counting and (2) an INDEPENDENT
Dirichlet-L-function route, plus the Kronecker-Hurwitz relation on the whole table.

Run: cd audit/k3t2_rigidity_v3/B-dyons && <venv-python> hurwitz.py 400
(committed run: DMAX = 400)  Writes hurwitz_results.json.

(1) H(D) = weighted number of reduced positive-definite forms (a,b,c), b^2-4ac=-D, imprimitive forms
    included; -a<b<=a<=c, b>=0 if a==c; weights 1/2 for a=c,b=0 and 1/3 for a=b=c. H(0)=-1/12 (declared input).
(2) Dirichlet route: write -D = f^2 d0, d0 the fundamental discriminant. For a fundamental d<-4
    h(d) = -(1/|d|) sum_{a=1}^{|d|-1} chi_d(a) a (chi_d Kronecker symbol); h(-3)=h(-4)=1. The weighted
    class number of the order of conductor g is h_w(g^2 d0) = h(d0)/(w0/2) * g * prod_{p|g}(1-chi_d0(p)/p)
    (w0 = 6, 4, 2 for d0=-3,-4, other) and H(D) = sum_{g | f} h_w(g^2 d0). No reduced form is counted here.
Both are exact (fractions). Hurwitz counting and Dirichlet route are compared for EVERY D<=DMAX with D=0,3 mod 4.
Negative control: a wrong weight (1 instead of 1/2 for a=c,b=0 forms) changes counted H and breaks both the
Dirichlet comparison and the Kronecker-Hurwitz relation.
"""
import json, sys
from fractions import Fraction as Fr
from sympy import factorint, jacobi_symbol
from common import HERE


def H_count(D, w_ac=Fr(1, 2)):
    if D == 0:
        return Fr(-1, 12)
    if D % 4 not in (0, 3):
        return Fr(0)
    tot = Fr(0)
    amax = int((D / 3) ** 0.5) + 2
    for a in range(1, amax + 1):
        for b in range(-a + 1, a + 1):
            num = b * b + D
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a or (c == a and b < 0):
                continue
            if b == 0 and a == c:
                tot += w_ac
            elif a == b == c:
                tot += Fr(1, 3)
            else:
                tot += 1
    return tot


def kronecker(d, a):
    """Kronecker symbol (d/a), a >= 1 integer, d a discriminant."""
    res = 1
    while a % 2 == 0:
        a //= 2
        if d % 2 == 0:
            return 0
        res *= 1 if d % 8 in (1, 7) else -1
    if a == 1:
        return res
    return res * jacobi_symbol(d % a, a)


def fund_disc_and_conductor(D):
    """-D = f^2 * d0, d0 fundamental discriminant (<0)."""
    d = -D
    best = None
    for f in range(1, int(D ** 0.5) + 2):
        if D % (f * f):
            continue
        d0 = d // (f * f)
        if d0 % 4 not in (0, 1):
            continue
        # fundamental test
        ok = False
        if d0 % 4 == 1:
            ok = all(e == 1 for e in factorint(-d0).values()) if -d0 > 1 else False
        else:
            m = d0 // 4
            if m % 4 in (2, 3):
                ok = all(e == 1 for e in factorint(-m).values())
        if ok:
            best = (d0, f)
    return best


def h_dirichlet_fund(d0):
    if d0 == -3 or d0 == -4:
        return Fr(1)
    s = sum(kronecker(d0, a) * a for a in range(1, -d0))
    return Fr(-s, -d0)


def w_of(d0):
    return 6 if d0 == -3 else 4 if d0 == -4 else 2


def h_w_order(d0, g):
    h0 = h_dirichlet_fund(d0)
    val = h0 / Fr(w_of(d0), 2) * g
    for p in factorint(g):
        val *= 1 - Fr(kronecker(d0, p), p)
    return val


def H_dirichlet(D):
    if D == 0:
        return Fr(-1, 12)
    if D % 4 not in (0, 3):
        return Fr(0)
    d0, f = fund_disc_and_conductor(D)
    return sum(h_w_order(d0, g) for g in range(1, f + 1) if f % g == 0)


def sigma1(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)


def min_sum(n):
    return sum(min(d, n // d) for d in range(1, n + 1) if n % d == 0)


if __name__ == "__main__":
    DMAX = int(sys.argv[1])
    tab, tabD, mism = {}, {}, []
    for D in range(0, DMAX + 1):
        if D != 0 and D % 4 not in (0, 3):
            continue
        tab[D] = H_count(D)
        tabD[D] = H_dirichlet(D)
        if tab[D] != tabD[D]:
            mism.append(D)
    # Kronecker-Hurwitz on whole table
    NMAX = DMAX // 4
    kh_fail = []
    for n in range(1, NMAX + 1):
        s = Fr(0)
        for t in range(-int((4 * n) ** 0.5) - 1, int((4 * n) ** 0.5) + 2):
            D = 4 * n - t * t
            if D >= 0:
                s += tab.get(D, Fr(0)) if D in tab else Fr(0)
        if s != 2 * sigma1(n) - min_sum(n):
            kh_fail.append(n)
    # negative control: wrong weight for a=c,b=0
    tab_bad = {D: H_count(D, w_ac=Fr(1)) for D in tab}
    n_bad_vs_dirichlet = sum(1 for D in tab if tab_bad[D] != tabD[D])
    kh_bad = 0
    for n in range(1, NMAX + 1):
        s = Fr(0)
        for t in range(-int((4 * n) ** 0.5) - 1, int((4 * n) ** 0.5) + 2):
            D = 4 * n - t * t
            if D >= 0 and D in tab_bad:
                s += tab_bad[D]
        if s != 2 * sigma1(n) - min_sum(n):
            kh_bad += 1
    # how many D are fundamental (pure Dirichlet-class-number-formula check)
    n_fund = 0
    for D in tab:
        if D > 4:
            d0, f = fund_disc_and_conductor(D)
            if f == 1:
                n_fund += 1
    out = {
        "DMAX": DMAX, "num_D_values": len(tab),
        "H_table": {str(D): str(v) for D, v in tab.items()},
        "counting_vs_dirichlet_mismatch_D": mism,
        "num_fundamental_discriminants_checked_with_pure_class_number_formula": n_fund,
        "kronecker_hurwitz_n_range": [1, NMAX], "kronecker_hurwitz_failures": kh_fail,
        "negative_control_wrong_weight_1_for_a_eq_c_b_eq_0": {
            "num_D_where_counting_differs_from_dirichlet": n_bad_vs_dirichlet,
            "num_n_where_kronecker_hurwitz_breaks": kh_bad},
        "input_H0": "-1/12 (declared input Hurwitz_H0)",
        "verified_by_exact_computation": (not mism) and (not kh_fail) and n_bad_vs_dirichlet > 0 and kh_bad > 0,
    }
    with open(HERE / "hurwitz_results.json", "w") as f:
        json.dump(out, f, indent=1)
    print({k: out[k] for k in ("DMAX", "num_D_values", "counting_vs_dirichlet_mismatch_D",
          "num_fundamental_discriminants_checked_with_pure_class_number_formula",
          "kronecker_hurwitz_failures", "negative_control_wrong_weight_1_for_a_eq_c_b_eq_0", "verified_by_exact_computation")})
