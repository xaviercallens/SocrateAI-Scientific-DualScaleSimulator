"""
Exact Laurent-series engine for Jacobi theta functions in the K3 x T^2
quarter-BPS dyon computation (DMVV / Dabholkar-Murthy-Zagier).

Convention (standard string-theory convention, q = e^{2 pi i tau}, y = e^{2 pi i z}):

theta1(tau,z) = -i q^{1/8} (y^{1/2}-y^{-1/2}) prod_{n>=1} (1-q^n)(1-q^n y)(1-q^n/y)
theta2(tau,z) =    q^{1/8} (y^{1/2}+y^{-1/2}) prod_{n>=1} (1-q^n)(1+q^n y)(1+q^n/y)
theta3(tau,z) =                               prod_{n>=1} (1-q^n)(1+q^{n-1/2} y)(1+q^{n-1/2}/y)
theta4(tau,z) =                               prod_{n>=1} (1-q^n)(1-q^{n-1/2} y)(1-q^{n-1/2}/y)

To avoid fractional exponents entirely, every series is carried in the
integer-lattice variables

    P := q^{1/8}   (so q = P^8, q^{1/8} = P, q^{n-1/2} = P^{8n-4})
    Y := y^{1/2}   (so y = Y^2, y^{1/2} = Y)

A Series is a dict {(a, b): Fraction} meaning coefficient of P^a * Y^b.
All products below are truncated at a maximum P-power PMAX (= 8 * QMAX,
QMAX being the truncation order in the true variable q) and a maximum
|Y|-power YMAX chosen generously (Y-power growth is bounded by the number
of factors touched, itself bounded by PMAX).
"""
from fractions import Fraction as Fr
from itertools import product as _iproduct

def const(c=1):
    return {(0, 0): Fr(c)}

def add(*series_list):
    out = {}
    for s in series_list:
        for k, v in s.items():
            out[k] = out.get(k, Fr(0)) + v
    return {k: v for k, v in out.items() if v != 0}

def scal(c, s):
    c = Fr(c)
    if c == 0:
        return {}
    return {k: c * v for k, v in s.items() if c * v != 0}

def mul(s1, s2, pmax, ymax=None):
    out = {}
    for (a1, b1), v1 in s1.items():
        if a1 > pmax:
            continue
        for (a2, b2), v2 in s2.items():
            a = a1 + a2
            if a > pmax:
                continue
            b = b1 + b2
            if ymax is not None and abs(b) > ymax:
                continue
            key = (a, b)
            out[key] = out.get(key, Fr(0)) + v1 * v2
    return {k: v for k, v in out.items() if v != 0}

def mul_many(series_iterable, pmax, ymax=None):
    out = const(1)
    for s in series_iterable:
        out = mul(out, s, pmax, ymax)
    return out

def truncate(s, pmax, ymax=None):
    out = {}
    for (a, b), v in s.items():
        if a > pmax:
            continue
        if ymax is not None and abs(b) > ymax:
            continue
        if v != 0:
            out[(a, b)] = v
    return out

def one_plus_x_pow(a, b, coeff, k, pmax, ymax=None):
    """
    Series for (1 + coeff * P^a * Y^b)^k, k any integer, truncated at P^pmax.
    Requires a > 0 (so the series in P terminates naturally within pmax).
    """
    assert a > 0, "one_plus_x_pow requires positive P-valuation for truncation to terminate"
    out = {(0, 0): Fr(1)}
    nmax = pmax // a
    # binomial coefficients C(k, n) for possibly negative integer k
    from math import comb
    for n in range(1, nmax + 1):
        pa = a * n
        if pa > pmax:
            break
        yb = b * n
        if ymax is not None and abs(yb) > ymax:
            continue
        if k >= 0:
            if n > k:
                continue
            binom = Fr(comb(k, n))
        else:
            # C(k,n) for negative k: (-1)^n * C(-k+n-1, n)
            binom = Fr((-1) ** n * comb(-k + n - 1, n))
        out[(pa, yb)] = out.get((pa, yb), Fr(0)) + binom * (coeff ** n)
    return {kk: v for kk, v in out.items() if v != 0}

def y_pow_to_l(s):
    """Convert a series keyed by (P-power a, Y-power b) with a multiple of 8
    and b even into a series keyed by (q-power n, y-power l) = (a/8, b/2)."""
    out = {}
    for (a, b), v in s.items():
        assert a % 8 == 0, f"non-integer q power encountered: P^{a}"
        assert b % 2 == 0, f"non-integer y power encountered: Y^{b}"
        out[(a // 8, b // 2)] = out.get((a // 8, b // 2), Fr(0)) + v
    return {k: v for k, v in out.items() if v != 0}


def eta24_over_q(qmax):
    """prod_{n>=1} (1-q^n)^24 as a series in q up to q^qmax (Fraction dict keyed by n)."""
    from math import comb
    s = {0: Fr(1)}
    for n in range(1, qmax + 1):
        # multiply by (1-q^n)^24
        factor = {}
        for j in range(0, 25):
            e = n * j
            if e > qmax:
                break
            factor[e] = Fr((-1) ** j * comb(24, j))
        new_s = {}
        for e1, v1 in s.items():
            for e2, v2 in factor.items():
                e = e1 + e2
                if e > qmax:
                    continue
                new_s[e] = new_s.get(e, Fr(0)) + v1 * v2
        s = new_s
    return s


def prod_1_minus_qn_power(qmax, power):
    """prod_{n>=1} (1-q^n)^power (power can be negative), q-power dict up to qmax."""
    from math import comb
    s = {0: Fr(1)}
    for n in range(1, qmax + 1):
        factor = {}
        nmax = qmax // n
        for j in range(0, nmax + 1):
            e = n * j
            if e > qmax:
                break
            if power >= 0:
                if j > power:
                    continue
                c = Fr(comb(power, j)) * (-1) ** j
            else:
                c = Fr((-1) ** j * comb(-power + j - 1, j))
            factor[e] = factor.get(e, Fr(0)) + c
        new_s = {}
        for e1, v1 in s.items():
            for e2, v2 in factor.items():
                e = e1 + e2
                if e > qmax:
                    continue
                new_s[e] = new_s.get(e, Fr(0)) + v1 * v2
        s = new_s
    return {k: v for k, v in s.items() if v != 0}
