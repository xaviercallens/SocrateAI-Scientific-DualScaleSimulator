"""
Independent (n, l)-keyed exact Laurent-series engine: coefficient of q^n y^l,
n >= 0 (power series in q), l any integer (Laurent in y).

This is a FRESH implementation (not a copy of LeanMaster's Lean code, not a
copy of the sibling dyons/series.py P,Y-power engine) built directly from the
defining formulas quoted in DualScaleDyons/DMVV.lean and
DualScaleDyons/ImmortalHigher.lean's doc comments (which are themselves
citations of Dabholkar-Murthy-Zagier arXiv:1208.4074, Tier L). All arithmetic
is exact (fractions.Fraction / Python int).
"""
from fractions import Fraction as Fr
from math import comb


def clean(s):
    return {k: v for k, v in s.items() if v != 0}


def add(*ss):
    out = {}
    for s in ss:
        for k, v in s.items():
            out[k] = out.get(k, Fr(0)) + v
    return clean(out)


def sub(s1, s2):
    return add(s1, scal(-1, s2))


def scal(c, s):
    c = Fr(c)
    if c == 0:
        return {}
    return clean({k: c * v for k, v in s.items()})


def mul(s1, s2, qmax):
    out = {}
    for (n1, l1), v1 in s1.items():
        if n1 > qmax:
            continue
        for (n2, l2), v2 in s2.items():
            n = n1 + n2
            if n > qmax:
                continue
            k = (n, l1 + l2)
            out[k] = out.get(k, Fr(0)) + v1 * v2
    return clean(out)


def mul_many(ss, qmax):
    out = {(0, 0): Fr(1)}
    for s in ss:
        out = mul(out, s, qmax)
    return out


def powr(s, e, qmax):
    out = {(0, 0): Fr(1)}
    for _ in range(e):
        out = mul(out, s, qmax)
    return out


def truncate(s, qmax):
    return clean({k: v for k, v in s.items() if k[0] <= qmax})


def from_1d(d1d):
    """A y-independent q-series {n: Fraction} -> (n,l) dict with l=0."""
    return clean({(n, 0): v for n, v in d1d.items()})


def one_minus_yq_pow(qmax, y_l, coeff_sign, k):
    """
    Series for (1 + coeff_sign * q^1 y^y_l)^k truncated at q^qmax, k any integer
    (only the q^1 y^{y_l} monomial is used per factor here; caller supplies n
    via repeated calls at shifted base, see one_minus_pow below for the general
    (1 - c q^n y^l)^k builder actually used).
    """
    raise NotImplementedError


def binom_pow_monomial(qmax, n, l, coeff, k):
    """
    Series for (1 + coeff * q^n y^l)^k, k any integer (pos or neg), truncated
    at q^qmax. Requires n > 0 so the q-series terminates/converges formally.
    """
    assert n > 0
    out = {(0, 0): Fr(1)}
    jmax = qmax // n
    for j in range(1, jmax + 1):
        qn = n * j
        if qn > qmax:
            break
        if k >= 0:
            if j > k:
                continue
            b = Fr(comb(k, j))
        else:
            b = Fr((-1) ** j * comb(-k + j - 1, j))
        out[(qn, l * j)] = out.get((qn, l * j), Fr(0)) + b * (Fr(coeff) ** j)
    return clean(out)


def neg_binom_coeff(e, j):
    """Coefficient of x^j in (1-x)^{-e}, e any integer (pos/neg/zero), j>=0.
    Exact: e(e+1)...(e+j-1)/j! (rising factorial over j!), matching Lean's
    negBinom. Always an integer for integer e (checked, not assumed)."""
    if j == 0:
        return 1
    num = 1
    for i in range(j):
        num *= (e + i)
    den = 1
    for i in range(1, j + 1):
        den *= i
    q, r = divmod(num, den)
    assert r == 0, f"negBinom({e},{j}) not exact: {num}/{den}"
    return q
