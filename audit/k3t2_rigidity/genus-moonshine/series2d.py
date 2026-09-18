"""
Exact bivariate Laurent series arithmetic for Jacobi theta / eta q,y-series.

Grid convention (fixed, chosen so every object we need lands on an integer
lattice with NO floating point anywhere):
  - q-exponent index i  :  actual q-exponent = i / 8   (theta functions and
    eta^3, eta^6 always have q-exponents that are multiples of 1/8)
  - y-exponent index j  :  actual y-exponent = j / 2   (theta_1, theta_2 carry
    y^{1/2} factors; everything else has integer y powers = even j)

A series is a plain dict {(i, j): Fraction}. i, j are Python ints (can be
negative for j). Multiplication truncates at a caller-supplied I_MAX
(inclusive) in the i (q) index -- this is the only truncation in the whole
pipeline, and it is applied consistently by every product-building routine
below so every reported result is checked for truncation-stability by
rerunning with a larger I_MAX (see run_all.py).
"""
from fractions import Fraction as Fr
from collections import defaultdict


def const(c, i=0, j=0):
    return {(i, j): Fr(c)}


def add(a, b):
    out = defaultdict(Fr)
    for k, v in a.items():
        out[k] += v
    for k, v in b.items():
        out[k] += v
    return {k: v for k, v in out.items() if v != 0}


def scal(a, c):
    c = Fr(c)
    if c == 0:
        return {}
    return {k: v * c for k, v in a.items()}


def mul(a, b, imax):
    out = defaultdict(Fr)
    for (i1, j1), v1 in a.items():
        if i1 > imax:
            continue
        for (i2, j2), v2 in b.items():
            i = i1 + i2
            if i > imax:
                continue
            out[(i, j1 + j2)] += v1 * v2
    return {k: v for k, v in out.items() if v != 0}


def mul_many(factors, imax):
    out = const(1)
    for f in factors:
        out = mul(out, f, imax)
    return out


def shift(a, di=0, dj=0):
    return {(i + di, j + dj): v for (i, j), v in a.items()}


def truncate(a, imax):
    return {k: v for k, v in a.items() if k[0] <= imax}


def y_to_1(a):
    """Substitute y=1 (sum over j for each i). Returns 1-D dict {i: Fraction}."""
    out = defaultdict(Fr)
    for (i, j), v in a.items():
        out[i] += v
    return {k: v for k, v in out.items() if v != 0}


def divide_scalar_series(num, den1d, imax):
    """
    Divide a 2D series `num` {(i,j):Fr} by a 1D (y-independent, i.e. j=0
    only) series `den1d` {i:Fr}, both living on q-index step `istep`
    (the spacing between consecutive nonzero i values in den1d -- e.g. 1 for
    theta_2/theta_2(0) [q-integer steps], 4 for theta_3/theta_3(0) or
    theta_4/theta_4(0) [q-half-integer steps, i.e. i steps of 4 since i=8q]).
    Requires den1d's minimal i (i0) to have a nonzero coefficient, and
    num's minimal i to be >= i0 and on the same residue class mod istep.
    Returns the quotient as a 2D dict, exact (Fraction arithmetic), valid
    up to i<=imax.
    """
    istep = _min_step(den1d)
    i0 = min(den1d)
    d0 = den1d[i0]
    # D(q) = q^{i0/8} * d0 * sum_k dlist[k] Q^k,  Q := q^{istep/8},  dlist[0]=1
    kmax = imax // istep
    dlist = [den1d.get(i0 + k * istep, Fr(0)) / d0 for k in range(kmax + 1)]
    # N(q) = q^{i0/8} * sum_k numg[k] Q^k  (numerator must start at the same
    # i0 -- i.e. num and den have the same leading q-power, as is the case
    # for every theta_i(z)/theta_i(0) ratio we build).
    numg = defaultdict(lambda: defaultdict(Fr))
    for (i, j), v in num.items():
        assert (i - i0) % istep == 0, f"numerator off q-grid: i={i}, i0={i0}, istep={istep}"
        k = (i - i0) // istep
        if k <= kmax:
            numg[k][j] += v
    # N/D = (1/d0) * [sum numg Q^k] / [sum dlist Q^k], a series in Q starting
    # at Q^0 (the q^{i0/8} factors cancel between N and D).
    # Let T(Q) := d0 * (N/D)(Q). Then T(Q)*dlist(Q) = numg(Q), so
    #   T[k] = numg[k] - sum_{m<k} T[m]*dlist[k-m]   (dlist[0]=1),
    # and the actual quotient is S[k] = T[k]/d0. The recursion MUST use the
    # undivided T[] (not S[]) on the right-hand side, or per-step division
    # compounds incorrectly.
    Tlist = {}
    for k in range(kmax + 1):
        rhs = defaultdict(Fr)
        for j, v in numg.get(k, {}).items():
            rhs[j] += v
        for m in range(k):
            tk = Tlist.get(m, {})
            dcoef = dlist[k - m]
            if dcoef == 0 or not tk:
                continue
            for j, v in tk.items():
                rhs[j] -= v * dcoef
        assert dlist[0] == 1
        Tlist[k] = {j: v for j, v in rhs.items() if v != 0}
    out = {}
    for k, jd in Tlist.items():
        i = k * istep
        for j, v in jd.items():
            s = v / d0
            if s != 0:
                out[(i, j)] = s
    return out


def reciprocal_1d(den1d, imax):
    """1/D(q) for a plain (y-independent) series D = {i: Fr}, D's minimal
    exponent i0 has nonzero coefficient d0. Returns {i: Fr} representing
    the reciprocal series (its minimal exponent is -i0; the result can and
    generally will have negative i keys -- that's fine, series2d.mul()
    handles negative i keys correctly, it only ever drops i > imax).
    Unlike divide_scalar_series, this places NO constraint on what a
    subsequent numerator's residue class (mod istep) must be: multiply the
    result by any 2D series via mul() to get an exact quotient.
    """
    istep = _min_step(den1d)
    i0 = min(den1d)
    d0 = den1d[i0]
    kmax = (imax + i0) // istep
    if kmax < 0:
        return {}
    dlist = [den1d.get(i0 + k * istep, Fr(0)) / d0 for k in range(kmax + 1)]
    assert dlist[0] == 1
    Tlist = [Fr(1)]
    for k in range(1, kmax + 1):
        acc = Fr(0)
        for m in range(k):
            if dlist[k - m] != 0:
                acc += Tlist[m] * dlist[k - m]
        Tlist.append(-acc)
    out = {}
    for k, t in enumerate(Tlist):
        s = t / d0
        if s != 0:
            out[-i0 + k * istep] = s
    return out


def _min_step(den1d):
    idxs = sorted(den1d)
    if len(idxs) < 2:
        return 1
    from math import gcd
    g = 0
    base = idxs[0]
    for x in idxs[1:]:
        g = gcd(g, x - base)
    return g if g != 0 else 1


def pretty(a, qden=8, yden=2, var='q', yvar='y'):
    """Human-readable string, sorted by i then j."""
    terms = []
    for (i, j) in sorted(a.keys()):
        v = a[(i, j)]
        if v == 0:
            continue
        qexp = Fr(i, qden)
        yexp = Fr(j, yden)
        terms.append(f"({v})*{var}^{qexp}*{yvar}^{yexp}")
    return " + ".join(terms) if terms else "0"
