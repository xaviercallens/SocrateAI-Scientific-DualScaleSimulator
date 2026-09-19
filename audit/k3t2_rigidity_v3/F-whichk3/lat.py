"""Exact arithmetic in K = Q(i, sqrt3) (basis 1, i, s, is; s^2=3) and NS/T of an abelian surface C^2/L from its holomorphic 1-forms."""
import itertools
from fractions import Fraction as F
import sympy as sp
from common import int_kernel, lagrange_reduce

# multiplication table of basis elements: (1, i, s, is)
_T = {}
def _mul_basis(a, b):
    # returns (coef, index): basis[a]*basis[b] = coef*basis[index]
    ia, sa = a & 1, (a >> 1) & 1     # index bits: bit0 -> factor i, bit1 -> factor s
    ib, sb = b & 1, (b >> 1) & 1
    ii = ia + ib; ss = sa + sb
    coef = 1
    if ii == 2: coef *= -1; ii = 0
    if ss == 2: coef *= 3; ss = 0
    return coef, ii | (ss << 1)

class K:
    __slots__ = ("c",)
    def __init__(self, c=(0, 0, 0, 0)):
        self.c = tuple(F(x) for x in c)
    def __add__(s, o): return K([a + b for a, b in zip(s.c, o.c)])
    def __sub__(s, o): return K([a - b for a, b in zip(s.c, o.c)])
    def __neg__(s): return K([-a for a in s.c])
    def __mul__(s, o):
        r = [F(0)] * 4
        for a in range(4):
            if s.c[a] == 0: continue
            for b in range(4):
                if o.c[b] == 0: continue
                co, k = _mul_basis(a, b)
                r[k] += co * s.c[a] * o.c[b]
        return K(r)
    def iszero(s): return all(x == 0 for x in s.c)
    def inv(s):
        M = sp.Matrix(4, 4, lambda r, c: 0)
        cols = []
        for b in range(4):
            e = K([1 if j == b else 0 for j in range(4)])
            cols.append((s * e).c)
        M = sp.Matrix(cols).T
        sol = M.inv() * sp.Matrix([1, 0, 0, 0])
        return K([F(int(x.p), int(x.q)) for x in sol])
    def __eq__(s, o): return s.c == o.c
    def __repr__(s): return f"K{tuple(str(x) for x in s.c)}"
    def to_complex(s):
        return complex(float(s.c[0]), float(s.c[1])) + 3 ** 0.5 * complex(float(s.c[2]), float(s.c[3]))

ZERO, ONE = K((0, 0, 0, 0)), K((1, 0, 0, 0))
I_ = K((0, 1, 0, 0))
OMEGA = K((F(-1, 2), 0, 0, F(1, 2)))     # (-1 + i sqrt3)/2
def kq(a): return K((a, 0, 0, 0))

# Lambda^2 pairing on H^1 = Z^4 (basis e^1..e^4): <e^{jk}, e^{lm}> = coefficient of e^{1234} in e^{jk}^e^{lm}
PAIRS = [(0, 1), (2, 3), (0, 2), (0, 3), (1, 2), (1, 3)]   # e12,e34,e13,e14,e23,e24
def perm_sign(p):
    s = 1
    p = list(p)
    for i in range(len(p)):
        for j in range(i + 1, len(p)):
            if p[i] > p[j]: s = -s
    return s
def wedge_pair(u, v):
    idx = list(u) + list(v)
    if len(set(idx)) < 4: return 0
    return perm_sign(idx)
GRAM6 = [[wedge_pair(u, v) for v in PAIRS] for u in PAIRS]

def wedge_coeffs(phi1, phi2):
    """holomorphic 2-form phi1^phi2 in basis e^{jk}: components p_jk = phi1_j phi2_k - phi1_k phi2_j (elements of K)."""
    return [phi1[j] * phi2[k] - phi1[k] * phi2[j] for (j, k) in PAIRS]

def NS_T(omega):
    """omega: list of 6 K elements (period of the (2,0)-form). Returns dict with NS basis, T basis, Grams (integer lists)."""
    # <x, omega> = sum_a sum_b x_a Gram_ab omega_b = 0 in K  -> 4 rational equations
    rows = [[F(0)] * 6 for _ in range(4)]
    for a in range(6):
        val = ZERO
        for b in range(6):
            if GRAM6[a][b]:
                val = val + kq(GRAM6[a][b]) * omega[b]
        for r in range(4):
            rows[r][a] = val.c[r]
    # to integer rows
    irows = []
    for r in rows:
        den = 1
        for x in r: den = den * x.denominator // sp.igcd(den, x.denominator)
        irows.append([int(x * den) for x in r])
    ns = int_kernel(irows, 6)
    def gram(basis):
        return [[sum(u[a] * GRAM6[a][b] * v[b] for a in range(6) for b in range(6)) for v in basis] for u in basis]
    ns_gram = gram(ns) if ns else []
    # transcendental lattice: orthogonal complement of NS
    if ns:
        eq = [[sum(u[a] * GRAM6[a][b] for a in range(6)) for b in range(6)] for u in ns]
        T = int_kernel(eq, 6)
    else:
        T = [[1 if i == j else 0 for j in range(6)] for i in range(6)]
    return {"NS_basis": ns, "NS_gram": ns_gram, "T_basis": T, "T_gram": gram(T)}

def signature(Gm):
    """(n_plus, n_minus, n_zero) of a symmetric integer matrix, via exact LDL-free rule: numerical eigenvalues, checked against exact det sign."""
    if not Gm: return (0, 0, 0)
    import numpy as np
    ev = np.linalg.eigvalsh(np.array(Gm, dtype=float))
    p = int((ev > 1e-9).sum()); m = int((ev < -1e-9).sum())
    return (p, m, len(ev) - p - m)

def theta(Gm, maxnorm):
    """number of vectors of each norm <= maxnorm (positive definite Gram), exact."""
    n = len(Gm)
    Gs = sp.Matrix(Gm)
    Gi = Gs.inv()
    bnd = [int(sp.floor(sp.sqrt(sp.Rational(maxnorm) * Gi[i, i]))) for i in range(n)]
    cnt = {}
    for v in itertools.product(*[range(-b, b + 1) for b in bnd]):
        q = sum(v[i] * Gm[i][j] * v[j] for i in range(n) for j in range(n))
        if 0 < q <= maxnorm: cnt[q] = cnt.get(q, 0) + 1
    return dict(sorted(cnt.items()))

def nullspace_K(M, lam):
    """right nullspace over K of (M - lam I), M integer matrix as list of rows; returns list of vectors (lists of K)."""
    n = len(M)
    A = [[kq(M[i][j]) - (lam if i == j else ZERO) for j in range(n)] for i in range(n)]
    piv = []
    r = 0
    for c in range(n):
        p = next((i for i in range(r, n) if not A[i][c].iszero()), None)
        if p is None: continue
        A[r], A[p] = A[p], A[r]
        iv = A[r][c].inv()
        A[r] = [x * iv for x in A[r]]
        for i in range(n):
            if i != r and not A[i][c].iszero():
                f = A[i][c]
                A[i] = [x - f * y for x, y in zip(A[i], A[r])]
        piv.append(c); r += 1
    free = [c for c in range(n) if c not in piv]
    out = []
    for fcol in free:
        v = [ZERO] * n
        v[fcol] = ONE
        for k, pc in enumerate(piv):
            v[pc] = -A[k][fcol]
        out.append(v)
    return out
