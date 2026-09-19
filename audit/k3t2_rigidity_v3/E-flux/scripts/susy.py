"""Exact checks of TT's supersymmetry conditions for a given flux (alpha_x, beta_x, alpha_y, beta_y) and
given (phi, tau).  Vectors live in the lattice with Gram matrix H (integer coordinates); the complex
inner product is the C-BILINEAR extension of the lattice form, as in TT (bar = complex conjugate).
"""
import sympy as sp
from sympy import I

def cdot(u, v, H):
    return sp.expand(sum(u[i] * H[i, j] * v[j] for i in range(len(u)) for j in range(len(v)) if H[i, j] != 0))

def conj(u):
    return [sp.conjugate(x) for x in u]

def iszero(x):
    x = sp.nsimplify(sp.simplify(x)) if False else sp.simplify(x)
    return x == 0

def tt_conditions(H, ax, bx, ay, by, phi, tau):
    """(3.5) nx = ax - phi bx, ny = ay - phi by ;  (3.16a-c);  (3.15)/(3.31)."""
    nx = [a - phi * b for a, b in zip(ax, bx)]
    ny = [a - phi * b for a, b in zip(ay, by)]
    Gzb = [nx[i] * tau - ny[i] for i in range(len(nx))]           # Gz-bar  (up to overall sign)  (3.8)
    Gz = [nx[i] * sp.conjugate(tau) - ny[i] for i in range(len(nx))]  # (3.7)
    c316a = cdot([sp.conjugate(nx[i]) * tau - sp.conjugate(ny[i]) for i in range(len(nx))], Gzb, H)
    c316b = cdot(Gz, Gzb, H)
    c316c = cdot(Gzb, Gzb, H)
    c331 = cdot(Gzb, conj(Gzb), H)   # (3.15)/(3.31): Gzb . conj(Gzb) > 0
    return {"3.16a": sp.simplify(c316a), "3.16b": sp.simplify(c316b), "3.16c": sp.simplify(c316c),
            "3.15_3.31_value": sp.simplify(c331), "Gzbar": [sp.simplify(x) for x in Gzb]}

def gram_of(H, vecs):
    return sp.Matrix(len(vecs), len(vecs), lambda i, j: cdot(vecs[i], vecs[j], H))

def sig_rank(Gm):
    """exact (n+, n-, n0) of a symmetric rational matrix by congruence diagonalisation over Q."""
    A = sp.Matrix(Gm).copy()
    n = A.shape[0]
    pos = neg = 0
    idx = list(range(n))
    A = sp.Matrix(A)
    while A.shape[0] > 0:
        m = A.shape[0]
        # find nonzero diagonal
        k = next((i for i in range(m) if A[i, i] != 0), None)
        if k is None:
            # find off-diagonal nonzero, add rows/cols to create a nonzero diagonal
            pair = next(((i, j) for i in range(m) for j in range(i + 1, m) if A[i, j] != 0), None)
            if pair is None:
                break
            i, j = pair
            E = sp.eye(m); E[i, j] = 1
            A = E * A * E.T
            continue
        d = A[k, k]
        pos += int(bool(d > 0)); neg += int(bool(d < 0))
        # eliminate
        keep = [i for i in range(m) if i != k]
        B = sp.Matrix(len(keep), len(keep), lambda a, b: A[keep[a], keep[b]] - A[keep[a], k] * A[k, keep[b]] / d)
        A = B
    zeros = n - pos - neg
    return int(pos), int(neg), int(zeros)

def int_kernel_basis(M):
    """Z-basis of {x in Z^n : M x = 0} for an integer matrix M (m x n), via column-style HNF with transform."""
    M = sp.Matrix(M)
    m, n = M.shape
    A = M.T.copy()            # n x m ; do row operations on [A | I_n]
    Uu = sp.eye(n)
    row = 0
    for col in range(m):
        # bring gcd of column col (rows >= row) to position row
        while True:
            nz = [i for i in range(row, n) if A[i, col] != 0]
            if len(nz) <= 1:
                break
            i0 = min(nz, key=lambda i: abs(A[i, col]))
            for i in nz:
                if i != i0:
                    q = A[i, col] // A[i0, col]
                    A[i, :] = A[i, :] - q * A[i0, :]
                    Uu[i, :] = Uu[i, :] - q * Uu[i0, :]
        nz = [i for i in range(row, n) if A[i, col] != 0]
        if nz:
            i0 = nz[0]
            A.row_swap(row, i0); Uu.row_swap(row, i0)
            row += 1
    # rows of Uu from index `row` on give kernel of M (rows of A that vanished)
    return [list(Uu[i, :]) for i in range(row, n)]

def saturation_basis(vecs, n):
    """Z-basis of (Q-span of vecs) ∩ Z^n."""
    F = sp.Matrix(vecs)                       # k x n
    ns = F.nullspace()                        # rational basis of orthogonal-complement-for-standard-dot
    if not ns:
        return [list(map(int, r)) for r in sp.eye(n).tolist()]
    Cm = sp.Matrix.hstack(*ns).T              # c x n  rational
    den = sp.ilcm(*[x.q for x in Cm]) if Cm.shape[0] else 1
    Ci = (Cm * den).applyfunc(int)
    return [[int(x) for x in r] for r in int_kernel_basis(Ci)]

def orbifold_hits(H, flux, Omega_c, n):
    """Is there a nonzero lattice vector v in V_flux ∩ Gamma with v·Omega = 0 (complex bilinear)?
    Omega_c: complex vector (sympy) proportional to Omega.  Works when the entries of Re/Im(Omega·b) lie in Q(sqrt(r)).
    Returns (bool, basis_used).  Exact."""
    B = saturation_basis([list(map(int, f)) for f in flux if any(f)], n)
    s = sp.symbols("s", positive=True)
    rows = []
    vals = [sp.expand(cdot([int(x) for x in b], Omega_c, H)) for b in B]
    eqs_re = [sp.re(v) for v in vals]; eqs_im = [sp.im(v) for v in vals]
    allv = eqs_re + eqs_im
    # decompose each real number a + b*sqrt(r): collect sqrt atoms
    sqrts = set()
    for v in allv:
        sqrts |= {a for a in sp.sympify(v).atoms(sp.Pow) if a.exp == sp.Rational(1, 2)}
    sqrts = sorted(sqrts, key=str)
    assert len(sqrts) <= 1, sqrts
    cols = []
    for v in allv:
        v = sp.expand(sp.simplify(v))
        if sqrts:
            q = sqrts[0]
            b0 = sp.simplify(v.subs(q, 0)); b1 = sp.simplify(sp.diff(v, q)) if False else sp.simplify((v - b0) / q)
            cols.append((b0, b1))
        else:
            cols.append((v, 0))
    k = len(B)
    rowsM = []
    for t in range(2):
        for part in range(2):
            rowsM.append([cols[part * k + j][t] for j in range(k)])
    Mm = sp.Matrix(rowsM)
    ns = Mm.nullspace()
    return (len(ns) > 0), B, len(ns)
