"""
Shared exact-arithmetic helpers for K3xT2 rigidity track (Track C).
No dependency on anything outside this file / the standard library / sympy.
All numbers are exact (fractions.Fraction or sympy Rational/Integer).
"""
from fractions import Fraction
import sympy as sp


def congruence_diagonalize(M):
    """
    Exact symmetric congruence diagonalization over Q.
    Input: M = list of lists of Fraction (symmetric, size n).
    Returns: list of diagonal entries D (length n, Fractions) such that
    there exists invertible P (over Q) with P^T M P = diag(D).
    By Sylvester's law of inertia, (#pos(D), #neg(D), #zero(D)) is an
    invariant of M -- the signature -- independent of the elimination path.
    Implementation: classic symmetric Gaussian elimination for congruence,
    with a hyperbolic-pair trick when the working diagonal is all zero but
    some off-diagonal entry is nonzero (needed for e.g. U = [[0,1],[1,0]]).
    """
    n = len(M)
    A = [[Fraction(x) for x in row] for row in M]
    diag = []
    active = list(range(n))  # indices of rows/cols still "live"

    while active:
        # find a nonzero diagonal entry among active indices
        piv = None
        for i in active:
            if A[i][i] != 0:
                piv = i
                break
        if piv is None:
            # all diagonal entries among active are zero; find nonzero off-diag
            found = None
            for idx_a, i in enumerate(active):
                for j in active[idx_a + 1:]:
                    if A[i][j] != 0:
                        found = (i, j)
                        break
                if found:
                    break
            if found is None:
                # remaining block is entirely zero
                for i in active:
                    diag.append(Fraction(0))
                break
            i, j = found
            # basis change: e_i -> e_i + e_j  (congruence: row_i += row_j, col_i += col_j)
            for k in range(n):
                A[i][k] = A[i][k] + A[j][k]
            for k in range(n):
                A[k][i] = A[k][i] + A[k][j]
            piv = i

        p = A[piv][piv]
        diag.append(p)
        # eliminate all other active rows/cols against pivot `piv`
        others = [k for k in active if k != piv]
        for k in others:
            factor = A[k][piv] / p
            if factor != 0:
                for m in range(n):
                    A[k][m] = A[k][m] - factor * A[piv][m]
                for m in range(n):
                    A[m][k] = A[m][k] - factor * A[m][piv]
        active = others

    return diag


def signature_from_diag(diag):
    pos = sum(1 for d in diag if d > 0)
    neg = sum(1 for d in diag if d < 0)
    zero = sum(1 for d in diag if d == 0)
    return pos, neg, zero


def signature_congruence(M):
    diag = congruence_diagonalize(M)
    return signature_from_diag(diag)


def signature_sturm(M):
    """
    Independent, second method: exact eigenvalue-sign counting via the
    characteristic polynomial's exact real roots WITH MULTIPLICITY
    (sympy Poly.real_roots(), which isolates real roots of each square-free
    factor via Sturm sequences and repeats each root by its multiplicity --
    unlike Poly.count_roots(), which counts DISTINCT roots only and silently
    undercounts when the matrix has repeated eigenvalues, e.g. several copies
    of the hyperbolic plane U). No floats anywhere: comparisons of CRootOf /
    Rational roots to 0 are exact.
    """
    n = len(M)
    Msym = sp.Matrix(M)
    x = sp.symbols('x_char')
    charpoly = Msym.charpoly(x).as_expr()
    p = sp.Poly(charpoly, x, domain='QQ')
    roots = p.real_roots()  # exact, WITH multiplicity (list length == total real mult.)
    pos = sum(1 for r in roots if r > 0)
    neg = sum(1 for r in roots if r < 0)
    zero = sum(1 for r in roots if r == 0)
    # any non-real roots would indicate a non-symmetric-matrix bug; symmetric
    # real matrices have all-real spectra, so this should always be n.
    assert pos + neg + zero == n, (
        f"real_roots() found {len(roots)} real roots (with multiplicity) "
        f"but matrix has size {n}; matrix may not be symmetric."
    )
    return int(pos), int(neg), int(zero)


def det_fraction_matrix(M):
    """Exact determinant of an integer/Fraction matrix via sympy Rational."""
    Msp = sp.Matrix([[sp.Rational(x) for x in row] for row in M])
    return Msp.det()


def leading_principal_minors(M):
    n = len(M)
    Msp = sp.Matrix([[sp.Integer(int(x)) for x in row] for row in M])
    return [Msp[:k, :k].det() for k in range(1, n + 1)]


def block_diag(*blocks):
    """Assemble a block-diagonal integer Gram matrix from square blocks (lists of lists)."""
    n = sum(len(b) for b in blocks)
    M = [[0] * n for _ in range(n)]
    offset = 0
    for b in blocks:
        k = len(b)
        for i in range(k):
            for j in range(k):
                M[offset + i][offset + j] = b[i][j]
        offset += k
    return M


U = [[0, 1], [1, 0]]  # hyperbolic plane Gram matrix
