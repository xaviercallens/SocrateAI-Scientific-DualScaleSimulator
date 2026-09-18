"""
Appell-Lerch sum mu(tau,z), handled WITHOUT ever expanding an infinite-tail
series: the only term that naively looks singular (n=0: 1/(1-y)) is combined
ANALYTICALLY, before any truncation, with the (1-y) factor that T1(tau,z)
carries identically at every q-order (T1 = y^{-1/2}(1-y) q^{1/8} PROD1(q,y)).
That exact cancellation is what makes every downstream object below a
finite, honest Laurent polynomial per q-order -- no regularization, no
infinite geometric-series tails anywhere.

Definition used (Zwegers-normalized, real form matching T1 = theta_1/i):
    S(tau,z) := sum_{n in Z} (-1)^n q^{n(n+1)/2} y^n / (1 - y q^n)
    mu(tau,z) := y^{1/2} S(tau,z) / T1(tau,z)

We never form S or mu directly (S's n=0 term alone is not a finite-per-order
Laurent series). Instead we form
    Psi(tau,z) := T1(tau,z) * S(tau,z)
term by term:
  - n = 0:      T1 * 1/(1-y) = y^{-1/2} q^{1/8} PROD1(q,y)   [exact, no division]
  - n >= 1:     T1 * (-1)^n q^{n(n+1)/2} * sum_k y^{n+k} q^{nk}
  - n = -m, m>=1: T1 * (-1)^m q^{m(m-1)/2} * (-1) * sum_k y^{-m-1-k} q^{m(k+1)}
each of which is finite per q-order once truncated at imax (verified by
truncation-stability in run_all.py).

Then:  mu * T1^2 = y^{1/2} * Psi(tau,z)     (finite; this is what "24*mu-term"
combines with in the decomposition once multiplied through by theta_1^2).
"""
from fractions import Fraction as Fr
from series2d import mul, add, shift
from thetas import prod1


def _shift_y(series, dj):
    return {(i, j + dj): v for (i, j), v in series.items()}


def s_nonzero_terms(N, imax):
    """sum_{n != 0} (-1)^n q^{n(n+1)/2} y^n/(1-y q^n), as a finite 2D series
    (Fraction coefficients), truncated at q-index imax. Returns the dict.
    """
    out = {}

    def add_term(i0, base_j, dj_per_k, sign, nstep_i, kmax):
        for k in range(kmax + 1):
            i = i0 + k * nstep_i
            if i > imax:
                break
            j = base_j + dj_per_k * k
            out[(i, j)] = out.get((i, j), Fr(0)) + sign

    # n >= 1: (-1)^n q^{n(n+1)/2} sum_k y^{n+k} q^{nk}   (j increases by 2 per k)
    n = 1
    while 4 * n * (n + 1) <= imax:  # i0 = 8*n(n+1)/2 = 4n(n+1)
        i0 = 4 * n * (n + 1)
        base_j = 2 * n
        sign = Fr(-1) ** n
        kmax = (imax - i0) // (8 * n)
        add_term(i0, base_j, 2, sign, 8 * n, kmax)
        n += 1

    # n = -m, m>=1: (-1)^m q^{m(m-1)/2} * (-1) * sum_k y^{-m-1-k} q^{m(k+1)}
    # (j DECREASES by 2 per k)
    m = 1
    while 4 * m * (m + 1) <= imax:  # i0 = 8*[m(m-1)/2+m] = 8*m(m+1)/2 = 4m(m+1)
        i0 = 4 * m * (m + 1)
        base_j = -2 * (m + 1)
        sign = -(Fr(-1) ** m)
        kmax = (imax - i0) // (8 * m)
        if kmax >= 0:
            add_term(i0, base_j, -2, sign, 8 * m, kmax)
        m += 1

    return {k: v for k, v in out.items() if v != 0}


def psi(N, imax, T1):
    """Psi(tau,z) = T1(tau,z) * S(tau,z), computed finitely (see module doc).
    N: product truncation length for PROD1 (>= imax/8 + 2 recommended).
    """
    # n = 0 piece: T1 * 1/(1-y) = y^{-1/2} q^{1/8} PROD1(q,y)
    p1 = prod1(N, imax)
    n0_piece = shift(_shift_y(p1, -1), di=1)  # * q^{1/8} (di=1) * y^{-1/2} (dj=-1)

    s_rest = s_nonzero_terms(N, imax)
    rest_piece = mul(T1, s_rest, imax)

    return add(n0_piece, rest_piece)
