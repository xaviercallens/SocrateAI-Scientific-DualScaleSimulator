"""
Jacobi theta functions theta_1..theta_4(tau,z) and eta(tau)^3, eta(tau)^6
as exact truncated Laurent series, on the series2d grid (i=8*qexp, j=2*yexp).

All product formulas are the standard triple-product expansions (Jacobi's
theta functions), e.g. Kac, "Infinite Dimensional Lie Algebras", or
Whittaker & Watson ch. 21. No numeric literature constants are used here --
only the defining product formulas. theta_1 = i * T1 with T1 real; we track
T1 and note the factor of i explicitly (theta_1^2 = -T1^2, real).
"""
from fractions import Fraction as Fr
from series2d import const, mul_many, add, scal, shift, y_to_1


def _one_minus_qn(n, imax):
    # (1 - q^n), i-index = 8n
    if 8 * n > imax:
        return const(1)
    return {(0, 0): Fr(1), (8 * n, 0): Fr(-1)}


def _pm_y_qn(sign_y, n, imax):
    # (1 + sign_y * y * q^n)  [+/- controls the y term's overall sign]
    if 8 * n > imax:
        return const(1)
    return {(0, 0): Fr(1), (8 * n, 2): Fr(sign_y)}


def _pm_yinv_qn(sign_y, n, imax):
    # (1 + sign_y * y^{-1} * q^n)
    if 8 * n > imax:
        return const(1)
    return {(0, 0): Fr(1), (8 * n, -2): Fr(sign_y)}


def _pm_y_qnhalf(sign_y, n, imax):
    # (1 + sign_y * y * q^{n-1/2}), i-index = 8n-4
    i = 8 * n - 4
    if i > imax:
        return const(1)
    return {(0, 0): Fr(1), (i, 2): Fr(sign_y)}


def _pm_yinv_qnhalf(sign_y, n, imax):
    i = 8 * n - 4
    if i > imax:
        return const(1)
    return {(0, 0): Fr(1), (i, -2): Fr(sign_y)}


def eta_pochhammer(N, imax):
    """prod_{n=1}^{N} (1 - q^n), plain series (i-index multiples of 8)."""
    facs = [_one_minus_qn(n, imax) for n in range(1, N + 1)]
    return mul_many(facs, imax)


def theta2(N, imax):
    """theta_2(tau,z) = q^{1/8}(y^{1/2}+y^{-1/2}) prod (1-q^n)(1+y q^n)(1+y^-1 q^n)."""
    prefac_q = const(1, i=1, j=0)
    prefac_y = {(0, 1): Fr(1), (0, -1): Fr(1)}
    facs = [prefac_q, prefac_y]
    for n in range(1, N + 1):
        facs.append(_one_minus_qn(n, imax))
        facs.append(_pm_y_qn(1, n, imax))
        facs.append(_pm_yinv_qn(1, n, imax))
    return mul_many(facs, imax)


def theta3(N, imax):
    facs = []
    for n in range(1, N + 1):
        facs.append(_one_minus_qn(n, imax))
        facs.append(_pm_y_qnhalf(1, n, imax))
        facs.append(_pm_yinv_qnhalf(1, n, imax))
    return mul_many(facs, imax)


def theta4(N, imax):
    facs = []
    for n in range(1, N + 1):
        facs.append(_one_minus_qn(n, imax))
        facs.append(_pm_y_qnhalf(-1, n, imax))
        facs.append(_pm_yinv_qnhalf(-1, n, imax))
    return mul_many(facs, imax)


def prod1(N, imax):
    """PROD1(q,y) = prod_{n=1}^N (1-q^n)(1-y q^n)(1-y^-1 q^n) (the product
    factor common to T1, with NO prefactor). PROD1(q,1) = prod(1-q^n)^3.
    """
    facs = []
    for n in range(1, N + 1):
        facs.append(_one_minus_qn(n, imax))
        facs.append(_pm_y_qn(-1, n, imax))
        facs.append(_pm_yinv_qn(-1, n, imax))
    return mul_many(facs, imax)


def theta1_over_i(N, imax):
    """T1 such that theta_1(tau,z) = i * T1(tau,z), T1 real:
    T1 = -(y^{1/2}-y^{-1/2}) q^{1/8} prod (1-q^n)(1-y q^n)(1-y^-1 q^n)
       =  y^{-1/2}(1-y) q^{1/8} PROD1(q,y).
    """
    prefac_q = const(1, i=1, j=0)
    prefac_y = {(0, 1): Fr(-1), (0, -1): Fr(1)}  # -(y^1/2 - y^-1/2)
    facs = [prefac_q, prefac_y, prod1(N, imax)]
    return mul_many(facs, imax)


def eta3(N, imax):
    """eta(tau)^3 = q^{1/8} prod(1-q^n)^3."""
    poch = eta_pochhammer(N, imax)
    from series2d import mul
    cube = mul(mul(poch, poch, imax), poch, imax)
    return shift(cube, di=1)


def eta6(N, imax):
    """eta(tau)^6 = q^{1/4} prod(1-q^n)^6."""
    poch = eta_pochhammer(N, imax)
    from series2d import mul
    e3 = mul(mul(poch, poch, imax), poch, imax)
    e6 = mul(e3, e3, imax)
    return shift(e6, di=2)
