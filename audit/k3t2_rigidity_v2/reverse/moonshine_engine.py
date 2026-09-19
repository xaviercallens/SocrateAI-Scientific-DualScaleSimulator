"""
Independent Python re-implementation of LeanMaster's DualScaleMoonshine engine
(QSeries.lean, Twining.lean, TwiningAll.lean, Characters.lean, CharactersAll.lean,
Decompositions.lean), built fresh from the SAME defining formulas (Cheng-Duncan-Harvey
arXiv:1204.2779, "1204_2779.txt" in this repo's papers/foundations/) rather than
transliterated from the Lean source. Every constant (character table entries, F_g
combinations of Lambda_N/eta-quotients, etc.) is a Tier-L literature input EXACTLY as it
is in Lean (this reverse pass does not re-derive the moonshine module from scratch -- that
is Gannon's theorem, out of scope); what is independent here is the ARITHMETIC ENGINE
(power series truncation, eta-quotient construction, the ring Z[b_n] inner-product
machinery) and the fact that it is run to N=20 instead of Lean's N=9.

Truncation order N is a free parameter throughout (Lean hard-codes 9 in every top-level
theorem; nothing here is capped at 9 except where explicitly noted).
"""
from fractions import Fraction as Fr

# ---------------------------------------------------------------------------
# Truncated integer q-series as lists of ints, index 0..N (matches QSeries.lean's
# convention: index n <-> coefficient of q^n, working with 24*D*H_g etc. to stay in Z).
# ---------------------------------------------------------------------------

def mul_trunc(N, a, b):
    out = [0] * (N + 1)
    for n in range(N + 1):
        s = 0
        for k in range(n + 1):
            s += (a[k] if k < len(a) else 0) * (b[n - k] if n - k < len(b) else 0)
        out[n] = s
    return out

def div_trunc(N, r, p):
    """r/p truncated, p[0] must be 1 (or -1)."""
    assert p[0] in (1, -1)
    inv0 = p[0]
    acc = []
    for n in range(N + 1):
        s = r[n] if n < len(r) else 0
        for k in range(n):
            s -= acc[k] * (p[n - k] if n - k < len(p) else 0)
        acc.append(s * inv0)
    return acc

def one_minus_qn(N, n):
    """`1 - q^n` truncated at q^N (n=0 gives the constant series `1`, matching Lean's
    `oneMinusQ`: the `k=0` branch takes priority over `k=n`)."""
    out = [0] * (N + 1)
    out[0] = 1
    if n != 0 and n <= N:
        out[n] = -1
    return out

def eta3(N):
    acc = one_minus_qn(N, 0)
    for i in range(N):
        f = one_minus_qn(N, i + 1)
        acc = mul_trunc(N, mul_trunc(N, mul_trunc(N, acc, f), f), f)
    return acc

def sigma1(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)

def f2_coeff(n):
    total = 0
    for s in range(1, 2 * n + 1):
        if (2 * n) % s != 0:
            continue
        r = 2 * n // s
        if not (s < r and (r - s) % 2 == 1):
            continue
        total += (-s if r % 2 == 1 else s)
    return total

def numer(N):
    out = [0] * (N + 1)
    out[0] = -2
    for n in range(1, N + 1):
        out[n] = 48 * (sigma1(n) + f2_coeff(n))
    return out

def hComputed(N):
    """q^{1/8} H^{(2)}, i.e. Cheng-Harrison's H, computed from -2E2+48F2 over eta^3."""
    return div_trunc(N, numer(N), eta3(N))

def add_s(a, b):
    return [x + y for x, y in zip(a, b)]

def scal_s(c, a):
    return [c * x for x in a]

def lambda24(N, Nlev, c):
    out = [0] * (N + 1)
    out[0] = c * Nlev * (Nlev - 1)
    for n in range(1, N + 1):
        term = sigma1(n) - (Nlev * sigma1(n // Nlev) if n % Nlev == 0 else 0)
        out[n] = c * 24 * Nlev * term
    return out

def twined24(N, chi, Nlev, c):
    lam = lambda24(N, Nlev, c)
    num = [chi * x for x in numer(N)]
    return div_trunc(N, add_s(num, lam), eta3(N))

# ---------------------------------------------------------------------------
# Eta quotients (TwiningAll.lean's etaQ/etaProd/etaPow), general N.
# ---------------------------------------------------------------------------

def eta_pow(N, k, e):
    acc = one_minus_qn(N, 0)
    i = 1
    while k * i <= N:
        f = one_minus_qn(N, k * i)
        for _ in range(e):
            acc = mul_trunc(N, acc, f)
        i += 1
    return acc

def eta_prod(N, fs):
    acc = one_minus_qn(N, 0)
    for k, e in fs:
        acc = mul_trunc(N, acc, eta_pow(N, k, e))
    return acc

def eta_weight(fs):
    return sum(k * e for k, e in fs)

def eta_q(N, num, den):
    s = (eta_weight(num) - eta_weight(den)) // 24
    numer_s = eta_prod(N, num)
    den_s = eta_prod(N, den)
    q = div_trunc(N, numer_s, den_s)
    return ([0] * s + q)[: N + 1]

def f11(N): return eta_q(N, [(1, 2), (11, 2)], [])
def f14(N): return eta_q(N, [(1, 1), (2, 1), (7, 1), (14, 1)], [])
def f15(N): return eta_q(N, [(1, 1), (3, 1), (5, 1), (15, 1)], [])
def f23a(N):
    t1 = eta_q(N, [(1, 3), (23, 3)], [(2, 1), (46, 1)])
    t2 = scal_s(3, eta_q(N, [(1, 2), (23, 2)], []))
    t3 = scal_s(4, eta_q(N, [(1, 1), (2, 1), (23, 1), (46, 1)], []))
    t4 = scal_s(4, eta_q(N, [(2, 2), (46, 2)], []))
    return add_s(add_s(t1, t2), add_s(t3, t4))
def f23b(N): return eta_q(N, [(1, 2), (23, 2)], [])

def data_all(N):
    """(D, chi_g, 24*D*F_g) for the 16 classes of TwiningAll.lean, at general N."""
    F11, F14, F15 = f11(N), f14(N), f15(N)
    F23a, F23b = f23a(N), f23b(N)
    return {
        "2B": (1, 0, add_s(lambda24(N, 2, 24), lambda24(N, 4, -8))),
        "3B": (1, 0, scal_s(-48, eta_q(N, [(1, 6)], [(3, 2)]))),
        "4A": (1, 0, add_s(add_s(lambda24(N, 2, -4), lambda24(N, 4, 6)), lambda24(N, 8, -2))),
        "4B": (1, 4, add_s(lambda24(N, 2, 4), lambda24(N, 4, -4))),
        "4C": (1, 0, scal_s(-48, eta_q(N, [(1, 4), (2, 2)], [(4, 2)]))),
        "6A": (1, 2, add_s(add_s(lambda24(N, 2, 2), lambda24(N, 3, 2)), lambda24(N, 6, -2))),
        "6B": (1, 0, scal_s(-48, eta_q(N, [(1, 2), (2, 2), (3, 2)], [(6, 2)]))),
        "8A": (1, 2, add_s(lambda24(N, 4, 1), lambda24(N, 8, -1))),
        "10A": (1, 0, scal_s(-48, eta_q(N, [(1, 3), (2, 1), (5, 1)], [(10, 1)]))),
        "11A": (5, 2, add_s(lambda24(N, 11, -2), scal_s(528, F11))),
        "12A": (1, 0, scal_s(-48, eta_q(N, [(1, 3), (4, 2), (6, 3)], [(2, 1), (3, 1), (12, 2)]))),
        "12B": (1, 0, scal_s(-48, eta_q(N, [(1, 4), (4, 1), (6, 1)], [(2, 1), (12, 1)]))),
        "14AB": (3, 1, add_s(add_s(add_s(lambda24(N, 2, 1), lambda24(N, 7, 1)), lambda24(N, 14, -1)),
                              scal_s(336, F14))),
        "15AB": (4, 1, add_s(add_s(add_s(lambda24(N, 3, 1), lambda24(N, 5, 1)), lambda24(N, 15, -1)),
                              scal_s(360, F15))),
        "21AB": (3, 0, add_s(scal_s(-168, eta_q(N, [(1, 3), (7, 3)], [(3, 1), (21, 1)])),
                              scal_s(24, eta_q(N, [(1, 6)], [(3, 2)])))),
        "23AB": (11, 1, add_s(add_s(lambda24(N, 23, -1), scal_s(552, F23a)), scal_s(1656, F23b))),
    }

def twined_d(N, d):
    D, chi, Fg24D = d
    num = scal_s(D * chi, numer(N))
    return div_trunc(N, add_s(num, Fg24D), eta3(N))

# ---------------------------------------------------------------------------
# Character table (CharactersAll.lean), colK/classSize, inner4 -- general N via
# computedSeries(N).
# ---------------------------------------------------------------------------

CLASS_NAMES = ["1A", "2A", "2B", "3A", "3B", "4A", "4B", "4C", "5A", "6A", "6B", "7A", "7B",
               "8A", "10A", "11A", "12A", "12B", "14A", "14B", "15A", "15B", "21A", "21B", "23A", "23B"]

COL_K = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 0, 0, 0, 0, 2, 2, 4, 4, 2, 2, 6, 6]

CENT = [244823040, 21504, 7680, 1080, 504, 384, 128, 96, 60, 24, 24, 42, 42, 16, 20, 11, 12, 12,
        14, 14, 15, 15, 21, 21, 23, 23]
ORDER_M24 = 244823040
CLASS_SIZE = [ORDER_M24 // c for c in CENT]

# CDH Table 8, this repo's M24RepDim order (transcribed independently from the same source
# text CharactersAll.lean cites, papers/foundations/1204_2779.txt Table 8; entries are
# (a,b) meaning a+b*omega, omega=b_n with n given by COL_K).
CHAR_TAB = [
 [(1,0)]*26,
 [(23,0),(7,0),(-1,0),(5,0),(-1,0),(-1,0),(3,0),(-1,0),(3,0),(1,0),(-1,0),(2,0),(2,0),(1,0),(-1,0),(1,0),(-1,0),(-1,0),(0,0),(0,0),(0,0),(0,0),(-1,0),(-1,0),(0,0),(0,0)],
 [(45,0),(-3,0),(5,0),(0,0),(3,0),(-3,0),(1,0),(1,0),(0,0),(0,0),(-1,0),(0,1),(-1,-1),(-1,0),(0,0),(1,0),(0,0),(1,0),(0,-1),(1,1),(0,0),(0,0),(0,1),(-1,-1),(-1,0),(-1,0)],
 [(45,0),(-3,0),(5,0),(0,0),(3,0),(-3,0),(1,0),(1,0),(0,0),(0,0),(-1,0),(-1,-1),(0,1),(-1,0),(0,0),(1,0),(0,0),(1,0),(1,1),(0,-1),(0,0),(0,0),(-1,-1),(0,1),(-1,0),(-1,0)],
 [(231,0),(7,0),(-9,0),(-3,0),(0,0),(-1,0),(-1,0),(3,0),(1,0),(1,0),(0,0),(0,0),(0,0),(-1,0),(1,0),(0,0),(-1,0),(0,0),(0,0),(0,0),(0,1),(-1,-1),(0,0),(0,0),(1,0),(1,0)],
 [(231,0),(7,0),(-9,0),(-3,0),(0,0),(-1,0),(-1,0),(3,0),(1,0),(1,0),(0,0),(0,0),(0,0),(-1,0),(1,0),(0,0),(-1,0),(0,0),(0,0),(0,0),(-1,-1),(0,1),(0,0),(0,0),(1,0),(1,0)],
 [(252,0),(28,0),(12,0),(9,0),(0,0),(4,0),(4,0),(0,0),(2,0),(1,0),(0,0),(0,0),(0,0),(0,0),(2,0),(-1,0),(1,0),(0,0),(0,0),(0,0),(-1,0),(-1,0),(0,0),(0,0),(-1,0),(-1,0)],
 [(253,0),(13,0),(-11,0),(10,0),(1,0),(-3,0),(1,0),(1,0),(3,0),(-2,0),(1,0),(1,0),(1,0),(-1,0),(-1,0),(0,0),(0,0),(1,0),(-1,0),(-1,0),(0,0),(0,0),(1,0),(1,0),(0,0),(0,0)],
 [(483,0),(35,0),(3,0),(6,0),(0,0),(3,0),(3,0),(3,0),(-2,0),(2,0),(0,0),(0,0),(0,0),(-1,0),(-2,0),(-1,0),(0,0),(0,0),(0,0),(0,0),(1,0),(1,0),(0,0),(0,0),(0,0),(0,0)],
 [(770,0),(-14,0),(10,0),(5,0),(-7,0),(2,0),(-2,0),(-2,0),(0,0),(1,0),(1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(-1,0),(1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,1),(-1,-1)],
 [(770,0),(-14,0),(10,0),(5,0),(-7,0),(2,0),(-2,0),(-2,0),(0,0),(1,0),(1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(-1,0),(1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(-1,-1),(0,1)],
 [(990,0),(-18,0),(-10,0),(0,0),(3,0),(6,0),(2,0),(-2,0),(0,0),(0,0),(-1,0),(0,1),(-1,-1),(0,0),(0,0),(0,0),(0,0),(1,0),(0,1),(-1,-1),(0,0),(0,0),(0,1),(-1,-1),(1,0),(1,0)],
 [(990,0),(-18,0),(-10,0),(0,0),(3,0),(6,0),(2,0),(-2,0),(0,0),(0,0),(-1,0),(-1,-1),(0,1),(0,0),(0,0),(0,0),(0,0),(1,0),(-1,-1),(0,1),(0,0),(0,0),(-1,-1),(0,1),(1,0),(1,0)],
 [(1035,0),(27,0),(35,0),(0,0),(6,0),(3,0),(-1,0),(3,0),(0,0),(0,0),(2,0),(-1,0),(-1,0),(1,0),(0,0),(1,0),(0,0),(0,0),(-1,0),(-1,0),(0,0),(0,0),(-1,0),(-1,0),(0,0),(0,0)],
 [(1035,0),(-21,0),(-5,0),(0,0),(-3,0),(3,0),(3,0),(-1,0),(0,0),(0,0),(1,0),(0,2),(-2,-2),(-1,0),(0,0),(1,0),(0,0),(-1,0),(0,0),(0,0),(0,0),(0,0),(0,-1),(1,1),(0,0),(0,0)],
 [(1035,0),(-21,0),(-5,0),(0,0),(-3,0),(3,0),(3,0),(-1,0),(0,0),(0,0),(1,0),(-2,-2),(0,2),(-1,0),(0,0),(1,0),(0,0),(-1,0),(0,0),(0,0),(0,0),(0,0),(1,1),(0,-1),(0,0),(0,0)],
 [(1265,0),(49,0),(-15,0),(5,0),(8,0),(-7,0),(1,0),(-3,0),(0,0),(1,0),(0,0),(-2,0),(-2,0),(1,0),(0,0),(0,0),(-1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(1,0),(1,0),(0,0),(0,0)],
 [(1771,0),(-21,0),(11,0),(16,0),(7,0),(3,0),(-5,0),(-1,0),(1,0),(0,0),(-1,0),(0,0),(0,0),(-1,0),(1,0),(0,0),(0,0),(-1,0),(0,0),(0,0),(1,0),(1,0),(0,0),(0,0),(0,0),(0,0)],
 [(2024,0),(8,0),(24,0),(-1,0),(8,0),(8,0),(0,0),(0,0),(-1,0),(-1,0),(0,0),(1,0),(1,0),(0,0),(-1,0),(0,0),(-1,0),(0,0),(1,0),(1,0),(-1,0),(-1,0),(1,0),(1,0),(0,0),(0,0)],
 [(2277,0),(21,0),(-19,0),(0,0),(6,0),(-3,0),(1,0),(-3,0),(-3,0),(0,0),(2,0),(2,0),(2,0),(-1,0),(1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(-1,0),(-1,0),(0,0),(0,0)],
 [(3312,0),(48,0),(16,0),(0,0),(-6,0),(0,0),(0,0),(0,0),(-3,0),(0,0),(-2,0),(1,0),(1,0),(0,0),(1,0),(1,0),(0,0),(0,0),(-1,0),(-1,0),(0,0),(0,0),(1,0),(1,0),(0,0),(0,0)],
 [(3520,0),(64,0),(0,0),(10,0),(-8,0),(0,0),(0,0),(0,0),(0,0),(-2,0),(0,0),(-1,0),(-1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(1,0),(1,0),(0,0),(0,0),(-1,0),(-1,0),(1,0),(1,0)],
 [(5313,0),(49,0),(9,0),(-15,0),(0,0),(1,0),(-3,0),(-3,0),(3,0),(1,0),(0,0),(0,0),(0,0),(-1,0),(-1,0),(0,0),(1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0)],
 [(5796,0),(-28,0),(36,0),(-9,0),(0,0),(-4,0),(4,0),(0,0),(1,0),(-1,0),(0,0),(0,0),(0,0),(0,0),(1,0),(-1,0),(-1,0),(0,0),(0,0),(0,0),(1,0),(1,0),(0,0),(0,0),(0,0),(0,0)],
 [(5544,0),(-56,0),(24,0),(9,0),(0,0),(-8,0),(0,0),(0,0),(-1,0),(1,0),(0,0),(0,0),(0,0),(0,0),(-1,0),(0,0),(1,0),(0,0),(0,0),(0,0),(-1,0),(-1,0),(0,0),(0,0),(1,0),(1,0)],
 [(10395,0),(-21,0),(-45,0),(0,0),(0,0),(3,0),(-1,0),(3,0),(0,0),(0,0),(0,0),(0,0),(0,0),(1,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(0,0),(-1,0),(-1,0)],
]

def qmul(k, x, y):
    return (x[0] * y[0] - k * x[1] * y[1], x[0] * y[1] + x[1] * y[0] - x[1] * y[1])

def qconj(x):
    return (x[0] - x[1], -x[1])

def inner4(r, s):
    """Sum_g classSize(g) * <r(g), conj(s(g))> in Z[omega], returned as (rational, w7, w15, w23)."""
    acc = [0, 0, 0, 0]
    for j in range(26):
        k = COL_K[j]
        c = CLASS_SIZE[j]
        p = qmul(k, r[j], qconj(s[j]))
        acc[0] += c * p[0]
        if k == 2:
            acc[1] += c * p[1]
        elif k == 4:
            acc[2] += c * p[1]
        elif k == 6:
            acc[3] += c * p[1]
    return tuple(acc)

def computed_series(N):
    """The 26-column computed twined series q^{1/8} H_g, one list per class, general N."""
    d = data_all(N)
    s2b, s3b, s4a, s4b_, s4c = twined_d(N, d["2B"]), twined_d(N, d["3B"]), twined_d(N, d["4A"]), twined_d(N, d["4B"]), twined_d(N, d["4C"])
    s6a, s6b_ = twined_d(N, d["6A"]), twined_d(N, d["6B"])
    s8a, s10a = twined_d(N, d["8A"]), twined_d(N, d["10A"])
    s11a = twined_d(N, d["11A"])
    s12a, s12b_ = twined_d(N, d["12A"]), twined_d(N, d["12B"])
    s14 = twined_d(N, d["14AB"]); s15 = twined_d(N, d["15AB"])
    s21 = twined_d(N, d["21AB"]); s23 = twined_d(N, d["23AB"])

    def norm_D(D_series, D):
        return [x // (24 * D) for x in D_series]

    def norm1(series):
        return [x // 24 for x in series]

    s1a = twined24(N, 24, 1, 0)
    s2a = twined24(N, 8, 2, -16)
    s3a = twined24(N, 6, 3, -6)
    s5a = twined24(N, 4, 5, -2)
    s7 = twined24(N, 3, 7, -1)

    out = [norm1(s1a), norm1(s2a), norm_D(s2b, d["2B"][0]), norm1(s3a), norm_D(s3b, d["3B"][0]),
           norm_D(s4a, d["4A"][0]), norm_D(s4b_, d["4B"][0]), norm_D(s4c, d["4C"][0]), norm1(s5a),
           norm_D(s6a, d["6A"][0]), norm_D(s6b_, d["6B"][0]), norm1(s7), norm1(s7),
           norm_D(s8a, d["8A"][0]), norm_D(s10a, d["10A"][0]), norm_D(s11a, d["11A"][0]),
           norm_D(s12a, d["12A"][0]), norm_D(s12b_, d["12B"][0]),
           norm_D(s14, d["14AB"][0]), norm_D(s14, d["14AB"][0]),
           norm_D(s15, d["15AB"][0]), norm_D(s15, d["15AB"][0]),
           norm_D(s21, d["21AB"][0]), norm_D(s21, d["21AB"][0]),
           norm_D(s23, d["23AB"][0]), norm_D(s23, d["23AB"][0])]
    return out

def moonshine_inner(N, n, i, series=None):
    if series is None:
        series = computed_series(N)
    row = [(series[j][n] if n < len(series[j]) else 0, 0) for j in range(26)]
    return inner4(row, CHAR_TAB[i])
