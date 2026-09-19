"""
Independent re-implementation of the DMVV product

    Sum_k G_k p^k = Prod_{r>=1, s>=0, t in Z} (1 - p^r q^s y^t)^{-c(4rs-t^2)}

(DMZ (5.13)-(5.15) / DualScaleDyons/DMVV.lean's `dmvv`), truncated at p^K, q^Q.

This mirrors the ALGORITHM of Lean's `mulFactor`/`dmvv` (the two are the same
defining product, so the bookkeeping is necessarily the same shape), but is
written fresh in Python with fractions.Fraction, reading c(D) from our OWN
theta-function computation (dyons/theta_forms_cache.json), not from Lean.
"""
import json
from fractions import Fraction as Fr
from nlseries import clean, neg_binom_coeff

CACHE_PATH = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/dyons/theta_forms_cache.json"


def _pf(s):
    if s is None:
        return None
    if "/" in s:
        a, b = s.split("/")
        return Fr(int(a), int(b))
    return Fr(int(s))


def load_cD(cache_path=CACHE_PATH):
    with open(cache_path) as f:
        cache = json.load(f)
    cD = {int(k): _pf(v) for k, v in cache["cD_table"].items()}
    dmax = max(cD.keys())
    return cD, dmax


def cK3(D, cD, dmax):
    """c(D): 0 for D < -1; else looked up from our independently-computed
    theta table; asserts D is within the table's reach (raises rather than
    silently returning 0 for an out-of-range D, matching the Lean file's own
    `dmvv_reachable` / `dmvv_unreachable_example` discipline)."""
    if D < -1:
        return Fr(0)
    if D not in cD:
        if D > dmax:
            raise ValueError(f"c({D}) not available: beyond theta cache max D={dmax}")
        return Fr(0)  # D of the wrong residue class (D%4 in {1,2}): genuinely zero
    return cD[D]


def mul_factor(G, K, Q, r, s, t, e):
    """Multiply Ser3 G (list of K+1 (n,l)-dicts) by (1 - p^r q^s y^t)^{-e}."""
    newG = [{} for R in range(K + 1)]
    for R in range(K + 1):
        acc = {}
        jmax = R // r
        # j=0 (identity) is included, contributing G[R] itself, since coeff(e,0)=1
        for j in range(0, jmax + 1):
            b = neg_binom_coeff(e, j)
            if b == 0:
                continue
            src = G[R - r * j]
            for (n, l), v in src.items():
                nn = n + s * j
                if nn > Q:
                    continue
                ll = l + t * j
                acc[(nn, ll)] = acc.get((nn, ll), Fr(0)) + b * v
        newG[R] = clean(acc)
    return newG


def dmvv(K, Q, cD, dmax):
    """Sum_{k<=K} G_k p^k, each G_k truncated at q^Q, exact Fraction coeffs."""
    G = [{(0, 0): Fr(1)} if R == 0 else {} for R in range(K + 1)]
    for r in range(1, K + 1):
        for s in range(0, Q + 1):
            tmax = int((4 * r * s + 1) ** 0.5) + 2
            for t in range(-tmax, tmax + 1):
                D = 4 * r * s - t * t
                if D < -1:
                    continue
                e = cK3(D, cD, dmax)
                if e == 0:
                    continue
                G = mul_factor(G, K, Q, r, s, t, int(e))
    return G


def p24_direct(nmax):
    """p24(n): coefficient of q^n in Prod_k (1-q^k)^{-24}, n=0..nmax, computed
    by direct exact convolution (24-fold), independent implementation."""
    # start with delta_0
    s = [Fr(0)] * (nmax + 1)
    s[0] = Fr(1)
    for k in range(1, nmax + 1):
        # multiply s by (1-q^k)^{-24} = sum_j C(24+j-1,j) q^{kj}
        factor = {}
        j = 0
        while k * j <= nmax:
            factor[k * j] = Fr(comb_pos(24, j))
            j += 1
        news = [Fr(0)] * (nmax + 1)
        for n1, v1 in enumerate(s):
            if v1 == 0:
                continue
            for e, v2 in factor.items():
                n = n1 + e
                if n > nmax:
                    continue
                news[n] += v1 * v2
        s = news
    return s


def comb_pos(e, j):
    """C(e+j-1, j) for e>0 (i.e. coefficient of x^j in (1-x)^{-e})."""
    from math import comb
    return comb(e + j - 1, j)
