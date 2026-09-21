#!/usr/bin/env python3
"""
Genus of X0(N) and of the Fricke quotient X0(N)+ = X0(N)/w_N.

WHY THIS EXISTS. Dolgachev (1996) Thm 7.1 -- literature, Tier L, quoted through
Stream 1 -- gives the coarse moduli space of M_n-polarized K3 surfaces as the
Fricke modular curve H/Gamma_0(n)+. A Hauptmodul (a single uniformizing
coordinate on that moduli space) exists exactly when that curve has genus 0.

⚠️ WHICH GROUP "Gamma_0(N)+" MEANS -- a limitation of the validation below.
This file computes the quotient of X0(N) by the FRICKE involution w_N alone,
which is Dolgachev's convention and the one the K3 moduli statement needs. Some
authors write Gamma_0(N)+ for Gamma_0(N) extended by ALL Atkin-Lehner
involutions, a group of order 2^omega(N). The two agree only when N has ONE prime
factor. The Ogg's-primes check in the self-test uses only PRIMES, where
omega(p) = 1 and the conventions coincide -- so it CANNOT distinguish them and
does not validate the composite-N rows. At this repository's N = 12 the
Atkin-Lehner group W(12) = {1, w_3, w_4, w_12} has order 4, so the distinction is
live exactly where it matters here. Treat the composite-N values as the
Fricke-only quotient by construction, not by test.

WHAT IT DOES NOT DO. It does NOT select a level. Both levels in play in this
programme -- Stream 1's N = 7 and this repository's implicit N = 12 -- have
genus 0, so the criterion passes both. It is included as a NEGATIVE CONTROL and
a sanity check on the level bookkeeping, not as a selection principle. See
audit/K3_SELECTION.md.

VALIDATION (run with --self-test):
  * genus_X0 against textbook values, including g(X0(37)) = 2, g(X0(49)) = 1,
    g(X0(50)) = 2, g(X0(121)) = 6;
  * class_number against 18 known class numbers, including h(-163) = 1,
    h(-71) = 7, h(-47) = 5, h(-44) = 3 (the last one requires counting only
    PRIMITIVE forms -- (2,2,6) has content 2 and must be excluded);
  * g(X0(N)+) is a NON-NEGATIVE INTEGER for every N in 5..300 (Riemann-Hurwitz
    would fail otherwise);
  * the PRIMES p with g(X0(p)+) = 0 come out as exactly
    {2,3,5,7,11,13,17,19,23,29,31,41,47,59,71} -- Ogg's supersingular primes,
    the primes dividing the order of the Monster. This is an independent check
    the author did not put in by hand.
"""
from math import gcd
from sympy import factorint, totient, divisors, isprime

__all__ = ["genus_X0", "class_number", "fricke_fixed_points", "genus_X0_plus"]


def _chi_minus1(p: int) -> int:
    """(-1/p): 0 at p = 2, +1 if p = 1 mod 4, -1 if p = 3 mod 4."""
    return 0 if p == 2 else (1 if p % 4 == 1 else -1)


def _chi_minus3(p: int) -> int:
    """(-3/p): 0 at p = 3, +1 if p = 1 mod 3, -1 if p = 2 mod 3."""
    return 0 if p == 3 else (1 if p % 3 == 1 else -1)


def genus_X0(N: int) -> int:
    """g(X0(N)) = 1 + mu/12 - nu2/4 - nu3/3 - nu_inf/2  (Diamond-Shurman Thm 3.1.1)."""
    if N == 1:
        return 0
    primes = list(factorint(N))
    mu = N
    for p in primes:
        mu = mu * (p + 1) // p
    nu2 = 0 if N % 4 == 0 else 1
    if nu2:
        for p in primes:
            nu2 *= 1 + _chi_minus1(p)
    nu3 = 0 if N % 9 == 0 else 1
    if nu3:
        for p in primes:
            nu3 *= 1 + _chi_minus3(p)
    nu_inf = sum(totient(gcd(d, N // d)) for d in divisors(N))
    g = float(1 + mu / 12 - nu2 / 4 - nu3 / 3 - nu_inf / 2)
    assert abs(g - round(g)) < 1e-9, (N, g)
    return int(round(g))


def class_number(D: int) -> int:
    """h(D) for D < 0: the number of PRIMITIVE reduced binary quadratic forms
    (a,b,c) with b^2 - 4ac = D, -a < b <= a <= c, and b >= 0 when a == c.

    Primitivity (gcd(a,b,c) = 1) is part of the definition and is load-bearing:
    without it h(-44) comes out 4 instead of 3, because (2,2,6) has content 2.
    """
    assert D < 0 and D % 4 in (0, 1), D
    h, a = 0, 1
    while 3 * a * a <= -D:
        for b in range(-a + 1, a + 1):
            num = b * b - D
            if num % (4 * a):
                continue
            c = num // (4 * a)
            if c < a or (a == c and b < 0):
                continue
            if gcd(gcd(a, abs(b)), c) != 1:
                continue
            h += 1
        a += 1
    return h


def fricke_fixed_points(N: int) -> int:
    """Number of fixed points of the Fricke involution w_N on X0(N), for N > 4."""
    f = class_number(-4 * N)
    if N % 4 == 3:
        f += class_number(-N)
    return f


def genus_X0_plus(N: int) -> int:
    """g(X0(N)+), from Riemann-Hurwitz on the degree-2 cover X0(N) -> X0(N)+:
    2g - 2 = 2(2g+ - 2) + f, so g+ = (2g + 2 - f)/4."""
    if N <= 4:
        # X0(N) already has genus 0 for N <= 4, and a quotient of a genus-0 curve
        # is genus 0. The fixed-point formula below is only valid for N > 4.
        assert genus_X0(N) == 0, N
        return 0
    val = float((2 * genus_X0(N) + 2 - fricke_fixed_points(N)) / 4)
    assert abs(val - round(val)) < 1e-9 and val >= -1e-9, (N, val)
    return int(round(val))


def _self_test() -> bool:
    ok = True
    known_g = {1:0,2:0,3:0,4:0,5:0,6:0,7:0,8:0,9:0,10:0,11:1,12:0,13:0,14:1,15:1,
               16:0,17:1,18:0,19:1,20:1,21:1,22:2,23:2,24:1,25:0,37:2,49:1,50:2,121:6}
    bad = {N: (genus_X0(N), g) for N, g in known_g.items() if genus_X0(N) != g}
    print(f"genus_X0 vs textbook      : {'PASS' if not bad else f'FAIL {bad}'}"); ok &= not bad

    known_h = {-3:1,-4:1,-7:1,-8:1,-11:1,-15:2,-19:1,-20:2,-23:3,-24:2,-43:1,
               -44:3,-47:5,-67:1,-71:7,-163:1,-84:4,-120:4}
    bad = {D: (class_number(D), h) for D, h in known_h.items() if class_number(D) != h}
    print(f"class_number vs known     : {'PASS' if not bad else f'FAIL {bad}'}"); ok &= not bad

    try:
        gp = {N: genus_X0_plus(N) for N in range(1, 301)}
        print("g(X0(N)+) integral & >= 0 : PASS")
    except AssertionError as e:
        print(f"g(X0(N)+) integral & >= 0 : FAIL {e}"); return False

    ogg = {2,3,5,7,11,13,17,19,23,29,31,41,47,59,71}
    got = {p for p in range(2, 301) if isprime(p) and genus_X0_plus(p) == 0}
    print(f"genus-0 primes = Ogg's set: {'PASS' if got == ogg else f'FAIL {sorted(got ^ ogg)}'}")
    ok &= got == ogg
    return ok


if __name__ == "__main__":
    import sys
    if "--self-test" in sys.argv:
        sys.exit(0 if _self_test() else 1)
    print(f"{'N':>4}  {'g(X0(N))':>9}  {'fix(w_N)':>9}  {'g(X0(N)+)':>10}")
    for N in (7, 11, 12, 37, 49, 71):
        print(f"{N:>4}  {genus_X0(N):>9}  {fricke_fixed_points(N):>9}  {genus_X0_plus(N):>10}")
