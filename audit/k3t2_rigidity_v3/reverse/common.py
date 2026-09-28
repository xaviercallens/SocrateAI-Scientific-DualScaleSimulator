"""
Shared exact-arithmetic primitives for the reverse pass (audit/k3t2_rigidity_v3/reverse/).

These are independent re-implementations (own code, not copy-pasted Lean) of the
binary-quadratic-form machinery that DualScaleDyons/WhichK3.lean and
DualScaleDyons/Immortal.lean define, so that the reverse-pass scripts can scan far
beyond the ranges those Lean theorems check without importing anything from Lean.
All arithmetic is over Python's arbitrary-precision `int` -- exact, no floats.
"""
from __future__ import annotations
from fractions import Fraction


def reduced_forms(D: int) -> list[tuple[int, int, int]]:
    """Reduced positive-definite binary quadratic forms (a, b, c), disc = b^2 - 4ac = -D,
    |b| <= a <= c, b >= 0 whenever |b| = a or a = c. Counts ALL classes (primitive and not),
    matching WhichK3.reducedForms / nForms."""
    if D <= 0:
        return []
    forms = []
    a = 1
    while 3 * a * a <= D:
        for b in range(-a + 1, a + 1):
            num = b * b + D
            if num % (4 * a) != 0:
                continue
            c = num // (4 * a)
            if c < a:
                continue
            if b < 0 and c == a:
                continue
            forms.append((a, b, c))
        a += 1
    return forms


def n_forms(D: int) -> int:
    return len(reduced_forms(D))


def h12(D: int) -> int:
    """12 * H(D), the Hurwitz class number, by the weighted reduced-form count
    (weight 6 for a(x^2+y^2)-type forms, weight 4 for a(x^2+xy+y^2)-type forms,
    weight 12 otherwise). Matches Immortal.h12's definition, reimplemented independently."""
    if D == 0:
        return -1
    if D < 0 or D % 4 == 1 or D % 4 == 2:
        return 0
    total = 0
    for (a, b, c) in reduced_forms(D):
        if b == 0 and c == a:
            w = 6
        elif b == a and c == a:
            w = 4
        else:
            w = 12
        total += w
    return total


def is_k_square(k: int, D: int) -> bool:
    """D = k*f^2 for some non-negative integer f."""
    f = 0
    while k * f * f <= D:
        if k * f * f == D:
            return True
        f += 1
    return False


def is_kummer_form(f: tuple[int, int, int]) -> bool:
    a, b, c = f
    return a % 2 == 0 and b % 2 == 0 and c % 2 == 0


def sigma(n: int) -> int:
    """Sum of positive divisors of n (n >= 1)."""
    if n <= 0:
        raise ValueError("sigma needs n >= 1")
    total = 0
    for d in range(1, n + 1):
        if n % d == 0:
            total += d
    return total
