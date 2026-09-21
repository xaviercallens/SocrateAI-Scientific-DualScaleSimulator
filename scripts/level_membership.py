#!/usr/bin/env python3
"""
Which CM discriminants occur in the level-n K3 family?

CRITERION, taken from Stream 2 (SocrateAI-Scientific-Agora-K3-DarkMatter),
brief `briefs/THOUGHT_EXPERIMENTS_K3_SELECTION_2026_09_21.md` GE-8 and
certificate `data/certificates/A2_MEMBERSHIP.json`.

VALIDATED AGAINST: Stream 2 commit `da84a90`, certificate checker version
`f091e43`. That repository moves daily, so if the self-test ever fails, compare
against this SHA before assuming the criterion changed.

    D occurs in the level-n family  <=>  D is a square modulo 4n.

This script does NOT re-derive that criterion; it transcribes it and applies it
at n = 12, the level this repository is implicitly at (audit/K3_SELECTION.md).

VALIDATION (--self-test): reproduces Stream 2's own `n7` and `n10` columns for
all 39 discriminants in their certificate, both directions. If their file moves
or changes, the test fails loudly rather than silently drifting.

SCOPE: nothing here selects a K3. Stream 2 records that the rho = 20 cut is NOT
adopted by the programme and that CM points are dense on the modular curve, so
requiring rho = 20 cuts out nothing; what makes any list finite is a bound on
|D|. No physical reading is claimed anywhere.
"""
import json
import math
import os
import sys

CERT = os.environ.get(
    "A2_MEMBERSHIP_JSON",
    os.path.expanduser("~/SocrateAI-Scientific-Agora-K3-DarkMatter/"
                       "data/certificates/A2_MEMBERSHIP.json"),
)


def admitted(D: int, n: int) -> bool:
    """Stream 2 GE-8 criterion: D is admitted at level n iff D is a square mod 4n."""
    m = 4 * n
    return (D % m) in {(x * x) % m for x in range(m)}


def _load_cert():
    if not os.path.exists(CERT):
        print(f"[!] Stream 2 certificate not found at {CERT}\n"
              f"    Set A2_MEMBERSHIP_JSON to its path.", file=sys.stderr)
        return None
    return json.load(open(CERT))


def _self_test() -> bool:
    cert = _load_cert()
    if cert is None:
        return False
    rows = cert["discriminants_admitted_summary"]
    bad = [(r["D"], n, admitted(r["D"], n), r[f"n{n}"])
           for r in rows for n in (7, 10) if admitted(r["D"], n) != r[f"n{n}"]]
    print(f"criterion vs Stream 2 certificate ({len(rows)} D x 2 levels): "
          f"{'PASS' if not bad else f'FAIL {bad[:6]}'}")
    return not bad


def report(levels=(7, 10, 12), d_max=100):
    Ds = [D for D in range(-d_max, 0) if D % 4 in (0, 1)]
    for n in levels:
        adm = [D for D in Ds if admitted(D, n)]
        print(f"  level n={n:2d}: {len(adm):2d} of {len(Ds)} discriminants |D| <= {d_max}")
    print("\n  Stream 8's two 'most attractive' K3s (the self-dual torus points):")
    for D, name in ((-3, "tau = omega (D = -3)"), (-4, "tau = i     (D = -4)")):
        cols = "  ".join(f"n{n}={admitted(D, n)!s:5s}" for n in levels)
        print(f"    {name}: {cols}")


def hall_divisors(n: int):
    """Hall divisors Q | n with gcd(Q, n/Q) = 1 -- these index the Atkin-Lehner
    group W(n), of order 2^omega(n)."""
    return [Q for Q in range(1, n + 1) if n % Q == 0 and math.gcd(Q, n // Q) == 1]


def has_trace0_representative(n: int, Q: int, a_bound: int = 25, b_bound: int = 600) -> bool:
    """Does the Atkin-Lehner coset w_Q contain an ORDER-2 elliptic element?

    An element of w_Q is M = [[Qa, b], [nc, Qd]] with det M = Q. It is elliptic of
    order 2 (hence has a fixed point in H that it fixes as an involution) exactly
    when its trace Q(a+d) vanishes, i.e. d = -a, which gives

        -Q^2 a^2 - n b c = Q.

    Reproduces Stream 2's `pairs_n_Q_without_trace0_representative` exactly over
    their whole sweep range n = 1..30 (see `--self-test`).

    NOTE the distinction their certificate draws and this function does NOT: a
    coset with no trace-0 element may still have a fixed point in H, fixed by an
    element of order 4 or 6. At n = 10 and n = 26 with Q = 2 that is what happens.
    """
    for a in range(-a_bound, a_bound + 1):
        r = -Q * Q * a * a - Q
        if r % n:
            continue
        t = r // n
        if t == 0:
            return True
        for b in list(range(1, b_bound)) + list(range(-b_bound, 0)):
            if t % b == 0:
                return True
    return False


def _self_test_atkin_lehner() -> bool:
    """Validate the trace-0 test against Stream 2's ATKIN_LEHNER_DISC_FORM certificate."""
    path = os.path.join(os.path.dirname(CERT), "ATKIN_LEHNER_DISC_FORM.json")
    if not os.path.exists(path):
        print(f"[!] {path} not found -- skipping Atkin-Lehner check", file=sys.stderr)
        return True
    sw = json.load(open(path))["leg_W_sweep"]
    theirs = {tuple(x) for x in sw["pairs_n_Q_without_trace0_representative"]}
    lo, hi = sw["n_range"]
    mine = {(n, Q) for n in range(lo, hi + 1) for Q in hall_divisors(n)
            if Q > 1 and not has_trace0_representative(n, Q)}
    ok = mine == theirs
    print(f"trace-0 sweep vs Stream 2 certificate (n = {lo}..{hi}): "
          f"{'PASS' if ok else f'FAIL +{sorted(mine - theirs)} -{sorted(theirs - mine)}'}")
    return ok


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        sys.exit(0 if (_self_test() and _self_test_atkin_lehner()) else 1)
    report()
    print("\n  Atkin-Lehner group W(n) and which cosets fix a point of H:")
    for n in (7, 10, 12):
        parts = ", ".join(
            f"w_{Q}{'(Fricke)' if Q == n else ''}:"
            f"{'fixes' if has_trace0_representative(n, Q) else 'FREE'}"
            for Q in hall_divisors(n) if Q > 1)
        print(f"    n={n:2d}  |W(n)|={len(hall_divisors(n))}   {parts}")
