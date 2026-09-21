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


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        sys.exit(0 if _self_test() else 1)
    report()
