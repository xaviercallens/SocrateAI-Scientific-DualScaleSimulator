"""
Reverse pass, predictions P3 and P4.

Lean theorems extended (all three test the SAME quantity, `polarDefect m coef`,
under different names in DualScaleDyons/DMVV.lean and ImmortalHigher.lean):
  * `polar_part_removes_pole`, `polar_coefficient_pinned`  (DMVV.lean): m in [1,2,3]
  * `immortal_exact_m23`                                   (ImmortalHigher.lean): m in {2,3}
all through q^2 (`dmvv 4 2`). "Not proved": m >= 4.

polarDefect(m, coef) := G_{m+1} - coef * A * A2m(m), with A2m(m) built from
A2m(m) = A2mRest(m) [the s!=0 strip terms, DMZ (9.55)] + the s=0 term
y/(1-y)^2 (whose product with A is `rPart`, per the doc comment: A*y/(1-y)^2
= R). `vanish2` checks the defect is 0 and has zero y-derivative at y=1
(cancels a double pole), coefficient-wise in q.

PREDICTION: for m = 4, 5 (beyond Lean's m<=3):
  (a) vanish2(polarDefect(m, p24(m+1))) holds through q^2 [extends
      polar_part_removes_pole / immortal_exact_m23]
  (b) vanish2(polarDefect(m, p24(m+1)+1)) and vanish2(polarDefect(m, p24(m+1)-1))
      FAIL [extends polar_coefficient_pinned's negative control]

VALIDATION: reproduce Lean's m=1,2,3 results first with this independent
(n,l)-Laurent-series engine (fresh from the defining formulas, not copied
from Lean or from the sibling P,Y-power engine in dyons/series.py).
"""
import json
from fractions import Fraction as Fr

from nlseries import add, sub, scal, mul, clean, truncate
from dmvv_product import load_cD, dmvv, p24_direct

CACHE = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/dyons/theta_forms_cache.json"


def pf(s):
    a_b = s.split("/") if isinstance(s, str) and "/" in s else None
    if a_b:
        return Fr(int(a_b[0]), int(a_b[1]))
    return Fr(int(s))


def load_A(qmax_needed):
    cache = json.load(open(CACHE))
    assert cache["QMAX"] >= qmax_needed
    A = {}
    for k, v in cache["A_series"].items():
        n, l = k.split(",")
        n = int(n)
        if n <= qmax_needed:
            A[(n, int(l))] = pf(v)
    return A


def r_part(Q):
    """Prod_{n=1..Q} (1 - y q^n)^2 (1 - y^-1 q^n)^2 (1 - q^n)^{-4}, truncated q^Q."""
    out = {(0, 0): Fr(1)}
    for n in range(1, Q + 1):
        f1 = _monomial_pow(Q, n, 1, -1, 2)   # (1 - q^n y)^2
        f2 = _monomial_pow(Q, n, -1, -1, 2)  # (1 - q^n y^-1)^2
        f3 = _monomial_pow(Q, n, 0, -1, -4)  # (1 - q^n)^-4
        for f in (f1, f2, f3):
            out = mul(out, f, Q)
    return out


def _monomial_pow(qmax, n, l, coeff, k):
    """(1 + coeff q^n y^l)^k truncated at q^qmax, k any integer."""
    from math import comb
    out = {(0, 0): Fr(1)}
    jmax = qmax // n
    for j in range(1, jmax + 1):
        qn = n * j
        if qn > qmax:
            break
        if k >= 0:
            if j > k:
                continue
            b = Fr(comb(k, j))
        else:
            b = Fr((-1) ** j * comb(-k + j - 1, j))
        out[(qn, l * j)] = out.get((qn, l * j), Fr(0)) + b * (Fr(coeff) ** j)
    return clean(out)


def a2m955(m, Q):
    """DMZ (9.55): the strip (|q|<|y|<1) Fourier expansion of A_{2,m} minus its
    s=0 term, Sum_{s!=0,t} 1/2(sgn t - sgn(s+eps)) t q^{ms^2-st} y^{2ms-t}."""
    S = {}
    for s in range(1, Q + 3):
        for t in range(1, Q + 3):
            qpow = m * s * s + s * t
            if qpow <= Q:
                S[(qpow, 2 * m * s + t)] = S.get((qpow, 2 * m * s + t), Fr(0)) + t
                S[(qpow, -2 * m * s - t)] = S.get((qpow, -2 * m * s - t), Fr(0)) + t
    return clean(S)


def a_times_a2m(m, Q, A):
    return add(r_part(Q), mul(A, a2m955(m, Q), Q))


def polar_defect(Gm1, m, coef, Q, A):
    return sub(truncate(Gm1, Q), scal(coef, a_times_a2m(m, Q, A)))


def vanish2(s, Q):
    """s.all fun p => p.eval1==0 && p.d1==0, coefficientwise in q, n=0..Q.
    HONESTY NOTE (matches DualScaleDyons/DMVV.lean's own doc comment under
    polar_part_removes_pole verbatim): polarDefect is built to be symmetric
    under y <-> y^-1 (A, A2m and the dmvv product all have that symmetry), so
    d1 = sum_l l*coeff(l) is an ODD functional of an even series and is 0
    identically, for every coef, including the +-1 negative controls. It never
    discriminates anything in this test. The only condition doing real work is
    eval1 (the q-independent pole coefficient equalling p24(m+1))."""
    per_n = {}
    for (n, l), v in s.items():
        if n > Q:
            continue
        e, d = per_n.get(n, (Fr(0), Fr(0)))
        per_n[n] = (e + v, d + l * v)
    ok = all(e == 0 and d == 0 for e, d in per_n.values())
    d1_ever_nonzero = any(d != 0 for e, d in per_n.values())
    bad = {n: (str(e), str(d)) for n, (e, d) in per_n.items() if not (e == 0 and d == 0)}
    return ok, bad, d1_ever_nonzero


def main():
    Q = 2
    cD, dmax = load_cD()
    A = load_A(Q)
    p24 = p24_direct(8)

    # need G_{m+1} for m=1..5, i.e. K up to 6
    G = dmvv(6, Q, cD, dmax)

    out = {"Q": Q, "validation_m1_2_3": {}, "P3_m4": {}, "P4_m5": {}}

    for m in (1, 2, 3):
        row = {}
        ok0, bad0, d1nz0 = vanish2(polar_defect(G[m + 1], m, p24[m + 1], Q, A), Q)
        okp, badp, d1nzp = vanish2(polar_defect(G[m + 1], m, p24[m + 1] + 1, Q, A), Q)
        okm, badm, d1nzm = vanish2(polar_defect(G[m + 1], m, p24[m + 1] - 1, Q, A), Q)
        row["vanish2_at_p24"] = {"holds": ok0, "bad": bad0}
        row["vanish2_at_p24_plus1_should_fail"] = {"holds_ie_BAD_if_true": okp, "bad": badp}
        row["vanish2_at_p24_minus1_should_fail"] = {"holds_ie_BAD_if_true": okm, "bad": badm}
        row["matches_lean_pattern"] = (ok0 is True) and (okp is False) and (okm is False)
        out["validation_m1_2_3"][f"m={m}"] = row

    all_validation_ok = all(out["validation_m1_2_3"][f"m={m}"]["matches_lean_pattern"] for m in (1, 2, 3))
    out["all_validation_matches_lean"] = all_validation_ok

    for m, key in ((4, "P3_m4"), (5, "P4_m5")):
        ok0, bad0, d1nz0 = vanish2(polar_defect(G[m + 1], m, p24[m + 1], Q, A), Q)
        okp, badp, d1nzp = vanish2(polar_defect(G[m + 1], m, p24[m + 1] + 1, Q, A), Q)
        okm, badm, d1nzm = vanish2(polar_defect(G[m + 1], m, p24[m + 1] - 1, Q, A), Q)
        out[key] = {
            "p24_m_plus_1": str(p24[m + 1]),
            "vanish2_at_p24": {"holds": ok0, "bad_coefficients": bad0},
            "vanish2_at_p24_plus1": {"holds": okp, "note": "predicted to be FALSE (negative control)"},
            "vanish2_at_p24_minus1": {"holds": okm, "note": "predicted to be FALSE (negative control)"},
            "prediction_confirmed": (ok0 is True) and (okp is False) and (okm is False),
            "lean_checked_range": "m in {1,2,3} only (polar_part_removes_pole / polar_coefficient_pinned / immortal_exact_m23)",
            "honesty_note_d1": (
                "d1 (the y-derivative-at-1 half of vanish2) was 0 identically in every case here "
                "(main and both negative controls: d1_ever_nonzero={}), by the y<->y^-1 symmetry "
                "of A, A2m and the dmvv product -- exactly the caveat DMVV.lean's own doc comment "
                "states under polar_part_removes_pole. The only condition doing real discriminating "
                "work is eval1: G_{{m+1}}(q^n,y=1) = p24(m+1) for every n<=Q, which the +-1 controls "
                "break. This is NOT an independently-verified double-pole cancellation; it is a "
                "single (first-order) pole-coefficient check."
            ).format(d1nz0 or d1nzp or d1nzm),
        }

    with open("/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/reverse/p3_p4_results.json", "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(json.dumps(out, indent=2, default=str))


if __name__ == "__main__":
    main()
