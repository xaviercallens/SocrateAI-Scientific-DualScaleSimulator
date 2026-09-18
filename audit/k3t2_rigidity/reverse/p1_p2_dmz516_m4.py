"""
Reverse pass, predictions P1 and P2.

P1 (Goettsche numbers beyond Lean's k<=5, DualScaleDyons/DMVV.lean `p24_values`
/`goettsche`): predict p24(6), p24(7), p24(8) two independent exact ways and
cross-check.

P2 (DMZ (5.16) at m=4, docstring-quoted in DualScaleDyons/DMVV.lean but only
Tier-A-checked through q^1 in `dmz_516_q1`; `dmz_516_q2` explicitly stops at
k<=4, i.e. m<=3, EXCLUDING the m=4/k=5 row): predict that
    72*G_5 = 51*B^5 + 155*E4*A^2*B^3 + 93*E6*A^3*B^2 + 102*E4^2*A^4*B + 31*E4*E6*A^5
holds through q^2, q^3, q^4 (beyond Lean's q^1-only coverage of this row),
where G_5 is the k=5 coefficient of the independently re-implemented DMVV
product (dmvv_product.py), and A, B, E4, E6 are loaded from our own
theta-function computation (dyons/theta_forms_cache.json), not from Lean.

VALIDATION step (run first, reported honestly): reproduce Lean's already-
Tier-A dmz_516_q1 (k<=5, q^1) and dmz_516_q2 (k<=4, q^2) with this
independent engine, as a negative-control-style sanity check before trusting
the k=5/q^2..q^4 extension.
"""
import json
from fractions import Fraction as Fr

from nlseries import add, sub, scal, mul, mul_many, powr, from_1d, truncate
from dmvv_product import load_cD, dmvv, p24_direct

CACHE = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/dyons/theta_forms_cache.json"


def pf(s):
    if s is None:
        return None
    if "/" in s:
        a, b = s.split("/")
        return Fr(int(a), int(b))
    return Fr(int(s))


def load_ABE4E6(qmax_needed):
    cache = json.load(open(CACHE))
    assert cache["QMAX"] >= qmax_needed, f"cache QMAX={cache['QMAX']} < needed {qmax_needed}"
    A, B = {}, {}
    for k, v in cache["A_series"].items():
        n, l = k.split(",")
        n = int(n)
        if n <= qmax_needed:
            A[(n, int(l))] = pf(v)
    for k, v in cache["B_series"].items():
        n, l = k.split(",")
        n = int(n)
        if n <= qmax_needed:
            B[(n, int(l))] = pf(v)
    E4_1d = {int(k): pf(v) for k, v in cache["E4_series"].items() if int(k) <= qmax_needed}
    E6_1d = {int(k): pf(v) for k, v in cache["E6_series"].items() if int(k) <= qmax_needed}
    return A, B, from_1d(E4_1d), from_1d(E6_1d)


def ser_eq(s1, s2, qmax):
    s1t, s2t = truncate(s1, qmax), truncate(s2, qmax)
    return s1t == s2t


def diff_report(s1, s2, qmax):
    s1t, s2t = truncate(s1, qmax), truncate(s2, qmax)
    keys = set(s1t) | set(s2t)
    mism = []
    for k in sorted(keys):
        v1, v2 = s1t.get(k, Fr(0)), s2t.get(k, Fr(0))
        if v1 != v2:
            mism.append({"n_l": list(k), "lhs": str(v1), "rhs": str(v2)})
    return mism


def main():
    QMAX_NEEDED = 6
    cD, dmax = load_cD()
    A, B, E4, E6 = load_ABE4E6(QMAX_NEEDED)

    out = {"P1_goettsche_k6_7_8": {}, "P2_dmz516_m4": {}, "validation": {}}

    # ---------------- P1: Goettsche numbers p24(6), p24(7), p24(8) ----------------
    p24_conv = p24_direct(8)
    G_y1 = dmvv(8, 0, cD, dmax)  # K=8, Q=0: only needs c(D<=32); well within dmax=88
    p24_dmvv = []
    for k in range(9):
        s = sum(v for (n, l), v in G_y1[k].items() if n == 0)
        p24_dmvv.append(s)
    # Negative control: the q^0,y=1 restriction of dmvv collapses (derived, not
    # assumed -- only s=0 factors survive q^Q at Q=0) to prod_r(1-p^r)^{-(c0+2c(-1))};
    # perturbing the exponent 24=c(0)+2c(-1) to 23 or 25 must give a DIFFERENT
    # table, showing the value 24 (hence c(0)=20, c(-1)=2) is load-bearing.
    # NOTE (honesty): method A and method B are the SAME convolution once the
    # exponent 24 is fixed -- method B additionally derives that the exponent
    # equals c(0)+2c(-1), so this is "one route plus an independent check that
    # the exponent is 24", not two independent routes to the numbers themselves.
    from dmvv_product import comb_pos

    def euler_numbers_exponent(exponent, nmax):
        s = [Fr(0)] * (nmax + 1)
        s[0] = Fr(1)
        for r in range(1, nmax + 1):
            factor = {}
            j = 0
            while r * j <= nmax:
                factor[r * j] = Fr(comb_pos(exponent, j))
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

    exp_correct = 24
    neg_ctrl = {
        "exponent_23": [str(x) for x in euler_numbers_exponent(23, 8)],
        "exponent_25": [str(x) for x in euler_numbers_exponent(25, 8)],
        "exponent_24_matches_p24": euler_numbers_exponent(24, 8) == p24_conv,
        "note": "c(0)=20, c(-1)=2 (Lean's c_first) give exponent c(0)+2c(-1)=24 exactly; "
                "23 or 25 give a different, wrong table -- the exponent is load-bearing.",
    }

    out["P1_goettsche_k6_7_8"] = {
        "route": "prod_k(1-q^k)^{-24}, coefficient k=6,7,8, by direct exact convolution",
        "cross_check": "the q^0,y=1 restriction of the independently re-implemented DMVV product "
                        "(dmvv_product.py) reduces, by derivation (only s=0 survives at Q=0), to the "
                        "SAME product with exponent c(0)+2c(-1); both give identical tables (see "
                        "two_methods_agree_k0_8) -- this is one route plus a structural derivation, "
                        "not two independent numerical routes",
        "method_A_direct_partition_convolution": {str(k): str(p24_conv[k]) for k in range(9)},
        "method_B_dmvv_q0_y1_restriction": {str(k): str(p24_dmvv[k]) for k in range(9)},
        "two_methods_agree_k0_8": p24_conv == p24_dmvv,
        "negative_control_exponent_perturbation": neg_ctrl,
        "lean_checked_range": "k<=5 (p24_values, goettsche)",
        "new_predictions": {
            "p24_6": str(p24_conv[6]),
            "p24_7": str(p24_conv[7]),
            "p24_8": str(p24_conv[8]),
        },
        "tier": "B (exact arithmetic, negative control on the exponent)",
    }

    # ---------------- Validation: reproduce dmz_516_q1 (k<=5, q^1) and dmz_516_q2 (k<=4, q^2) ----------------
    G = dmvv(5, 4, cD, dmax)  # K=5 (covers m up to 4), Q=4 (covers q^1..q^4)

    def P(xs, qmax):
        return truncate(mul_many(xs, qmax), qmax)

    rows = []
    for k, (d, rhs_builder) in enumerate([
        (1, lambda qm: {(0, 0): Fr(1)}),
        (1, lambda qm: scal(2, B)),
        (4, lambda qm: add(scal(9, P([B, B], qm)), scal(3, P([E4, A, A], qm)))),
        (27, lambda qm: add(add(scal(50, P([B, B, B], qm)), scal(48, P([E4, A, A, B], qm))),
                             scal(10, P([E6, A, A, A], qm)))),
        (384, lambda qm: add(add(scal(475, P([B, B, B, B], qm)), scal(886, P([E4, A, A, B, B], qm))),
                              add(scal(360, P([E6, A, A, A, B], qm)), scal(199, P([E4, E4, A, A, A, A], qm))))),
    ]):
        for qm in (1, 2):
            lhs = scal(d, G[k])
            rhs = rhs_builder(qm)
            ok = ser_eq(lhs, rhs, qm)
            rows.append({"k": k, "q_order": qm, "holds": ok,
                         "lean_theorem": "dmz_516_q1" if qm == 1 else "dmz_516_q2",
                         "lean_covers_this_k_at_this_order": (qm == 1) or (qm == 2 and k <= 4)})
    out["validation"]["dmz516_k0_4_reproduced"] = rows
    out["validation"]["all_hold"] = all(r["holds"] for r in rows)

    # ---------------- P2: DMZ (5.16) m=4 (k=5) beyond Lean's q^1-only coverage ----------------
    G5 = G[5]
    lhs = scal(72, G5)
    coeffs = [51, 155, 93, 102, 31]
    term_specs = [[B, B, B, B, B], [E4, A, A, B, B, B], [E6, A, A, A, B, B],
                  [E4, E4, A, A, A, A, B], [E4, E6, A, A, A, A, A]]

    def rhs_with(coeffs, qm):
        terms = [P(spec, qm) for spec in term_specs]
        r = {}
        for c, t in zip(coeffs, terms):
            r = add(r, scal(c, t))
        return r

    for qm in (1, 2, 3, 4):
        rhs = rhs_with(coeffs, qm)
        ok = ser_eq(lhs, rhs, qm)
        mism = [] if ok else diff_report(lhs, rhs, qm)
        entry = {
            "holds": ok,
            "lean_checked_this_order_for_m4": (qm == 1),
            "mismatches": mism[:20],
        }
        if qm == 2:
            # Negative control at the order actually beyond Lean's coverage:
            # perturb each of the 5 RHS coefficients by +1 and confirm the
            # identity breaks (mirrors dyons/part2_psi_m.py's c(D)+1 control).
            neg = []
            for i in range(5):
                c2 = list(coeffs)
                c2[i] += 1
                mism_i = diff_report(lhs, rhs_with(c2, qm), qm)
                neg.append({"perturbed_coeff_index": i, "orig": coeffs[i], "perturbed": c2[i],
                            "n_mismatched_coeffs": len(mism_i), "identity_breaks": len(mism_i) > 0})
            entry["negative_control_perturb_each_coeff_by_plus1"] = neg
        out["P2_dmz516_m4"][f"q^{qm}"] = entry
    out["P2_dmz516_m4"]["provenance"] = ("RHS formula for m=4 (72 Delta psi_4 = 51 A^-1 B^5 + 155 E4 A B^3 "
        "+ 93 E6 A^2 B^2 + 102 E4^2 A^3 B + 31 E4 E6 A^4) is Tier L, quoted verbatim from DMZ (5.16) via "
        "DualScaleDyons/DMVV.lean's doc comment (arXiv:1208.4074); Lean's own dmz_516_q1 checks it at q^1 "
        "(Tier A), dmz_516_q2 explicitly EXCLUDES k=5/m=4 at q^2. The q^2..q^4 check here is Tier B "
        "(exact arithmetic against an independently re-implemented DMVV product), not a re-derivation "
        "of the RHS ansatz itself.")

    with open("/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/reverse/p1_p2_results.json", "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(json.dumps(out, indent=2, default=str)[:6000])


if __name__ == "__main__":
    main()
