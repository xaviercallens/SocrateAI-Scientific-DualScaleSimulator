"""
Build A = phi_{-2,1}, B = phi_{0,1}, E4, E6, Delta as exact q,y Laurent series
(q-truncated at QMAX, y within +-YMAX), and extract c(D) = coefficients of 2B.

Run: python theta_forms.py [QMAX]
Writes theta_forms_cache.json with all series + the c(D) table + self-checks.
"""
import sys, json
from fractions import Fraction as Fr
from series import (one_plus_x_pow, mul, mul_many, add, scal, const,
                     y_pow_to_l, eta24_over_q, prod_1_minus_qn_power, truncate)

def frac(x):
    return f"{x.numerator}/{x.denominator}" if x.denominator != 1 else str(x.numerator)

def series_to_jsonable(s):
    return {f"{n},{l}": frac(v) for (n, l), v in sorted(s.items())}


def build(QMAX):
    PMAX = 8 * QMAX
    YMAX = 4 * QMAX + 8  # generous cap on |Y|-power

    # ---- ratio3 = theta3(z)/theta3(0) = prod_n (1+P^{8n-4}Y^2)(1+P^{8n-4}Y^-2)/(1+P^{8n-4})^2
    factors3 = []
    for n in range(1, QMAX + 2):
        a = 8 * n - 4
        if a > PMAX:
            break
        factors3.append(one_plus_x_pow(a, 2, Fr(1), 1, PMAX, YMAX))
        factors3.append(one_plus_x_pow(a, -2, Fr(1), 1, PMAX, YMAX))
        factors3.append(one_plus_x_pow(a, 0, Fr(1), -2, PMAX, YMAX))
    ratio3 = mul_many(factors3, PMAX, YMAX)

    # ---- ratio4 = theta4(z)/theta4(0) = prod_n (1-P^{8n-4}Y^2)(1-P^{8n-4}Y^-2)/(1+P^{8n-4})^{-2}... wait sign
    factors4 = []
    for n in range(1, QMAX + 2):
        a = 8 * n - 4
        if a > PMAX:
            break
        factors4.append(one_plus_x_pow(a, 2, Fr(-1), 1, PMAX, YMAX))
        factors4.append(one_plus_x_pow(a, -2, Fr(-1), 1, PMAX, YMAX))
        factors4.append(one_plus_x_pow(a, 0, Fr(-1), -2, PMAX, YMAX))
    ratio4 = mul_many(factors4, PMAX, YMAX)

    # ---- ratio2 = theta2(z)/theta2(0) = (Y+Y^-1)/2 * prod_n (1+P^{8n}Y^2)(1+P^{8n}Y^-2)/(1+P^{8n})^2
    prefactor2 = {(0, 1): Fr(1, 2), (0, -1): Fr(1, 2)}
    factors2 = []
    for n in range(1, QMAX + 1):
        a = 8 * n
        if a > PMAX:
            break
        factors2.append(one_plus_x_pow(a, 2, Fr(1), 1, PMAX, YMAX))
        factors2.append(one_plus_x_pow(a, -2, Fr(1), 1, PMAX, YMAX))
        factors2.append(one_plus_x_pow(a, 0, Fr(1), -2, PMAX, YMAX))
    ratio2 = mul(prefactor2, mul_many(factors2, PMAX, YMAX), PMAX, YMAX)

    ratio2_sq = mul(ratio2, ratio2, PMAX, YMAX)
    ratio3_sq = mul(ratio3, ratio3, PMAX, YMAX)
    ratio4_sq = mul(ratio4, ratio4, PMAX, YMAX)

    # NOTE: the individual ratio_i^2 need NOT separately have integer q-power (theta3/theta4
    # individually mix half-integer P-powers with the q=e^{2 pi i tau} convention); only the
    # SUM 4*(ratio2^2+ratio3^2+ratio4^2) = phi_{0,1} is guaranteed integer-power in q,y.
    # That cancellation (odd-P-power terms cancelling in the sum) is itself the negative
    # control, checked below via y_pow_to_l's assertion on the summed series B_PY.
    B_PY = scal(4, add(ratio2_sq, ratio3_sq, ratio4_sq))
    B = y_pow_to_l(B_PY)

    # ---- A = phi_{-2,1} = - theta1^2/eta^6 (sign fixed so that A|q^0 = y-2+1/y; see note)
    factorsA = []
    for n in range(1, QMAX + 1):
        a = 8 * n
        if a > PMAX:
            break
        factorsA.append(one_plus_x_pow(a, 2, Fr(-1), 2, PMAX, YMAX))
        factorsA.append(one_plus_x_pow(a, -2, Fr(-1), 2, PMAX, YMAX))
        factorsA.append(one_plus_x_pow(a, 0, Fr(-1), -4, PMAX, YMAX))
    prefactorA_raw = {(0, 2): Fr(1), (0, 0): Fr(-2), (0, -2): Fr(1)}  # Y^2-2+Y^-2 = y-2+1/y
    A_raw_PY = mul(scal(-1, prefactorA_raw), mul_many(factorsA, PMAX, YMAX), PMAX, YMAX)  # theta1^2/eta^6
    A_std_PY = scal(-1, A_raw_PY)  # phi_{-2,1} := -theta1^2/eta^6  (standard DMZ convention)
    for (a, b) in A_std_PY:
        assert a % 8 == 0 and b % 2 == 0, "A: fractional power leaked"
    A = y_pow_to_l(A_std_PY)
    A_raw = y_pow_to_l(A_raw_PY)

    # ---- E4, E6 (divisor-sum definitions; exact integers)
    def sigma(k, n):
        return sum(d ** k for d in range(1, n + 1) if n % d == 0)
    E4 = {0: Fr(1)}
    E6 = {0: Fr(1)}
    for n in range(1, QMAX + 1):
        E4[n] = Fr(240 * sigma(3, n))
        E6[n] = Fr(-504 * sigma(5, n))

    # ---- Delta = q * prod (1-q^n)^24
    eta24 = eta24_over_q(QMAX - 1 if QMAX >= 1 else 0)
    Delta = {n + 1: v for n, v in eta24.items() if n + 1 <= QMAX}

    return {
        "QMAX": QMAX, "PMAX": PMAX, "YMAX": YMAX,
        "A": A, "A_theta1sq_over_eta6_raw": A_raw, "B": B, "E4": E4, "E6": E6, "Delta": Delta,
    }


def index1_check(name, s2d, dmax):
    """Check that coefficient depends only on D = 4n - l^2 (weak-Jacobi-form index-1 structure).
    Returns dict D -> {value, agree, samples}."""
    by_D = {}
    for (n, l), v in s2d.items():
        D = 4 * n - l * l
        by_D.setdefault(D, set()).add(v)
    report = {}
    ok = True
    for D, vals in sorted(by_D.items()):
        if D > dmax or D < -8:
            continue
        agree = (len(vals) == 1)
        ok = ok and agree
        report[D] = {"value": frac(next(iter(vals))) if agree else None,
                     "n_distinct_values_seen": len(vals),
                     "all_values": sorted(frac(v) for v in vals) if not agree else None}
    return ok, report


if __name__ == "__main__":
    QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 14
    data = build(QMAX)
    A, B = data["A"], data["B"]

    okA, repA = index1_check("A", A, QMAX)
    okB, repB = index1_check("B", B, QMAX)

    # c(D) := coefficients of 2B
    twoB = scal(2, B)
    cD = {}
    for (n, l), v in twoB.items():
        D = 4 * n - l * l
        cD.setdefault(D, set()).add(v)
    cD_clean = {}
    cD_conflict = False
    for D, vals in cD.items():
        if len(vals) == 1:
            cD_clean[D] = frac(next(iter(vals)))
        else:
            cD_conflict = True
            cD_clean[D] = None

    checks = {
        "A_index1_structure_holds": okA,
        "B_index1_structure_holds": okB,
        "A_q0_leading_term_y_minus2_plus_yinv": {
            "y^1": frac(A.get((0, 1), Fr(0))),
            "const": frac(A.get((0, 0), Fr(0))),
            "y^-1": frac(A.get((0, -1), Fr(0))),
            "expected_pattern": "1, -2, 1 (i.e. A|q0 = y - 2 + 1/y)",
        },
        "e_K3_from_c0_plus_2cm1": {
            "c(0)": cD_clean.get(0), "c(-1)": cD_clean.get(-1),
            "note": "e(K3) = c(0) + 2*c(-1); expect 24 (Euler characteristic of K3) but this is COMPUTED not assumed",
        },
        "cD_no_index1_conflict": (not cD_conflict),
    }

    out = {
        "QMAX": QMAX,
        "checks": checks,
        "A_index1_report_sample": {str(k): v for k, v in list(repA.items())[:12]},
        "B_index1_report_sample": {str(k): v for k, v in list(repB.items())[:12]},
        "cD_table": cD_clean,
        "A_series": series_to_jsonable(A),
        "A_raw_theta1sq_over_eta6": series_to_jsonable(data["A_theta1sq_over_eta6_raw"]),
        "B_series": series_to_jsonable(B),
        "E4_series": {str(k): frac(v) for k, v in data["E4"].items()},
        "E6_series": {str(k): frac(v) for k, v in data["E6"].items()},
        "Delta_series": {str(k): frac(v) for k, v in data["Delta"].items()},
    }
    with open("theta_forms_cache.json", "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps(checks, indent=1))
    print("cD table (D -> c(D)):")
    for D in sorted(cD_clean):
        print(f"  D={D:4d}  c={cD_clean[D]}")
