"""
Part 1: Euler numbers of Hilb^k(K3) from the DMVV product EVALUATED AT z=0 (y=1) by direct expansion,
compared with Goettsche's formula with chi from Track D's exports (NOT circular: the two sides use
different inputs; v2 compared Goettsche with itself).

Run: cd audit/k3t2_rigidity_v3/B-dyons && <venv-python> part1_euler.py 40 8 4
  args: QMAX (theta cache Q<QMAX>), KMAX (Hilb^k for k<=KMAX), QQ (q-truncation of the product)
Committed run: 40 8 4.  Writes part1_euler_results.json.

Method. c(D) = k_norm * cB(D), k_norm = declared input K3_elliptic_genus_factor_k (inputs.json), cB from
our own theta construction of phi_{0,1}. The product is
    prod_{m=1..KMAX} prod_{n=0..QQ} (1 - p^m q^n)^(-e(m,n)),  e(m,n) = sum_t c(4mn - t^2)
(the y=1 specialisation of prod (1-p^m q^n y^t)^(-c(4mn-t^2))). e(m,n) is computed by summing OUR c table
over t; nothing about chi or Goettsche enters. Expansion is exact (integer arithmetic).
Checks: (a) q^j coefficients of p^k for j>=1 vanish (Z(tau,0) is constant); (b) the q^0 coefficient equals
Goettsche prod (1-p^m)^(-chi_D) with chi_D read from Track D exports; (c) e(1,0)=c(0)+2c(-1) vs chi_D.
Selection of the factor: solve k*cB_sum = chi_D exactly (cB_sum = cB(0)+2cB(-1) computed) -> classification
NORMALISATION (the condition contains the target chi); the Euler numbers themselves are
CONDITIONAL_ON_INPUT (declared k). Neighbouring k=1,3 are HYPOTHETICAL controls (no surface with
Z=phi_{0,1}, chi=12 in this DMVV class is claimed; not counted as discrimination).
"""
import json, sys
from fractions import Fraction as Fr
from common import HERE, load_cache, load_cB, load_inputs, find_trackD_exports


def factor_series(e, m, n, KMAX, QQ):
    """(1 - x)^(-e), x = p^m q^n, as dict {(a,b): coeff}, exact."""
    out = {(0, 0): Fr(1)}
    c = Fr(1)
    j = 1
    while m * j <= KMAX and n * j <= QQ:
        c = c * (e + j - 1) / j
        if c != 0:
            out[(m * j, n * j)] = c
        j += 1
        if n == 0 and m * j > KMAX:
            break
    return out


def mul(s1, s2, KMAX, QQ):
    out = {}
    for (a1, b1), v1 in s1.items():
        for (a2, b2), v2 in s2.items():
            a, b = a1 + a2, b1 + b2
            if a > KMAX or b > QQ:
                continue
            out[(a, b)] = out.get((a, b), 0) + v1 * v2
    return {k: v for k, v in out.items() if v != 0}


def product_at_y1(cD, KMAX, QQ):
    """cD: dict D -> exact c(D); requires c(D)=0 for D<=-2 (checked by caller)."""
    S = {(0, 0): Fr(1)}
    exps = {}
    for m in range(1, KMAX + 1):
        for n in range(0, QQ + 1):
            top = 4 * m * n
            e = Fr(0)
            for t in range(-int((top + 1) ** 0.5) - 1, int((top + 1) ** 0.5) + 2):
                e += cD.get(top - t * t, Fr(0))
            exps[(m, n)] = e
            if e != 0:
                S = mul(S, factor_series(e, m, n, KMAX, QQ), KMAX, QQ)
    return S, exps


def goettsche(chi, KMAX):
    """[p^k] prod_m (1-p^m)^(-chi)."""
    s = {0: Fr(1)}
    for m in range(1, KMAX + 1):
        f = {0: Fr(1)}
        c = Fr(1)
        j = 1
        while m * j <= KMAX:
            c = c * (chi + j - 1) / j
            f[m * j] = c
            j += 1
        new = {}
        for a, v in s.items():
            for b, w in f.items():
                if a + b <= KMAX:
                    new[a + b] = new.get(a + b, 0) + v * w
        s = new
    return [s.get(k, Fr(0)) for k in range(KMAX + 1)]


if __name__ == "__main__":
    QMAX, KMAX, QQ = int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])
    inp = load_inputs()
    knorm = Fr(inp["K3_elliptic_genus_factor_k"]["value"])
    cache = load_cache(QMAX)
    cB = load_cB(cache)
    assert 4 * KMAX * QQ <= 4 * QMAX - 3, "theta cache too small for requested product range"
    small_D = {D: v for D, v in cB.items() if D <= -2 and v != 0}
    cB_sum = cB[0] + 2 * cB[-1]
    tag, dpath, dexp = find_trackD_exports()
    chi_key = "chi" if "chi" in dexp else "chi_K3_resolved"
    chi_D = Fr(dexp[chi_key])

    def run(k):
        cD = {D: k * v for D, v in cB.items()}
        S, exps = product_at_y1(cD, KMAX, QQ)
        q0 = [S.get((kk, 0), Fr(0)) for kk in range(KMAX + 1)]
        qpos = {f"p{kk}q{j}": str(S.get((kk, j), Fr(0))) for kk in range(KMAX + 1) for j in range(1, QQ + 1) if S.get((kk, j), Fr(0)) != 0}
        return q0, qpos, exps

    q0, qpos, exps = run(knorm)
    gt = goettsche(chi_D, KMAX)
    match = (q0 == gt)
    # exponents e(m,n): only n=0 should be nonzero
    nonzero_exps_n_pos = {f"{m},{n}": str(v) for (m, n), v in exps.items() if n > 0 and v != 0}
    # exact solve of k from chi: k * cB_sum = chi_D
    k_solved = chi_D / cB_sum
    # controls (HYPOTHETICAL for k not realised): every integer k=1..6, uncentred scan
    scan = {}
    for k in range(1, 7):
        q0k, _, _ = run(Fr(k))
        scan[str(k)] = {"e_Hilb_k": [str(x) for x in q0k],
                         "matches_goettsche_with_chi_D": q0k == gt,
                         "chi_implied_e(Hilb^1)": str(q0k[1]) if KMAX >= 1 else None,
                         "label": "HYPOTHETICAL control (only k=2 is a declared K3 input)" if k != int(knorm) else "declared input"}
    out = {
        "args": {"QMAX": QMAX, "KMAX": KMAX, "QQ": QQ},
        "trackD_exports_source": {"version_dir": tag, "path": dpath, "chi_key_used": chi_key, "chi": str(chi_D), "trackD_status_field": dexp.get("status"), "trackD_source_field": dexp.get("source")},
        "c_D_le_-2_nonzero (must be empty)": {str(D): str(v) for D, v in small_D.items()},
        "cB(0),cB(-1),cB_sum=cB(0)+2cB(-1)": [str(cB[0]), str(cB[-1]), str(cB_sum)],
        "k_used (declared input K3_elliptic_genus_factor_k)": str(knorm),
        "chi_from_own_series = k*cB_sum": str(knorm * cB_sum),
        "k_solved_exactly_from_chi_D = chi_D/cB_sum": str(k_solved),
        "declared_k_equals_solved_k": knorm == k_solved,
        "e(m,n) at y=1 nonzero for n>=1 (must be empty)": nonzero_exps_n_pos,
        "e(m,0) for m<=KMAX": {str(m): str(exps[(m, 0)]) for m in range(1, KMAX + 1)},
        "q^j (j>=1) coefficients of p^k, k<=KMAX (must be empty)": qpos,
        "euler_numbers_from_DMVV_at_y1_q0": [str(x) for x in q0],
        "goettsche_with_chi_from_trackD": [str(x) for x in gt],
        "product_equals_goettsche": match,
        "k_scan_1_to_6_all_integers": scan,
        "expected_field_filled_after_computing": {
            "e(Hilb^k(K3)) k=0..4": "1, 24, 324, 3200, 25650 (Goettsche; FROM MEMORY, literature check only)",
        },
        "shared_inputs": ["K3_elliptic_genus_factor_k", "DMVV_product_formula", "Jacobi_theta_definitions",
                          "Sym_k_Euler_numbers_comparison_formula", "Track_D_chi_K3"],
        "rigidity": {
            "k": {"parameter_inserted": "overall factor k in Z=k*phi_{0,1}",
                  "selecting_condition": "k*(cB(0)+2cB(-1)) = chi_K3 from Track D (contains the target chi)",
                  "condition_uses_true_value": True, "solution_set": f"{{{k_solved}}} exact solve",
                  "classification": "NORMALISATION",
                  "control": "k=1,3,4,5,6 change every e(Hilb^k); HYPOTHETICAL (not realizable surfaces), not counted as discrimination",
                  "control_perturbs_same_parameter": True},
            "euler_numbers": {"classification": "CONDITIONAL_ON_INPUT", "input_it_depends_on": "K3_elliptic_genus_factor_k, Track D chi",
                              "note": "given DMVV and Z=2phi_{0,1}, the whole sequence is computed; agreement with Goettsche(chi_D) is a consistency check of DMVV expansion + chi_D"},
        },
    }
    with open(HERE / "part1_euler_results.json", "w") as f:
        json.dump(out, f, indent=1)
    print(json.dumps({k: out[k] for k in ("declared_k_equals_solved_k", "product_equals_goettsche",
          "euler_numbers_from_DMVV_at_y1_q0", "q^j (j>=1) coefficients of p^k, k<=KMAX (must be empty)",
          "e(m,n) at y=1 nonzero for n>=1 (must be empty)", "trackD_exports_source")}, indent=1))
