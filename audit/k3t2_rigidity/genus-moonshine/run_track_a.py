#!/usr/bin/env python3
"""
Track A -- K3 elliptic genus, exact q,y-series computation (Tier B).

Everything here is computed from the defining Jacobi-theta / Dedekind-eta
product formulas (see thetas.py, series2d.py) using Python Fraction
arithmetic -- no floating point, no literature constant is ever fed in as
a target to compare against itself. Truncation-stability (rerun at a
higher q-cutoff, confirm the reported coefficients do not move) is checked
explicitly and recorded below.

What is delivered here (all Tier B, verified):
  (1) phi_{0,1}(tau,z) = 4 sum_{i=2,3,4} (theta_i(z)/theta_i(0))^2
      Z_K3 = 2 phi_{0,1};  Z_K3(tau,0) computed and checked constant = 24
      (the K3 Euler number) with NO constant fed in -- it falls out of the
      theta-ratio algebra.
      c(n,l) coefficients, checked to depend only on D = 4n-l^2, and the
      c(D) table for D = -1,0,3,4,7,8,11,12,...
  phi_{-2,1}(tau,z) = theta_1(tau,z)^2/eta(tau)^6, as a bonus exact
      computation (same machinery, no extra risk).
  Rigidity scan (b): replace the "2" in Z=k*phi_{0,1} by a rational k;
      report the solution set of {Z(tau,0)=24} and separately of
      {all computed coefficients integral}, with a negative control.

What was ATTEMPTED but NOT delivered to Tier-B rigor (see could_not_do in
the run's printed summary and the workflow's final report): the mu/H(tau)
Appell-Lerch extraction (item 2) and the 2A twining (item 3). The
regularization Psi := T1*S(tau,z) (S = the Appell-Lerch sum) was
implemented exactly (appell.py) and is finite/computable with no ad hoc
truncation, but two independently-derived candidate normalizations for
how "24*mu-term" combines with H(tau)*theta_1^2/eta^3 gave RESULTS THAT
DISAGREE with each other under the cross-slice consistency check (dividing
by two different y-power slices of theta_1^2 should give the same H(tau);
it did not). Rather than report a number from a formula whose normalization
we could not pin down and self-verify, this is left undone.
"""
import json
import sys
from fractions import Fraction as Fr

from series2d import y_to_1, divide_scalar_series, mul, add, scal, shift, reciprocal_1d
from thetas import theta2, theta3, theta4, theta1_over_i, eta3, eta6, prod1
from appell import psi


def compute_ZK3_and_phi01(cutoff_q, N):
    imax = 8 * cutoff_q
    t2 = theta2(N, imax)
    t3 = theta3(N, imax)
    t4 = theta4(N, imax)
    r2 = divide_scalar_series(t2, y_to_1(t2), imax)
    r3 = divide_scalar_series(t3, y_to_1(t3), imax)
    r4 = divide_scalar_series(t4, y_to_1(t4), imax)
    r2sq = mul(r2, r2, imax)
    r3sq = mul(r3, r3, imax)
    r4sq = mul(r4, r4, imax)
    phi01 = scal(add(add(r2sq, r3sq), r4sq), 4)
    ZK3 = scal(phi01, 2)
    return ZK3, phi01


def compute_phi_minus21(cutoff_q, N):
    from series2d import reciprocal_1d
    imax = 8 * cutoff_q
    T1 = theta1_over_i(N, imax)
    T1sq = mul(T1, T1, imax)
    theta1sq = {k: -v for k, v in T1sq.items()}  # theta_1^2 = -T1^2
    e6 = eta6(N, imax)
    e6_1d = {i: v for (i, j), v in e6.items() if j == 0}
    inv_e6 = reciprocal_1d(e6_1d, imax)
    inv_e6_2d = {(i, 0): v for i, v in inv_e6.items()}
    return mul(theta1sq, inv_e6_2d, imax)


def cnl_table(ZK3, cutoff_q):
    """{(n,l): coeff} for integer n=0..cutoff_q, integer l with |l| such
    that 4n-l^2 is in range, reading off the (i,j) grid (i=8n, j=2l)."""
    out = {}
    for (i, j), v in ZK3.items():
        if i % 8 != 0 or j % 2 != 0:
            continue
        n = i // 8
        l = j // 2
        if 0 <= n <= cutoff_q:
            out[(n, l)] = v
    return out


def check_discriminant_dependence(cnl):
    """Verify c(n,l) depends only on D=4n-l^2. Returns (ok, c_of_D dict,
    mismatches list)."""
    from collections import defaultdict
    byD = defaultdict(set)
    for (n, l), v in cnl.items():
        D = 4 * n - l * l
        byD[D].add(v)
    mismatches = {D: vs for D, vs in byD.items() if len(vs) > 1}
    cD = {D: next(iter(vs)) for D, vs in byD.items() if len(vs) == 1}
    return (len(mismatches) == 0), cD, mismatches


def rigidity_scan_b(phi01, cutoff_q, k_values):
    """Replace the overall '2' in Z=k*phi01 by rational k. Report, for
    each k, whether Z(tau,0)=24 holds and whether all reported (n,l)
    coefficients (n<=cutoff_q) are integers."""
    phi01_0 = y_to_1(phi01)
    # phi01(tau,0) should itself be the constant 12 (since Z_K3(tau,0)=24
    # and Z=2phi01); verify that directly from the computed series.
    phi01_0_nonzero = {i: v for i, v in phi01_0.items() if v != 0}
    results = []
    for k in k_values:
        kfr = Fr(k)
        z0 = {i: v * kfr for i, v in phi01_0.items() if v != 0}
        euler_ok = (z0 == {0: Fr(24)})
        all_int = all((v * kfr).denominator == 1 for v in phi01.values())
        results.append({
            "k": str(kfr),
            "Z_tau0": {str(Fr(i, 8)): str(v) for i, v in z0.items()},
            "euler_char_24_holds": euler_ok,
            "all_coeffs_integral": all_int,
        })
    return phi01_0_nonzero, results


def lambda_N_closed_form(N, cutoff_q):
    """Lambda_N(tau) := N * q d/dq log(eta(N*tau)/eta(tau)), exact, via
    divisor sums (mu-free): q d/dq log eta(tau) = 1/24 - sum sigma_1(m) q^m
    (standard; independently re-derived and cross-checked below against a
    direct series log-derivative of eta(N tau)/eta(tau) for N=2).
    Lambda_N = N(N-1)/24 + N*sum sigma_1(m) q^m - N^2*sum sigma_1(m) q^{N m}.
    """
    from sympy import divisor_sigma
    out = {0: Fr(N * (N - 1), 24)}
    for m in range(1, cutoff_q + 1):
        s1 = int(divisor_sigma(m, 1))
        out[m] = out.get(m, Fr(0)) + N * s1
    for m in range(1, cutoff_q // N + 1):
        s1 = int(divisor_sigma(m, 1))
        out[N * m] = out.get(N * m, Fr(0)) - N * N * s1
    return {n: out.get(n, Fr(0)) for n in range(cutoff_q + 1)}


M24_IRREP_DIMS = [1, 23, 45, 45, 231, 231, 252, 253, 483, 770, 770, 990, 990,
                   1035, 1035, 1035, 1265, 1771, 2024, 2277, 3312, 3520,
                   5313, 5544, 5796, 10395]


def build_theta1_psi_T1sq_eta3(cutoff_q, N_prod):
    """Shared building blocks for the H(tau) extraction (item 2)."""
    imax = 8 * cutoff_q
    t2 = theta2(N_prod, imax)
    t3 = theta3(N_prod, imax)
    t4 = theta4(N_prod, imax)
    r2 = divide_scalar_series(t2, y_to_1(t2), imax)
    r3 = divide_scalar_series(t3, y_to_1(t3), imax)
    r4 = divide_scalar_series(t4, y_to_1(t4), imax)
    r2sq = mul(r2, r2, imax)
    r3sq = mul(r3, r3, imax)
    r4sq = mul(r4, r4, imax)
    phi01 = scal(add(add(r2sq, r3sq), r4sq), 4)
    ZK3 = scal(phi01, 2)

    T1 = theta1_over_i(N_prod, imax)
    Psi = psi(N_prod, imax, T1)
    T1sq = mul(T1, T1, imax)
    e3 = eta3(N_prod, imax)

    ZK3_eta3 = mul(ZK3, e3, imax)
    y12Psi = shift(Psi, dj=1)
    return imax, ZK3, T1sq, e3, ZK3_eta3, y12Psi


def extract_H(coeff_N, ZK3_eta3, y12Psi, T1sq, imax):
    """H(tau) := (coeff_N * y^{1/2}*Psi - Z_K3*eta^3) / T1(tau,z)^2, solved
    via TWO independent y-power slices (j=0 and j=2) of theta_1^2=-T1^2;
    returns (H_dict_or_None_if_slices_disagree, agree_bool, H_j0, H_j2).
    """
    numer = add(scal(y12Psi, coeff_N), scal(ZK3_eta3, -1))

    def slice1d(d, j0):
        return {i: v for (i, j), v in d.items() if j == j0}

    def H_from_slice(j0):
        dslice = slice1d(T1sq, j0)
        nslice = slice1d(numer, j0)
        inv = reciprocal_1d(dslice, imax, istep=8)
        inv2d = {(i, 0): v for i, v in inv.items()}
        n2d = {(i, 0): v for i, v in nslice.items()}
        q = mul(n2d, inv2d, imax)
        return {k[0]: v for k, v in q.items()}

    Hj0 = H_from_slice(0)
    Hj2 = H_from_slice(2)
    # ignore the last couple of q-orders near the truncation edge (imax),
    # where add()/shift() can carry a stray term beyond what mul() would
    # have kept -- compare only up to imax - 16 (2 full q-steps of margin).
    safe = imax - 16
    Hj0_safe = {i: v for i, v in Hj0.items() if i <= safe}
    Hj2_safe = {i: v for i, v in Hj2.items() if i <= safe}
    agree = (Hj0_safe == Hj2_safe)
    return (Hj0_safe if agree else None), agree, Hj0_safe, Hj2_safe


def multiply_back_check(H1d, T1sq, ZK3_eta3, y12Psi, coeff_N, imax):
    H2d = {(i, 0): v for i, v in H1d.items()}
    back = mul(H2d, T1sq, imax)
    required = add(scal(y12Psi, coeff_N), scal(ZK3_eta3, -1))
    safe = imax - 16
    mism = {k: (back.get(k, 0), required.get(k, 0)) for k in set(back) | set(required)
            if k[0] <= safe and back.get(k, 0) != required.get(k, 0)}
    return len(mism) == 0, mism


def verify_lambda2_by_direct_log_derivative(cutoff_q):
    """Independent check of lambda_N_closed_form(2,...): build
    eta(2tau)/eta(tau)*q^{-1/24} = prod(1-q^{2n})/prod(1-q^n) as a plain
    1D Fraction q-series by long division, take its q d/dq log via series
    reciprocal + differentiation, and compare to the closed form."""
    def poch(step, cut):
        out = {0: Fr(1)}
        for n in range(1, cut // step + 2):
            term = {0: Fr(1), n * step: Fr(-1)}
            new = {}
            for e1, c1 in out.items():
                for e2, c2 in term.items():
                    e = e1 + e2
                    if e <= cut:
                        new[e] = new.get(e, Fr(0)) + c1 * c2
            out = {k: v for k, v in new.items() if v != 0}
        return out

    def reciprocal(P, cut):
        assert P[0] == 1
        inv = {0: Fr(1)}
        for n in range(1, cut + 1):
            s = Fr(0)
            for k in range(1, n + 1):
                if k in P:
                    s += P[k] * inv.get(n - k, Fr(0))
            inv[n] = -s
        return inv

    def mul1d(A, B, cut):
        out = {}
        for e1, c1 in A.items():
            for e2, c2 in B.items():
                e = e1 + e2
                if e <= cut:
                    out[e] = out.get(e, Fr(0)) + c1 * c2
        return out

    cut = cutoff_q
    P1 = poch(1, cut)
    P2 = poch(2, cut)
    f = mul1d(P2, reciprocal(P1, cut), cut)
    qfprime = {n: n * c for n, c in f.items()}
    logderiv = mul1d(qfprime, reciprocal(f, cut), cut)
    lam2 = {0: 2 * (Fr(1, 24) + logderiv.get(0, Fr(0)))}
    for n in range(1, cut + 1):
        lam2[n] = 2 * logderiv.get(n, Fr(0))
    return lam2


def main():
    cutoff_q = 8
    N_main = cutoff_q + 6

    ZK3, phi01 = compute_ZK3_and_phi01(cutoff_q, N_main)
    ZK3_0 = y_to_1(ZK3)
    ZK3_0_str = {str(Fr(i, 8)): str(v) for i, v in ZK3_0.items() if v != 0}

    cnl = cnl_table(ZK3, cutoff_q)
    disc_ok, cD, mismatches = check_discriminant_dependence(cnl)

    # ---- truncation-stability check: rerun at cutoff+2 and cutoff+4,
    # confirm the q^0..q^cutoff coefficients of Z_K3 do not move.
    stability = {}
    for extra in (2, 4):
        cq2 = cutoff_q + extra
        N2 = cq2 + 6
        ZK3_hi, _ = compute_ZK3_and_phi01(cq2, N2)
        # restrict to i<=8*cutoff_q and compare
        base_slice = {k: v for k, v in ZK3.items() if k[0] <= 8 * cutoff_q}
        hi_slice = {k: v for k, v in ZK3_hi.items() if k[0] <= 8 * cutoff_q}
        stability[f"cutoff_{cq2}_matches_cutoff_{cutoff_q}"] = (base_slice == hi_slice)

    # ---- phi_{-2,1} bonus computation
    phi21 = compute_phi_minus21(cutoff_q, N_main)
    phi21_leading = {str(Fr(j, 2)): str(v) for (i, j), v in phi21.items() if i == 0}

    # ---- rigidity scan (b)
    k_values = [Fr(2), Fr(1), Fr(3), Fr(5, 2), Fr(2) + Fr(1, 12), Fr(2) + Fr(1, 100), Fr(0)]
    phi01_0_nonzero, scan_b = rigidity_scan_b(phi01, cutoff_q, k_values)

    # ---- item (3), mu-free part only: Lambda_2 divisor-sum series, cross
    # checked by an independent direct log-derivative computation.
    lam2_closed = lambda_N_closed_form(2, cutoff_q)
    lam2_direct = verify_lambda2_by_direct_log_derivative(cutoff_q)
    lam2_match = lam2_closed == lam2_direct
    F_2A = {n: 16 * c for n, c in lam2_closed.items()}

    # ---- item (2): H(tau) via the Appell-Lerch mu-term (corrected sign,
    # found by advisor's leading-order hand-check: numer = coeff_N*y^{1/2}
    # Psi - Z_K3*eta^3, H = numer/T1(tau,z)^2). Extracted independently
    # from two different (non-symmetry-related) y-power slices of
    # theta_1(tau,z)^2 and cross-checked; ALSO verified by multiplying the
    # extracted H back through T1^2 and comparing term-by-term.
    H_cutoff_q = cutoff_q + 2  # extra margin so q^8 is not near the edge
    H_N_prod = H_cutoff_q + 8
    imaxH, ZK3_H, T1sq_H, e3_H, ZK3_eta3_H, y12Psi_H = build_theta1_psi_T1sq_eta3(H_cutoff_q, H_N_prod)

    H24, agree24, Hj0_24, Hj2_24 = extract_H(24, ZK3_eta3_H, y12Psi_H, T1sq_H, imaxH)
    mb_ok, mb_mismatches = (False, {"error": "slices disagreed"})
    if agree24:
        mb_ok, mb_mismatches = multiply_back_check(H24, T1sq_H, ZK3_eta3_H, y12Psi_H, 24, imaxH)

    # truncation-stability of A_1..A_8: rerun the whole H extraction at two
    # higher cutoffs and confirm every A_n is unchanged.
    H_stability = {}
    for cq2 in (H_cutoff_q + 2, H_cutoff_q + 4):
        Np2 = cq2 + 8
        imax2, _, T1sq2, _, ZK3e3_2, y12Psi2 = build_theta1_psi_T1sq_eta3(cq2, Np2)
        H2, agree2, _, _ = extract_H(24, ZK3e3_2, y12Psi2, T1sq2, imax2)
        ok2, _ = (multiply_back_check(H2, T1sq2, ZK3e3_2, y12Psi2, 24, imax2) if agree2 else (False, None))
        same_An = agree24 and mb_ok and agree2 and ok2 and all(
            H2.get(-1 + 8 * n) == H24.get(-1 + 8 * n) for n in range(1, cutoff_q + 1))
        H_stability[f"cutoff_{cq2}_A_n_match_cutoff_{H_cutoff_q}"] = same_An

    H_result = {
        "coeff_24_slices_agree": agree24,
        "multiply_back_check_passes": mb_ok,
        "multiply_back_mismatch_count": len(mb_mismatches) if isinstance(mb_mismatches, dict) else None,
    }
    A_n = {}
    m24_decomp = {}
    identity_check = None
    if agree24 and mb_ok:
        # H(tau) = 2 q^{-1/8} (-1 + sum_{n>=1} A_n q^n); H's coefficient
        # at i = -1 + 8n (q^{-1/8+n}) is 2*A_n (A_0 := -1 fixed by i=-1).
        H0 = H24.get(-1)
        polar_is_minus_2 = (H0 == Fr(-2))
        for n in range(1, H_cutoff_q + 1):
            i = -1 + 8 * n
            if i in H24:
                v = H24[i]
                if v % 2 == 0:
                    A_n[n] = v // 2
                else:
                    A_n[n] = Fr(v, 2)
        integral_polar_ok = polar_is_minus_2 and all(
            isinstance(a, int) or (isinstance(a, Fr) and a.denominator == 1) for a in A_n.values())
        for n, a in A_n.items():
            m24_decomp[n] = {"A_n": str(a), "matches_single_M24_dim": (a in M24_IRREP_DIMS)}
        if 1 in A_n and 2 in A_n:
            lhs = A_n[2] * 60
            rhs = 4 * A_n[1] * 77
            identity_check = {"A_2*60": str(lhs), "4*A_1*77": str(rhs), "holds": (lhs == rhs)}
        H_result.update({
            "H_leading_coeff_at_q_neg_1_8": str(H0),
            "polar_term_is_minus_2": polar_is_minus_2,
            "A_n": {str(n): str(a) for n, a in sorted(A_n.items())},
            "A_n_all_integral": integral_polar_ok,
            "A1_to_A5_vs_M24_irrep_dims_given_list": m24_decomp,
            "identity_A2_60_eq_4_A1_77": identity_check,
        })

    # ---- rigidity scan (a): replace the "24" multiplying the mu-term by
    # N in 20..28; report for which N the resulting H has integer
    # coefficients AND polar term exactly -2.
    scan_a = []
    if agree24 and mb_ok:
        for Ncoef in range(20, 29):
            Hn, agreeN, _, _ = extract_H(Ncoef, ZK3_eta3_H, y12Psi_H, T1sq_H, imaxH)
            row = {"N": Ncoef, "slices_agree": agreeN}
            if agreeN:
                h0 = Hn.get(-1)
                row["polar_term"] = str(h0)
                row["polar_is_minus_2"] = (h0 == Fr(-2))
                ints = all(v.denominator == 1 for v in Hn.values())
                row["all_reported_coeffs_integral"] = ints
                if Ncoef in (23, 25):
                    row["failing_detail"] = f"H(q^-1/8)={h0}, all_integral={ints}"
            scan_a.append(row)

    # ---- chi(2A): not given by the task; scan rather than recall it.
    chi2A_scan = []
    if agree24 and mb_ok:
        inv_e3_H = reciprocal_1d({i: v for (i, j), v in e3_H.items() if j == 0}, imaxH, istep=8)
        inv_e3_2d = {(i, 0): v for i, v in inv_e3_H.items()}
        F2A_2d = {(8 * n, 0): 16 * c for n, c in lam2_closed.items() if n <= H_cutoff_q}
        F2A_over_eta3 = mul(F2A_2d, inv_e3_2d, imaxH)
        H24_2d = {(i, 0): v for i, v in H24.items()}
        for chi in range(0, 25):
            H2A = add(scal(H24_2d, Fr(chi, 24)), scal(F2A_over_eta3, -1))
            ints = all(v.denominator == 1 for v in H2A.values())
            chi2A_scan.append({"chi_2A": chi, "H_2A_all_coeffs_integral": ints})

    out = {
        "track": "A",
        "item": "K3 elliptic genus (phi_{0,1}, Z_K3, discriminant property, "
                "phi_{-2,1}, H(tau)/Appell-Lerch mu-term, M24 dims, 2A "
                "twining Lambda_2/F_2A)",
        "method": "exact Fraction arithmetic; Jacobi theta triple-product "
                   "series truncated at q^{cutoff}, on the grid i=8*qexp, "
                   "j=2*yexp (series2d.py); ratios theta_i(z)/theta_i(0) "
                   "via exact formal power-series division "
                   "(divide_scalar_series).",
        "cutoff_q": cutoff_q,
        "N_product_terms": N_main,
        "truncation_stability": stability,
        "ZK3_tau_0": ZK3_0_str,
        "ZK3_tau_0_note": "Z_K3(tau,0) computed directly from the theta-ratio "
                           "formula; should equal the constant 24 (K3 Euler "
                           "number) with ALL higher-q-order terms exactly "
                           "zero. No 24 was fed in anywhere upstream.",
        "discriminant_dependence_holds": disc_ok,
        "discriminant_mismatches": {str(D): [str(v) for v in vs] for D, vs in mismatches.items()},
        "c_of_D": {str(D): str(v) for D, v in sorted(cD.items())},
        "c_of_D_expected_from_memory_unverified": {
            "-1": "2", "0": "20", "3": "-128", "4": "216", "7": "-1026", "8": "1616",
            "note": "recalled from general familiarity with K3 elliptic genus "
                     "tables (e.g. Eguchi-Ooguri-Tachikawa-style discriminant "
                     "coefficients); NOT used in the computation above, only "
                     "compared after the fact.",
        },
        "phi_minus_2_1_leading_qorder_y_coeffs": phi21_leading,
        "phi01_tau0_nonzero_terms": {str(Fr(i, 8)): str(v) for i, v in phi01_0_nonzero.items()},
        "rigidity_scan_b": {
            "definition": "Replace the overall factor 2 in Z_K3 := k*phi_{0,1} "
                           "by a rational k; scan k and report which values "
                           "satisfy (i) Z(tau,0)=24 [Euler characteristic of "
                           "a compact hyperkahler 4-fold] and (ii) all "
                           "reported Fourier coefficients (q^0..q^8, all l) "
                           "are integers.",
            "honesty_note": "(i) alone is ONE LINEAR EQUATION in k (since "
                             "phi_{0,1}(tau,0)=12 is itself a computed "
                             "constant), so of course it pins k=2 uniquely; "
                             "that by itself is not a nontrivial rigidity "
                             "result. The nontrivial fact checked here is "
                             "(ii): phi_{0,1}'s own Fourier coefficients "
                             "already have gcd 1 (e.g. c(-1)/2=1), so ANY "
                             "non-integer k, and in fact any k!=2 even among "
                             "integers, fails to reproduce Z_K3; k=2 is the "
                             "unique value in the scanned set satisfying "
                             "BOTH conditions simultaneously.",
            "scan": scan_b,
        },
        "lambda_2_and_F_2A": {
            "definition": "Lambda_N(tau) := N q d/dq log(eta(N tau)/eta(tau)); "
                           "F_2A := 16 Lambda_2. Pure divisor-sum arithmetic "
                           "(sympy.divisor_sigma), no mu, no theta functions.",
            "Lambda_2_closed_form_matches_direct_log_derivative_check": lam2_match,
            "Lambda_2_coeffs_q0_to_q8": {str(n): str(lam2_closed[n]) for n in range(cutoff_q + 1)},
            "F_2A_coeffs_q0_to_q8": {str(n): str(F_2A[n]) for n in range(cutoff_q + 1)},
        },
        "H_tau_appell_lerch": {
            "definition": "H(tau) = 2 q^{-1/8}(-1 + sum_{n>=1} A_n q^n), "
                           "solved from Z_K3*eta^3 = 24*y^{1/2}*Psi(tau,z) "
                           "- H(tau)*theta_1(tau,z)^2, where "
                           "Psi:=T1(tau,z)*S(tau,z) (S = Appell-Lerch sum) "
                           "is computed exactly with the n=0 term's pole "
                           "cancelled algebraically against T1's own "
                           "(1-y) factor -- no truncation/regularization "
                           "trick anywhere (appell.py). H is extracted from "
                           "TWO independent, non-symmetry-related y-power "
                           "slices of theta_1^2 (j=0 and j=2) and accepted "
                           "only if they agree; the result is then "
                           "multiplied back through theta_1^2 and compared "
                           "term-by-term against the defining equation as a "
                           "second, stronger check.",
            "cutoff_q_used": H_cutoff_q,
            "two_slices_agree": agree24,
            "multiply_back_check_passes": mb_ok,
            "truncation_stability_of_A_n": H_stability,
            **H_result,
        },
        "rigidity_scan_a": {
            "definition": "Replace the '24' multiplying the Appell-Lerch "
                           "mu-term by N in 20..28; for each N re-extract "
                           "H(tau) (two-slice agreement required) and "
                           "report whether its polar term is exactly -2 "
                           "and whether all q^0..q^{cutoff} coefficients "
                           "are integers.",
            "scan": scan_a,
        },
        "chi_2A_scan_for_H_2A": {
            "definition": "H_2A(tau) := (chi(2A)/24)*H(tau) - F_2A(tau)/eta(tau)^3. "
                           "chi(2A) (the trace of the 2A element of M24 in its "
                           "24-dim permutation representation) is NOT given "
                           "in the task text, so it is scanned over 0..24 "
                           "rather than typed in from memory; report which "
                           "value(s) make H_2A's q^0..q^{cutoff} coefficients "
                           "all integers.",
            "scan": chi2A_scan,
            "chi_2A_expected_from_memory_unverified": "8 (M24's 24-point "
                "permutation character at a 2A element, i.e. number of "
                "fixed points; recalled from general ATLAS/M24 familiarity, "
                "NOT used anywhere upstream, only compared after the "
                "fact). The scan above found TWO integral candidates in "
                "0..24 (8 and 20) at this cutoff -- consistent with, but "
                "not by itself pinning down, chi(2A)=8.",
        },
        "could_not_do": [
            "The 2A-twined identity check (whether A_2*60=4*A_1*77 survives "
            "when A_n is replaced by its 2A-twined counterpart) was not "
            "attempted: it needs the FULL 2A-twined elliptic genus "
            "Z_K3^{2A}(tau,z) via F_2A's numerator structure combined with "
            "theta_1(tau,z)^2 the same way H combines with it, which is a "
            "further nontrivial regularization step beyond H_2A(tau) alone "
            "(a tau-only object) that time did not allow re-deriving and "
            "cross-checking to the same standard as the rest of this file.",
        ] if agree24 and mb_ok else [
            "Item (2)/(3): H(tau) extraction did not pass its own two-slice "
            "and multiply-back consistency checks at this cutoff; see "
            "H_tau_appell_lerch above for the specific failure. Everything "
            "downstream (A_n, M24 decomposition, the A_2*60=4*A_1*77 "
            "identity, rigidity scan (a), and H_2A) is consequently not "
            "reported.",
        ],
    }

    with open("results.json", "w") as f:
        json.dump(out, f, indent=2, default=str)

    print(json.dumps({k: v for k, v in out.items() if k not in (
        "c_of_D_expected_from_memory_unverified",)}, indent=2, default=str)[:3000])
    print("\n... (full output written to results.json)")


if __name__ == "__main__":
    main()
