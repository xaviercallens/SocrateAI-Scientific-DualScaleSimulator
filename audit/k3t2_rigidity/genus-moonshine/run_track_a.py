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

from series2d import y_to_1, divide_scalar_series, mul, add, scal
from thetas import theta2, theta3, theta4, theta1_over_i, eta6, prod1


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

    out = {
        "track": "A",
        "item": "K3 elliptic genus (phi_{0,1}, Z_K3, discriminant property, phi_{-2,1})",
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
        "could_not_do": [
            "Item (2): H(tau) extraction via the Appell-Lerch mu-term, the "
            "M24 decomposition of A_1..A_8, and the A_2*60=4*A_1*77 "
            "identity check. What WAS done exactly and finitely (appell.py): "
            "Psi(tau,z) := T1(tau,z)*S(tau,z), where S is the Appell-Lerch "
            "sum, with the n=0 term's naive 1/(1-y) pole cancelled "
            "ALGEBRAICALLY against T1's own exact factor of (1-y) -- no "
            "infinite-tail truncation or regularization prescription "
            "anywhere. A candidate normalization with self-consistent q,y "
            "gradings was found (Z_K3*eta^3 = -24*y^{1/2}*Psi + H*theta_1^2, "
            "i.e. theta_1^2/eta^3 multiplies the WHOLE bracket [24*mu+H], "
            "not H alone) and its leading term, read off one y-power slice "
            "of theta_1^2, landed on H_0=-2 -- matching the task's stated "
            "polar term. But the DEFINITIVE check (multiply the candidate "
            "H back through theta_1^2 and compare term-by-term against "
            "-[Z_K3*eta^3+24*y^{1/2}*Psi], not just one y-slice) shows real "
            "mismatches starting at the very leading q-order (e.g. at "
            "q^{1/8}*y^0: 4 from the candidate vs -44 required). This is a "
            "genuine normalization/formula error, not a coding bug that "
            "was left unfixed: two DIFFERENT y-power slices (j=0 vs j=2) "
            "of theta_1^2, which should give the identical y-independent "
            "H(tau) if the formula were right, disagree outright (+22 vs "
            "-2 at the leading order; j=2 and j=-2 agreeing with each "
            "other is a trivial consequence of theta_1^2's y<->1/y "
            "symmetry, not an independent check). No H(tau)/A_n sequence "
            "is reported because none was found that passes this check.",
            "Item (3), the H-dependent half: the 2A-twined series H_2A and "
            "the survival of the A_2*60=4*A_1*77 identity under twining, "
            "and rigidity scan (a) on the '24' in the mu-term -- all "
            "depend on the item-(2) extraction above and were not "
            "attempted once that failed its own consistency check. The "
            "mu-FREE half of item (3), Lambda_2 (needed for F_2A = "
            "16*Lambda_2), WAS computed exactly and cross-checked by an "
            "independent direct log-derivative computation (see "
            "lambda_2_and_F_2A above) -- it just cannot be turned into "
            "H_2A without H.",
        ],
    }

    with open("results.json", "w") as f:
        json.dump(out, f, indent=2, default=str)

    print(json.dumps({k: v for k, v in out.items() if k not in (
        "c_of_D_expected_from_memory_unverified",)}, indent=2, default=str)[:3000])
    print("\n... (full output written to results.json)")


if __name__ == "__main__":
    main()
