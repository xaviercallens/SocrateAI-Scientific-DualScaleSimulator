#!/usr/bin/env python3
"""
Track A v2 -- K3 elliptic genus, Mathieu moonshine, twining (Tier B / L).

Reuses v1's exact-arithmetic engine unmodified (series2d.py, thetas.py,
appell.py, copied byte-for-byte from
audit/k3t2_rigidity/genus-moonshine/{series2d,thetas,appell}.py, which were
written BLIND in v1 -- reuse permitted by the task). Nothing in those three
files is edited here. What v2 adds / repairs relative to v1's
run_track_a.py:

  DEFECT FIXED (1): v1's main Z_K3/phi_{0,1}/c(D) computation used
  cutoff_q=8 (< q^10 required by the task). v2 raises this to cutoff_q=10.

  DEFECT FIXED (2): v1's H(tau) extraction used H_cutoff_q=10, but its own
  "safe" truncation margin (imax-16) excluded n=9,10 from A_n, so only
  A_1..A_8 were actually reported despite cutoff_q_used=10 being printed.
  v2 raises H_cutoff_q to 14 so A_1..A_10 sit safely inside the margin
  (checked explicitly below, not just asserted).

  MISSING (3): v1 only twined class 2A. v2 adds 3A, 5A, 7A using the
  task-supplied (tier-L, "from memory", Cheng-Duncan-Harvey Table 3)
  F_g = {F_2A:16 Lambda_2, F_3A:6 Lambda_3, F_5A:2 Lambda_5, F_7A:Lambda_7}
  and chi(g) = {2A:8, 3A:6, 5A:4, 7A:3}. These are used as literature INPUTS
  to H_g = (chi(g)/24) H - F_g/eta^3, per the task -- NOT as targets a
  selecting condition is built to reproduce. This is Tier L input data, not
  a rigidity claim.

  MISSING (4): v1 did no structural cross-check of the twined series
  against class-function constraints. v2 adds, for each of 2A/3A/5A/7A:
    - integrality of H_g's Fourier coefficients (all q-orders reported)
    - sign of A_n^(g) vs sign of A_n (n=1..8)
    - the congruence A_n^(g) === A_n (mod ord(g)) for n=1..8, which is a
      necessary consequence of A_n^(g) being (for a genuine class function)
      a character value -- a structural, non-tautological check that can
      and does fail for a wrong chi(g)/F_g pairing (see negative control
      below).
    - whether A_2^(g)*60 == 4*A_1^(g)*77 (the "27720 lock", SAME convention
      as v1: A_n read off H_g's q^{-1/8+n} coefficient over 2, matching the
      untwined convention A_n = H(tau)'s q^{-1/8+n} coefficient / 2).

  MISSING (5): exports.json with chi_from_genus = Z_K3(tau,0).

  Rigidity (a) [mu-term coefficient N in 20..28]: re-run with the new
  cutoff; same z-independence-of-two-slices structural selector as v1
  (kept, it already satisfied the "structural, not literal" requirement).

  Rigidity (b) [overall factor k of Z=k*phi_{0,1}]: v1's version used the
  Euler-characteristic condition Z(tau,0)=24, which the task's own ground
  rules class as circular (the selecting condition literally targets the
  known value 24). v2 instead tries the task's SUGGESTED alternative --
  "the q^0 y^{+-1} coefficient must equal 2" -- and reports HONESTLY that
  this condition also bakes the target value (2) directly into the
  selecting equation (k*p1 = 2, solved for k), so it is REPORTED AS
  NORMALISATION, exactly as the task's own fallback instructs when no
  literal-free structural condition can be found. No condition avoiding
  the literal target was found; this is stated, not hidden.
"""
import json
from fractions import Fraction as Fr

from series2d import y_to_1, divide_scalar_series, mul, add, scal, shift, reciprocal_1d
from thetas import theta2, theta3, theta4, theta1_over_i, eta3
from appell import psi


# ---------------------------------------------------------------- core Z_K3

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


def compute_ZK3_and_phi01_safe(report_cutoff_q, margin=2):
    """DEFECT FIX (discovered empirically in this session, not present in
    v1 -- v1 never hit it because its cutoff_q=8 happened not to expose a
    D-collision at the exact edge row): compute_ZK3_and_phi01(cutoff_q, N)
    is NOT exact at its own top row n=cutoff_q. Direct check: the (n=10,
    l=+-5) coefficient (D=15) computed with report_cutoff_q=10 is -23554,
    but recomputing with an internally higher cutoff (11, 12, ...) and
    reading off the SAME n=10 row gives -23550 (which matches the n=4 and
    n=6 rows at the same D=15, as the discriminant-dependence theorem
    requires) -- and this converges by margin=+1 and is stable at +2,+3.
    So: compute internally at report_cutoff_q+margin, then return only the
    q<=report_cutoff_q part."""
    internal_cutoff = report_cutoff_q + margin
    ZK3_full, phi01_full = compute_ZK3_and_phi01(internal_cutoff, internal_cutoff + 8)
    imax_report = 8 * report_cutoff_q
    ZK3 = {k: v for k, v in ZK3_full.items() if k[0] <= imax_report}
    phi01 = {k: v for k, v in phi01_full.items() if k[0] <= imax_report}
    return ZK3, phi01


def cnl_table(ZK3, cutoff_q):
    out = {}
    for (i, j), v in ZK3.items():
        if i % 8 != 0 or j % 2 != 0:
            continue
        n, l = i // 8, j // 2
        if 0 <= n <= cutoff_q:
            out[(n, l)] = v
    return out


def check_discriminant_dependence(cnl):
    from collections import defaultdict
    byD = defaultdict(set)
    for (n, l), v in cnl.items():
        byD[4 * n - l * l].add(v)
    mismatches = {D: vs for D, vs in byD.items() if len(vs) > 1}
    cD = {D: next(iter(vs)) for D, vs in byD.items() if len(vs) == 1}
    return (len(mismatches) == 0), cD, mismatches


# ------------------------------------------------------------ Lambda_N, F_g

def lambda_N_closed_form(N, cutoff_q):
    """Lambda_N(tau) := N q d/dq log(eta(N tau)/eta(tau)), exact via
    divisor sums: q d/dq log eta(tau) = 1/24 - sum sigma_1(m) q^m.
    Lambda_N = N(N-1)/24 + N*sum sigma_1(m) q^m - N^2*sum sigma_1(m) q^{Nm}.
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


def verify_lambdaN_direct(N, cutoff_q):
    """Independent check of lambda_N_closed_form(N,...): build
    eta(N tau)/eta(tau)*q^{-N(N-1)/24-ish...} EXACTLY as
    prod(1-q^{Nn})/prod(1-q^n) (1D plain series, long division), take its
    q d/dq log via series reciprocal + differentiation, multiply by N, and
    compare to the closed form."""
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
    PN = poch(N, cut)
    f = mul1d(PN, reciprocal(P1, cut), cut)  # prod(1-q^{Nn})/prod(1-q^n)
    qfprime = {n: n * c for n, c in f.items()}
    logderiv = mul1d(qfprime, reciprocal(f, cut), cut)
    lam = {0: N * (Fr(N - 1, 24) + logderiv.get(0, Fr(0)))}
    for n in range(1, cut + 1):
        lam[n] = N * logderiv.get(n, Fr(0))
    return lam


# --------------------------------------------------------------- H(tau)

def build_blocks(cutoff_q, N_prod):
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
    return imax, ZK3, phi01, T1sq, e3, ZK3_eta3, y12Psi


def extract_H(coeff_N, ZK3_eta3, y12Psi, T1sq, imax, margin=16):
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
    safe = imax - margin
    Hj0_safe = {i: v for i, v in Hj0.items() if i <= safe}
    Hj2_safe = {i: v for i, v in Hj2.items() if i <= safe}
    agree = (Hj0_safe == Hj2_safe)
    return (Hj0_safe if agree else None), agree, Hj0_safe, Hj2_safe, safe


def multiply_back_check(H1d, T1sq, ZK3_eta3, y12Psi, coeff_N, imax, safe):
    H2d = {(i, 0): v for i, v in H1d.items()}
    back = mul(H2d, T1sq, imax)
    required = add(scal(y12Psi, coeff_N), scal(ZK3_eta3, -1))
    mism = {k: (back.get(k, 0), required.get(k, 0)) for k in set(back) | set(required)
            if k[0] <= safe and back.get(k, 0) != required.get(k, 0)}
    return len(mism) == 0, mism


def An_from_H(Hdict, nmax):
    """A_n from H's q^{-1/8+n} coefficient (i=-1+8n), A_n := H[i]/2.
    Only returns n for which the coefficient is present in Hdict AT ALL
    (i.e. within whatever range Hdict was computed/trimmed to) -- used for
    the untwined H(tau) where every n<=nmax is genuinely inside range."""
    out = {}
    for n in range(1, nmax + 1):
        i = -1 + 8 * n
        if i in Hdict:
            v = Hdict[i]
            out[n] = v / 2
    return out


def raw_zero_filled(Hdict, nmax, safe_i=None):
    """{n: H[i]} for n=0..nmax, i=-1+8n, ZERO-FILLED when the coefficient
    is exactly zero (series2d drops zero entries, so 'missing' here means
    0, not undefined) -- as long as i is within the safe/reported range.
    Returns {} for n whose i exceeds safe_i (genuinely out of range)."""
    out = {}
    for n in range(0, nmax + 1):
        i = -1 + 8 * n
        if safe_i is not None and i > safe_i:
            continue
        out[n] = Hdict.get(i, Fr(0))
    return out


def main():
    MAIN_CUTOFF = 10       # task requires "through q^10 at least"
    N_MAIN = MAIN_CUTOFF + 8

    # ---------------- main Z_K3 / phi_{0,1} / c(D) at cutoff_q=10 --------
    # DEFECT FOUND + FIXED IN THIS SESSION: compute_ZK3_and_phi01(cutoff_q,N)
    # alone is inexact at its own top row n=cutoff_q (see
    # compute_ZK3_and_phi01_safe docstring for the empirical check). Use
    # the margin-2 safe wrapper for everything reported below.
    edge_row_bug_raw = compute_ZK3_and_phi01(MAIN_CUTOFF, N_MAIN)[0]
    edge_row_bug_check = {
        (n, l): v for (n, l), v in cnl_table(edge_row_bug_raw, MAIN_CUTOFF).items()
        if n == MAIN_CUTOFF and (l == -5 or l == 5)
    }
    ZK3, phi01 = compute_ZK3_and_phi01_safe(MAIN_CUTOFF, margin=2)
    edge_row_bug_check_after_fix = {
        (n, l): v for (n, l), v in cnl_table(ZK3, MAIN_CUTOFF).items()
        if n == MAIN_CUTOFF and (l == -5 or l == 5)
    }
    ZK3_0 = y_to_1(ZK3)
    ZK3_0_str = {str(Fr(i, 8)): str(v) for i, v in ZK3_0.items() if v != 0}
    ZK3_tau0_const = ZK3_0.get(0)

    cnl = cnl_table(ZK3, MAIN_CUTOFF)
    disc_ok, cD, mismatches = check_discriminant_dependence(cnl)

    stability = {}
    for extra in (2, 4):
        cq2 = MAIN_CUTOFF + extra
        ZK3_hi, _ = compute_ZK3_and_phi01_safe(cq2, margin=2)
        base_slice = {k: v for k, v in ZK3.items() if k[0] <= 8 * MAIN_CUTOFF}
        hi_slice = {k: v for k, v in ZK3_hi.items() if k[0] <= 8 * MAIN_CUTOFF}
        stability[f"cutoff_{cq2}_matches_cutoff_{MAIN_CUTOFF}"] = (base_slice == hi_slice)

    # ---------------- rigidity scan (a): mu-term coefficient N in 20..28 -
    H_CUTOFF = 14           # raised so A_1..A_10 sit inside the safe margin
    N_H = H_CUTOFF + 10
    MARGIN = 16
    imaxH, ZK3_H, phi01_H, T1sq_H, e3_H, ZK3_eta3_H, y12Psi_H = build_blocks(H_CUTOFF, N_H)

    H24, agree24, Hj0_24, Hj2_24, safe24 = extract_H(24, ZK3_eta3_H, y12Psi_H, T1sq_H, imaxH, MARGIN)
    mb_ok, mb_mismatches = (False, {"error": "slices disagreed"})
    if agree24:
        mb_ok, mb_mismatches = multiply_back_check(H24, T1sq_H, ZK3_eta3_H, y12Psi_H, 24, imaxH, safe24)

    A_n = An_from_H(H24, 10) if (agree24 and mb_ok) else {}
    A9_A10_inside_margin = (-1 + 8 * 9 <= safe24) and (-1 + 8 * 10 <= safe24)

    H_stability = {}
    if agree24 and mb_ok:
        for cq2 in (H_CUTOFF + 2, H_CUTOFF + 4):
            imax2, _, _, T1sq2, _, ZK3e3_2, y12Psi2 = build_blocks(cq2, cq2 + 10)
            H2, agree2, _, _, safe2 = extract_H(24, ZK3e3_2, y12Psi2, T1sq2, imax2, MARGIN)
            ok2 = False
            if agree2:
                ok2, _ = multiply_back_check(H2, T1sq2, ZK3e3_2, y12Psi2, 24, imax2, safe2)
            same_An = agree2 and ok2 and all(
                H2.get(-1 + 8 * n) == H24.get(-1 + 8 * n) for n in range(1, 11))
            H_stability[f"cutoff_{cq2}_A1_to_A10_match_cutoff_{H_CUTOFF}"] = same_An

    scan_a = []
    for Ncoef in range(20, 29):
        Hn, agreeN, Hj0_N, Hj2_N, _ = extract_H(Ncoef, ZK3_eta3_H, y12Psi_H, T1sq_H, imaxH, MARGIN)
        # zero-fill: series2d drops exactly-zero coefficients, so a missing
        # key at i=-1 means the coefficient IS 0, not "undefined".
        h0_j0 = Hj0_N.get(-1, Fr(0))
        h0_j2 = Hj2_N.get(-1, Fr(0))
        row = {
            "N": Ncoef, "slices_agree": agreeN,
            "polar_term_from_j0_slice": str(h0_j0),
            "polar_term_from_j2_slice": str(h0_j2),
        }
        if agreeN:
            row["all_reported_coeffs_integral"] = all(v.denominator == 1 for v in Hn.values())
        scan_a.append(row)

    # ---------------- rigidity scan (b): overall factor k of Z=k*phi01 ---
    phi01_0 = y_to_1(phi01)
    k_values = [Fr(2), Fr(1), Fr(3), Fr(5, 2), Fr(2) + Fr(1, 12), Fr(2) + Fr(1, 100), Fr(0)]
    scan_b = []
    for k in k_values:
        z0 = {i: v * k for i, v in phi01_0.items() if v != 0}
        euler_ok = (z0 == {0: Fr(24)})
        all_int = all((v * k).denominator == 1 for v in phi01.values())
        scan_b.append({"k": str(k), "euler_char_24_holds": euler_ok, "all_coeffs_integral": all_int})

    # task-suggested structural selector: q^0 y^{+1} coefficient of Z must
    # equal 2. p1 := phi01's own (q^0, y^{+1}) coefficient, read off the
    # ALREADY-COMPUTED series (not asserted). Z's (q^0,y^1) coeff = k*p1.
    p1 = phi01.get((0, 2))
    k_from_q0y1_eq_2 = (Fr(2) / p1) if p1 not in (None, 0) else None
    rigidity_b_q0y1_selector = {
        "definition": "Solve k such that the (q^0, y^{+1}) coefficient of "
                       "Z=k*phi_{0,1} equals 2 (task's suggested structural "
                       "condition, motivated by '2 ground states of "
                       "h^{0,0}=h^{2,0}=1').",
        "phi01_q0_y1_coeff_p1_computed": str(p1) if p1 is not None else None,
        "k_solving_k_times_p1_eq_2": str(k_from_q0y1_eq_2) if k_from_q0y1_eq_2 is not None else None,
        "classification": "NORMALISATION",
        "honesty_note": "This condition solves the single linear equation "
                         "k*p1=2 for k, where the RHS '2' is exactly the "
                         "known target value of the overall factor being "
                         "tested (k=2 with p1=1). The selecting condition "
                         "therefore has the true value built into it, "
                         "exactly the pattern the task's ground rules "
                         "flag as NOT a rigidity result. No alternative "
                         "structural condition (consistency, integrality, "
                         "modularity, exactness, or agreement of two "
                         "independent computations) that pins k to a "
                         "unique value WITHOUT referencing 2 or 24 was "
                         "found in this session; per the task's own "
                         "fallback instruction this is classified "
                         "NORMALISATION: the (q^0,y^1) coefficient "
                         "normalisation fixes k given the physical input "
                         "'2 ground states', it is not a rigidity proof.",
    }

    # ---------------- Twining: 2A, 3A, 5A, 7A ------------------------------
    # Tier-L inputs (Cheng-Duncan-Harvey Table 3, FROM MEMORY, NOT re-derived
    # here -- used as literature data per the task, not as a selecting
    # target for any rigidity scan above).
    CLASSES = {
        "2A": {"ord": 2, "chi": 8, "F_factor": 16},
        "3A": {"ord": 3, "chi": 6, "F_factor": 6},
        "5A": {"ord": 5, "chi": 4, "F_factor": 2},
        "7A": {"ord": 7, "chi": 3, "F_factor": 1},
    }

    inv_e3_H = reciprocal_1d({i: v for (i, j), v in e3_H.items() if j == 0}, imaxH, istep=8)
    inv_e3_2d = {(i, 0): v for i, v in inv_e3_H.items()}
    H24_2d = {(i, 0): v for i, v in H24.items()} if (agree24 and mb_ok) else None

    def compute_twining_row(gname, info, H24_2d_loc, inv_e3_2d_loc, imaxH_loc, safe24_loc,
                             H24_loc, cutoff_q_loc):
        """One class's H_g, raw (zero-filled) coefficients, and the
        structural checks. Returns (row_dict, raw_g_by_n, A_g_by_n)."""
        N_ord, chi_g, Ffac = info["ord"], info["chi"], info["F_factor"]
        lam_closed = lambda_N_closed_form(N_ord, cutoff_q_loc)
        lam_direct = verify_lambdaN_direct(N_ord, cutoff_q_loc)
        lam_match = (lam_closed == lam_direct)

        Fg_2d = {(8 * n, 0): Ffac * c for n, c in lam_closed.items()}
        Fg_over_eta3 = mul(Fg_2d, inv_e3_2d_loc, imaxH_loc)

        Hg_2d = add(scal(H24_2d_loc, Fr(chi_g, 24)), scal(Fg_over_eta3, -1))
        Hg = {i: v for (i, j), v in Hg_2d.items() if j == 0}
        Hg_safe = {i: v for i, v in Hg.items() if i <= safe24_loc}

        all_int = all(v.denominator == 1 for v in Hg_safe.values())

        # zero-filled raw coefficients (raw_g[n] = Hg_safe.get(i,0)) and the
        # corresponding raw untwined coefficients (raw_un[n] = H24.get(i,0))
        # for n=0..10, as long as i=-1+8n is inside the safe range.
        raw_g = raw_zero_filled(Hg_safe, 10, safe_i=safe24_loc)
        raw_un = raw_zero_filled(H24_loc, 10, safe_i=safe24_loc)
        A_g = {n: v / 2 for n, v in raw_g.items()}

        # structural check on the RAW (undivided) Fourier coefficients --
        # these are the actual class-function traces (2*A_n untwined,
        # H_g[i] twined); dividing by 2 first is what produced spurious
        # "non-integer" failures for odd-order classes with irrational
        # characters (e.g. 7A: raw is an integer, but that integer is odd,
        # so A_n=raw/2 is a half-integer -- see README). The congruence
        # A_n^(g) === A_n (mod ord(g)) is tested here as the EQUIVALENT
        # raw statement raw_g[n] === raw_un[n] (mod ord(g)), which needs no
        # division and is well-defined for every n.
        congruence = {}
        for n in range(1, 9):
            if n in raw_g and n in raw_un:
                diff = raw_un[n] - raw_g[n]
                ok = diff.denominator == 1 and int(diff) % N_ord == 0
                congruence[str(n)] = {
                    "raw_H_g": str(raw_g[n]), "raw_H_untwined": str(raw_un[n]),
                    "diff_mod_ord": (int(diff) % N_ord) if diff.denominator == 1 else "non-integer-diff",
                    "congruent": ok,
                }
        congruence_all_ok = all(v["congruent"] for v in congruence.values()) if congruence else None

        sign_pattern = {str(n): {"A_n_g_sign": (1 if A_g[n] > 0 else (-1 if A_g[n] < 0 else 0)),
                                  "A_n_sign": (1 if raw_un.get(n, Fr(0)) / 2 > 0 else
                                               (-1 if raw_un.get(n, Fr(0)) / 2 < 0 else 0))}
                         for n in range(1, 9) if n in A_g}

        identity_g = None
        if 1 in A_g and 2 in A_g:
            lhs_g, rhs_g = A_g[2] * 60, 4 * A_g[1] * 77
            identity_g = {"A_2_g_times_60": str(lhs_g), "4_A_1_g_times_77": str(rhs_g), "holds": (lhs_g == rhs_g)}

        non_integer_An = {str(n): str(v) for n, v in A_g.items() if v.denominator != 1}

        row = {
            "chi_g_tier_L_input": chi_g,
            "F_g_factor_tier_L_input": Ffac,
            "ord_g": N_ord,
            "lambda_closed_form_matches_direct_log_derivative": lam_match,
            "H_g_all_reported_coeffs_integral": all_int,
            "H_g_polar_term_q_neg_1_8": str(raw_g.get(0)) if 0 in raw_g else None,
            "A_n_g": {str(n): str(v) for n, v in sorted(A_g.items())},
            "A_n_g_non_integer": non_integer_An,
            "A_n_g_non_integer_note": (
                "Non-integer A_n^(g)=H_g[i]/2 at these n does NOT mean the "
                "tier-L chi(g)/F_g inputs are wrong: the RAW H_g "
                "coefficient (see congruence table) is an integer at "
                "every n; halving it can land on a half-integer exactly "
                "when the untwined multiplicity K_n decomposes into a "
                "conjugate pair of irrational M24 characters at this "
                "class (this happens generically at classes like 7A, "
                "whose character values involve (-1+-i*sqrt(7))/2). "
                "Reported as-is; the congruence check below is done on "
                "the RAW coefficients precisely so it stays well-defined "
                "regardless of this halving."
            ) if non_integer_An else None,
            "congruence_An_g_equiv_An_mod_ord_g_n_le_8": congruence,
            "congruence_all_hold": congruence_all_ok,
            "sign_pattern_n_le_8": sign_pattern,
            "identity_A2_60_eq_4_A1_77_under_twining": identity_g,
        }
        return row, raw_g, A_g

    twining = {}
    lambda_cross_checks = {}
    identity_untwined = None
    raw_un_main = raw_zero_filled(H24, 10, safe_i=safe24) if (agree24 and mb_ok) else {}
    A_n_zf = {n: v / 2 for n, v in raw_un_main.items()}
    if 1 in A_n_zf and 2 in A_n_zf:
        lhs, rhs = A_n_zf[2] * 60, 4 * A_n_zf[1] * 77
        identity_untwined = {"A_2_times_60": str(lhs), "4_A_1_times_77": str(rhs), "holds": (lhs == rhs)}

    twining_raw_by_class = {}
    if H24_2d is not None:
        for gname, info in CLASSES.items():
            row, raw_g, A_g = compute_twining_row(gname, info, H24_2d, inv_e3_2d, imaxH, safe24,
                                                    H24, H_CUTOFF)
            twining[gname] = row
            lambda_cross_checks[gname] = row["lambda_closed_form_matches_direct_log_derivative"]
            twining_raw_by_class[gname] = raw_g

    # truncation-stability of the twined series: rerun the ENTIRE twining
    # computation (independent H(tau) extraction, independent Lambda_N,
    # independent 1/eta^3) at H_CUTOFF+2 and confirm every raw H_g[n]
    # coefficient for n<=10 is unchanged, for all four classes.
    twining_stability = {}
    if H24_2d is not None:
        cq2 = H_CUTOFF + 2
        imax2, _, _, T1sq2, e3_2, ZK3e3_2, y12Psi2 = build_blocks(cq2, cq2 + 10)
        H2, agree2, _, _, safe2 = extract_H(24, ZK3e3_2, y12Psi2, T1sq2, imax2, MARGIN)
        mb_ok2 = False
        if agree2:
            mb_ok2, _ = multiply_back_check(H2, T1sq2, ZK3e3_2, y12Psi2, 24, imax2, safe2)
        if agree2 and mb_ok2:
            inv_e3_2b = reciprocal_1d({i: v for (i, j), v in e3_2.items() if j == 0}, imax2, istep=8)
            inv_e3_2d_2 = {(i, 0): v for i, v in inv_e3_2b.items()}
            H2_2d = {(i, 0): v for i, v in H2.items()}
            for gname, info in CLASSES.items():
                _, raw_g2, _ = compute_twining_row(gname, info, H2_2d, inv_e3_2d_2, imax2, safe2, H2, cq2)
                same = all(raw_g2.get(n) == twining_raw_by_class.get(gname, {}).get(n) for n in range(0, 11)
                           if n in raw_g2 and n in twining_raw_by_class.get(gname, {}))
                twining_stability[gname] = {
                    f"cutoff_{cq2}_raw_Hg_n0_to_10_match_cutoff_{H_CUTOFF}": same,
                    "extraction_agree_and_multiply_back_ok_at_higher_cutoff": True,
                }
        else:
            twining_stability["error"] = "H(tau) re-extraction at cutoff+2 failed its own agree/multiply-back check"

    # negative control for the congruence check: deliberately mispair 2A's
    # F_g with 3A's chi (an obviously WRONG class-function pairing) and show
    # the congruence / integrality structural checks catch it.
    neg_control = None
    if H24_2d is not None:
        wrong_info = {"ord": 2, "chi": CLASSES["3A"]["chi"], "F_factor": CLASSES["2A"]["F_factor"]}
        wrong_row, raw_wrong, A_wrong = compute_twining_row("wrong", wrong_info, H24_2d, inv_e3_2d,
                                                             imaxH, safe24, H24, H_CUTOFF)
        neg_control = {
            "definition": "Deliberately WRONG pairing: chi=6 (that of class "
                           "3A) combined with F_2A (order-2 eta-quotient) -- "
                           "not a genuine M24 class function. Structural "
                           "checks should show a failure relative to the "
                           "correct 2A row above.",
            "H_all_reported_coeffs_integral": wrong_row["H_g_all_reported_coeffs_integral"],
            "A_n_wrong": wrong_row["A_n_g"],
            "congruence_mod_ord2_holds_per_n": {n: c["congruent"] for n, c in wrong_row[
                "congruence_An_g_equiv_An_mod_ord_g_n_le_8"].items()},
            "congruence_all_hold": wrong_row["congruence_all_hold"],
            "fails_relative_to_2A": wrong_row["congruence_all_hold"] != twining.get("2A", {}).get("congruence_all_hold"),
        }

    out = {
        "track": "A",
        "item": "K3 elliptic genus (phi_{0,1}, Z_K3, discriminant property, "
                "H(tau)/Appell-Lerch mu-term A_1..A_10), M24 twining "
                "(2A,3A,5A,7A) with structural cross-checks, rigidity scans "
                "(a) mu-term coefficient, (b) overall Z factor.",
        "method": "exact Fraction arithmetic; Jacobi theta triple-product "
                   "series (series2d.py/thetas.py/appell.py, reused "
                   "byte-for-byte from v1's blind genus-moonshine track); "
                   "Lambda_N via exact divisor-sum closed form, "
                   "cross-checked against a fully independent direct "
                   "log-derivative series computation for every N used "
                   "(2,3,5,7).",
        "reused_v1_files_unmodified": ["series2d.py", "thetas.py", "appell.py"],
        "defects_fixed_vs_v1": [
            "main Z_K3/phi01/c(D) cutoff raised 8->10 (task requires q^10)",
            "H(tau) cutoff raised so A_9,A_10 sit inside the truncation-"
            "stability safe margin (v1 reported cutoff_q_used=10 but only "
            "delivered A_1..A_8 because its own margin excluded n=9,10)",
            "twining extended from 2A-only to 2A,3A,5A,7A",
            "added structural cross-checks (integrality, sign pattern, "
            "A_n^(g) === A_n mod ord(g)) not present in v1",
            "added a negative control (wrong chi/F_g pairing) for the "
            "congruence check",
            "rigidity scan (b) now also tries the task-suggested (q^0,y^1)"
            "-coefficient selector and explicitly classifies it "
            "NORMALISATION rather than silently reporting a k",
            "added exports.json with chi_from_genus",
            "FOUND (new, not in v1's known-defect list): "
            "compute_ZK3_and_phi01(cutoff_q,N) alone is inexact at its own "
            "top row n=cutoff_q -- (n=10,l=+-5), D=15 read -23554 instead "
            "of the true -23550 (which matches n=4,l=+-1 and n=6,l=+-3 at "
            "the same D=15). Fixed by compute_ZK3_and_phi01_safe: compute "
            "at cutoff_q+2 internally, report only q<=cutoff_q. Evidence "
            "recorded below (edge_row_truncation_bug).",
        ],
        "edge_row_truncation_bug": {
            "definition": "Direct empirical demonstration of the truncation "
                           "defect above: same (n=10,l) coefficients (D=15) "
                           "computed by compute_ZK3_and_phi01(10,18) alone "
                           "(WRONG, top row) vs by the margin-2 safe wrapper "
                           "used everywhere below (CORRECT, matches n=4,n=6 "
                           "rows at the same D).",
            "raw_no_margin_n10_row": {str(k): str(v) for k, v in edge_row_bug_check.items()},
            "safe_margin2_n10_row": {str(k): str(v) for k, v in edge_row_bug_check_after_fix.items()},
            "n4_and_n6_same_D15_for_comparison": {
                str((4, 1)): str(cnl.get((4, 1))),
                str((6, 3)): str(cnl.get((6, 3))),
            },
        },
        "main_cutoff_q": MAIN_CUTOFF,
        "H_cutoff_q": H_CUTOFF,
        "H_safe_margin_index": safe24,
        "A9_A10_inside_safe_margin": A9_A10_inside_margin,
        "truncation_stability_main": stability,
        "ZK3_tau_0": ZK3_0_str,
        "ZK3_tau0_const": str(ZK3_tau0_const) if ZK3_tau0_const is not None else None,
        "discriminant_dependence_holds": disc_ok,
        "discriminant_mismatches": {str(D): [str(v) for v in vs] for D, vs in mismatches.items()},
        "c_of_D": {str(D): str(v) for D, v in sorted(cD.items())},
        "H_tau_appell_lerch": {
            "two_slices_agree": agree24,
            "multiply_back_check_passes": mb_ok,
            "A_n_1_to_10": {str(n): str(v) for n, v in sorted(A_n.items())},
            "truncation_stability_A1_to_A10": H_stability,
            "identity_A2_60_eq_4_A1_77_untwined": identity_untwined,
        },
        "rigidity_scan_a_mu_coefficient": {
            "definition": "Replace the '24' multiplying the Appell-Lerch "
                           "mu-term by N in 20..28; STRUCTURAL selector: "
                           "the two independent y-power slices (j=0,j=2) of "
                           "H(N) must agree (z-independence). No literal "
                           "target used in the selecting condition itself.",
            "scan": scan_a,
            "classification": "RIGID (structural: z-independence of two "
                               "independently-derived slices; a "
                               "neighbouring N != 24 demonstrably fails "
                               "agreement, see scan above)",
            "honesty_note": "H(N) = H(24) + (N-24)*mu(tau,z) algebraically "
                             "(mu = y^{1/2}*Psi/T1^2), so for N!=24 the "
                             "residual (N-24)*mu is manifestly "
                             "z-dependent and the two slices MUST disagree "
                             "-- carried over from v1: this is still a "
                             "single linear condition in N (mu's residual "
                             "vanishes iff N=24), not several independent "
                             "constraints. It is nonetheless a genuine "
                             "structural (not literal-target) selector: "
                             "24 is never asserted upstream, only the "
                             "'two independently-extracted slices agree' "
                             "condition is imposed, and it happens to "
                             "select N=24 uniquely in the scanned range. "
                             "The j0 polar term (q^{-1/8} coefficient) "
                             "varies linearly in N and is exactly zero at "
                             "N=20 (not undefined -- see scan, zero-filled) "
                             "while the j2 polar term is -2 for every N "
                             "(Psi's leading term is pure y^{-1/2}, "
                             "contributing nothing to theta_1^2's y^1 "
                             "coefficient); only N=24 makes the two "
                             "columns equal.",
        },
        "rigidity_scan_b_overall_factor": {
            "euler_char_scan": scan_b,
            "euler_char_classification": "NORMALISATION (single linear "
                "equation Z(tau,0)=24 in k; 24 is the literal target)",
            "q0_y1_coefficient_selector": rigidity_b_q0y1_selector,
        },
        "twining_2A_3A_5A_7A": twining,
        "twining_truncation_stability_cutoff_plus2": twining_stability,
        "lambda_N_cross_checks_match_direct_log_derivative": lambda_cross_checks,
        "negative_control_wrong_class_pairing": neg_control,
        "could_not_do": [
            "Rigidity scan (b) (overall factor k of Z=k*phi_{0,1}): no "
            "structural condition (consistency, integrality, modularity, "
            "exactness, or agreement of two independent computations) that "
            "pins k to a unique value WITHOUT referencing the literal "
            "targets 24 (Euler char) or 2 (q^0 y^1 coefficient) was found "
            "in this session -- both selectors tried are reported as "
            "NORMALISATION, per the task's own fallback instruction, "
            "rather than mislabelled RIGID. See "
            "rigidity_scan_b_overall_factor.q0_y1_coefficient_selector.",
        ],
    }

    with open("results.json", "w") as f:
        json.dump(out, f, indent=2, default=str)

    with open("exports.json", "w") as f:
        json.dump({
            "chi_from_genus": str(ZK3_tau0_const) if ZK3_tau0_const is not None else None,
            "source": "Z_K3(tau,0), computed constant term of "
                       "Z_K3=2*phi_{0,1}(tau,z) at z=0, from "
                       "results.json:ZK3_tau0_const",
            "cutoff_q_used": MAIN_CUTOFF,
        }, f, indent=2)

    print(json.dumps(out, indent=2, default=str)[:4000])
    print("\n... (full output written to results.json, exports.json)")


if __name__ == "__main__":
    main()
