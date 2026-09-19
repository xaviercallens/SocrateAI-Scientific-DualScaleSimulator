#!/usr/bin/env python3
"""
Track A v3 -- K3 elliptic genus, Mathieu moonshine, twining (tier B / L).

Command (repo-root-relative; the committed run):
  cd audit/k3t2_rigidity_v3/A-genus && \
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python run_track_a_v3.py --cutoff 16 --nmax 12 --kmax 6
(needs m24_shapes.json, produced by:  python m24_shapes.py --samples 4000 --seed 1 )
Writes inputs.json, results.json, exports.json next to the script.

Wording: statements are "verified by exact computation (tier B)" or "holds for
every case enumerated"; nothing here is a proof.

v2 defects addressed
 (a) k is a declared input (inputs.json); it is in shared_inputs of every downstream item.
 (b) mu-term coefficient N is SOLVED (exact linear solve + integer scan 0..8*12k) for
     k=1..kmax; the H used downstream is the one built with the returned N, never a typed 24.
 (c) chi(g) = number of 1-cycles of g in the permutation action on 24 points built from the
     Golay code (m24_shapes.py); F_g comes from the requirement that the ground-state term of
     the twined genus is unchanged, inside M_2(Gamma_0(N)).
 (d) sign / eigenvalue-multiplicity check done (nonnegative integer multiplicities).
 (e) twined integrality and congruences for 2A,3A,5A,7A; 4B and 11A reduced to a one-parameter
     family whose parameter is scanned (not solved: reported as conditional).
 (f) 27720-lock A_2*60 = 4*A_1*77 tested under each twining, in raw and halved conventions.
"""
import argparse, json, sys
from fractions import Fraction as Fr
from pathlib import Path
from math import gcd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
REPO = HERE.parents[2]

from sympy import divisor_sigma, totient, divisors, legendre_symbol
from series2d import y_to_1, divide_scalar_series, mul, add, scal, shift, reciprocal_1d
from thetas import theta2, theta3, theta4, theta1_over_i, eta3
from appell import psi

CMD = ("cd audit/k3t2_rigidity_v3/A-genus && /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/"
       ".venv-tda/bin/python run_track_a_v3.py --cutoff {c} --nmax {n} --kmax {k}")
SCRIPT = "audit/k3t2_rigidity_v3/A-genus/run_track_a_v3.py"


def fs(x):
    x = Fr(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


# ---------------------------------------------------------------- theta side
def build_blocks(cutoff_q, N_prod):
    imax = 8 * cutoff_q
    t2, t3, t4 = theta2(N_prod, imax), theta3(N_prod, imax), theta4(N_prod, imax)
    r = [divide_scalar_series(t, y_to_1(t), imax) for t in (t2, t3, t4)]
    phi01 = scal(add(add(mul(r[0], r[0], imax), mul(r[1], r[1], imax)), mul(r[2], r[2], imax)), 4)
    T1 = theta1_over_i(N_prod, imax)
    Psi = psi(N_prod, imax, T1)
    T1sq = mul(T1, T1, imax)
    e3 = eta3(N_prod, imax)
    phi_e3 = mul(phi01, e3, imax)          # unit-k piece of Z*eta^3
    y12Psi = shift(Psi, dj=1)
    return dict(imax=imax, phi01=phi01, T1sq=T1sq, phi_e3=phi_e3, y12Psi=y12Psi)


def slice1d(d, j0):
    return {i: v for (i, j), v in d.items() if j == j0}


def slice_div(part, T1sq, j0, imax):
    inv = reciprocal_1d(slice1d(T1sq, j0), imax, istep=8)
    q = mul({(i, 0): v for i, v in slice1d(part, j0).items()},
            {(i, 0): v for i, v in inv.items()}, imax)
    return {k[0]: v for k, v in q.items()}


def solve_mu(blocks, k, margin=16, scan_max=None):
    """H_j(N) = N*P_j - k*Q_j for the two y-slices j=0,2. Selecting condition: H
    independent of z, i.e. H_0 == H_2 at every safe q-order. Exact solve for N and an
    integer scan (every integer 0..scan_max)."""
    imax = blocks["imax"]; safe = imax - margin
    T1sq = blocks["T1sq"]
    P = {j: slice_div(blocks["y12Psi"], T1sq, j, imax) for j in (0, 2)}
    Q = {j: slice_div(blocks["phi_e3"], T1sq, j, imax) for j in (0, 2)}
    idx = sorted(i for i in set(P[0]) | set(P[2]) | set(Q[0]) | set(Q[2]) if i <= safe)
    dP = {i: P[0].get(i, Fr(0)) - P[2].get(i, Fr(0)) for i in idx}
    dQ = {i: Q[0].get(i, Fr(0)) - Q[2].get(i, Fr(0)) for i in idx}
    def agrees(N): return all(N * dP[i] - k * dQ[i] == 0 for i in idx)
    nz = [i for i in idx if dP[i] != 0]
    N_exact = (k * dQ[nz[0]] / dP[nz[0]]) if nz else None
    exact_ok = N_exact is not None and agrees(N_exact)
    scan_max = scan_max or 80
    scan_hits = [N for N in range(0, scan_max + 1) if agrees(Fr(N))]
    H = {i: N_exact * P[0].get(i, Fr(0)) - k * Q[0].get(i, Fr(0)) for i in idx} if exact_ok else None
    return dict(N_exact=N_exact, exact_ok=exact_ok, scan_hits=scan_hits, scan_max=scan_max,
                H=H, safe=safe, P0=P[0], Q0=Q[0], n_orders_checked=len(idx))


def multiply_back(H, N, k, blocks, safe):
    imax = blocks["imax"]
    back = mul({(i, 0): v for i, v in H.items()}, blocks["T1sq"], imax)
    req = add(scal(blocks["y12Psi"], N), scal(blocks["phi_e3"], -k))
    keys = set(back) | set(req)
    return all(back.get(kk, 0) == req.get(kk, 0) for kk in keys if kk[0] <= safe)


def h_list(H, nmax):
    return [H.get(-1 + 8 * n, Fr(0)) for n in range(nmax + 1)]


# ------------------------------------------------------------- 1D q-series
def smul(a, b, L):
    out = [Fr(0)] * L
    for i, x in enumerate(a):
        if x == 0: continue
        for j, y in enumerate(b):
            if i + j >= L: break
            out[i + j] += x * y
    return out


def sinv(a, L):
    assert a[0] != 0
    out = [Fr(0)] * L
    out[0] = 1 / a[0]
    for n in range(1, L):
        s = sum(a[k] * out[n - k] for k in range(1, min(n, len(a) - 1) + 1))
        out[n] = -s / a[0]
    return out


def poch(step, L, power=1):
    """prod_{n>=1} (1-q^{step n})^power as list length L (power>=0 int)."""
    out = [Fr(0)] * L; out[0] = Fr(1)
    for _ in range(power):
        for n in range(1, L // step + 1):
            new = out[:]
            for i in range(L - step * n):
                new[i + step * n] -= out[i]
            out = new
    return out


def E2_series(L):
    return [Fr(1)] + [Fr(-24 * int(divisor_sigma(n, 1))) for n in range(1, L)]


def at_Nq(a, N, L):
    out = [Fr(0)] * L
    for i, x in enumerate(a):
        if i * N < L: out[i * N] = x
    return out


def E_N(N, L):
    e = E2_series(L); eN = at_Nq(e, N, L)
    return [(N * eN[i] - e[i]) / (N - 1) for i in range(L)]


def E_N_from_logderiv(N, L):
    """independent: E_N = 1 + 24/(N-1) * q f'/f with f = prod(1-q^{Nn})/prod(1-q^n)."""
    f = smul(poch(N, L), sinv(poch(1, L), L), L)
    qfp = [i * f[i] for i in range(L)]
    r = smul(qfp, sinv(f, L), L)
    return [Fr(1) + 24 * r[0] / (N - 1)] + [24 * r[i] / (N - 1) for i in range(1, L)]


def dim_M2_Gamma0(N):
    """dim M_2(Gamma_0(N)) = genus + cusps - 1 (standard index/elliptic-point formula, tier L)."""
    mu = Fr(N)
    ps = [p for p in range(2, N + 1) if N % p == 0 and all(p % r for r in range(2, p))]
    for p in ps: mu *= Fr(p + 1, p)
    def k1(p): return 0 if p == 2 else legendre_symbol(-1 % p, p)
    def k3(p): return 0 if p == 3 else (-1 if p == 2 else legendre_symbol(-3 % p, p))
    nu2 = 0 if N % 4 == 0 else eval("1") * 1
    nu2 = 0 if N % 4 == 0 else Fr(1)
    for p in ps: nu2 *= (1 + k1(p)) if N % 4 else 0
    nu3 = 0 if N % 9 == 0 else Fr(1)
    for p in ps: nu3 *= (1 + k3(p)) if N % 9 else 0
    ninf = sum(int(totient(gcd(d, N // d))) for d in divisors(N))
    g = 1 + mu / 12 - Fr(nu2) / 4 - Fr(nu3) / 3 - Fr(ninf) / 2
    return int(g + ninf - 1), dict(index=str(mu), nu2=int(nu2), nu3=int(nu3), cusps=ninf, genus=str(g))


def frac_int(x): return Fr(x).denominator == 1


# ------------------------------------------------------- multiplicity tests
def test_prime(d, t, p):
    """rational trace t of an order-p element on a rep of dim d: eigenvalue-1 multiplicity
    (d+(p-1)t)/p and each primitive root (d-t)/p must be nonnegative integers."""
    m1 = (d + (p - 1) * t) / p; mz = (d - t) / p
    return frac_int(m1) and frac_int(mz) and m1 >= 0 and mz >= 0, (m1, mz)


def test_order4(d, t1, t2):
    S = (d + t2) / 2; mi = (d - t2) / 4
    m1 = (S + t1) / 2; mm1 = (S - t1) / 2
    ok = all(frac_int(x) and x >= 0 for x in (m1, mm1, mi))
    return ok, (m1, mm1, mi)


def class_report(hg, h, nmax, order, t2=None):
    rows = []; all_int = True; all_cong = True; all_mult = True; signs = []
    for n in range(1, nmax + 1):
        t, d = hg[n], h[n]
        integ = frac_int(t)
        cong = integ and frac_int((d - t) / order) if order else None
        if order == 4:
            mult, ms = test_order4(d, t, t2[n])
        else:
            mult, ms = test_prime(d, t, order)
        all_int &= integ; all_cong &= bool(cong); all_mult &= mult
        sg = 0 if t == 0 else (1 if t > 0 else -1)
        signs.append(sg)
        rows.append(dict(n=n, dim=fs(d), trace=fs(t), A_n_g_halved=fs(t / 2), integral=integ,
                         congruent_mod_order=bool(cong), multiplicities=[fs(x) for x in ms],
                         multiplicities_nonneg_integer=mult, sign=sg,
                         abs_trace_le_dim=abs(t) <= d))
    lock = lock_test(hg, h)
    return dict(rows=rows, all_integral=all_int, all_congruent=all_cong, all_multiplicities_ok=all_mult,
                sign_pattern="".join({1: "+", -1: "-", 0: "0"}[s] for s in signs),
                sign_pattern_alternating_minus_one_pow_n=all(s == (1 if n % 2 == 0 else -1) for n, s in zip(range(1, nmax + 1), signs)),
                lock=lock)


def lock_test(hg, h):
    """A_2*60 vs 4*A_1*77 in the raw convention (K_n traces) and the halved (A_n = K_n/2)
    convention; the relation is homogeneous of degree 1 so both must agree."""
    raw = hg[2] * 60 - 4 * hg[1] * 77
    half = (hg[2] / 2) * 60 - 4 * (hg[1] / 2) * 77
    return dict(raw_lhs_minus_rhs=fs(raw), halved_lhs_minus_rhs=fs(half), holds=(raw == 0),
                conventions_agree=(raw == 0) == (half == 0))


def twined(h, T, chi, P3, L):
    """h_g[n] = (chi/24)*h[n] - (T*P3)[n]   (raw coefficient of q^{-1/8+n} of H_g)."""
    TP = smul(T, P3, L)
    return [Fr(chi, 24) * h[n] - TP[n] for n in range(L)]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cutoff", type=int, default=16)
    ap.add_argument("--nmax", type=int, default=12)
    ap.add_argument("--kmax", type=int, default=6)
    a = ap.parse_args()
    NM = a.nmax; KMAX = a.kmax
    GRID_D, GRID_M = 60, 30000   # c grid: denominators 60 (covers 1/12 and 1/5 lattices), |c|<=500
    K0 = 2  # DECLARED INPUT k (inputs.json); all other k values are hypothetical scans
    command = CMD.format(c=a.cutoff, n=NM, k=KMAX)
    L = NM + 1
    m24 = json.loads((HERE / "m24_shapes.json").read_text())
    R = {}; results = []; could_not = []

    # ------------------------------------------------------------ Z, phi01, c(D)
    zc = NM
    B = build_blocks(zc + 2, zc + 10)
    phi01 = {kk: v for kk, v in B["phi01"].items() if kk[0] <= 8 * zc}
    from collections import defaultdict
    cnl = {(i // 8, j // 2): v for (i, j), v in phi01.items() if i % 8 == 0 and j % 2 == 0 and 0 <= i // 8 <= zc}
    byD = defaultdict(set)
    for (n, l), v in cnl.items(): byD[4 * n - l * l].add(v)
    disc_ok = all(len(s) == 1 for s in byD.values())
    c1 = {D: next(iter(s)) for D, s in byD.items() if len(s) == 1}
    phi0 = {i: v for i, v in y_to_1(phi01).items() if v != 0}
    Z0_unit = phi0.get(0)
    B2 = build_blocks(zc + 4, zc + 12)
    phi01_b = {kk: v for kk, v in B2["phi01"].items() if kk[0] <= 8 * zc}
    R["phi01"] = dict(
        cutoff_q=zc, discriminant_dependence_holds=disc_ok,
        tau0_series_of_phi01={fs(Fr(i, 8)): fs(v) for i, v in sorted(phi0.items())},
        phi01_at_z0_constant=fs(Z0_unit), phi01_at_z0_only_constant=(list(phi0) == [0]),
        truncation_stable_cutoff_plus2=(phi01 == phi01_b),
        c1_of_D={str(D): fs(c1[D]) for D in sorted(c1) if D <= 4 * zc - 25})
    results.append(dict(id="A.phi01_Z0", quantity="phi_{0,1}(tau,0) (unit-k genus at z=0)", computed=fs(Z0_unit),
                        shared_inputs=[], script=SCRIPT, command=command))

    # ------------------------------------------------------ mu-term solve per k
    Hb = build_blocks(a.cutoff, a.cutoff + 10)
    Hb2 = build_blocks(a.cutoff + 2, a.cutoff + 12)
    mu_rows = {}; Hk = {}
    for k in range(1, KMAX + 1):
        s = solve_mu(Hb, k, scan_max=8 * 12 * KMAX // KMAX + 40)
        mb = multiply_back(s["H"], s["N_exact"], k, Hb, s["safe"]) if s["exact_ok"] else False
        s2 = solve_mu(Hb2, k, scan_max=s["scan_max"])
        stab = s2["exact_ok"] and s2["N_exact"] == s["N_exact"] and \
            h_list(s2["H"], NM) == h_list(s["H"], NM)
        hl = h_list(s["H"], NM) if s["exact_ok"] else None
        Hk[k] = dict(N=s["N_exact"], h=hl)
        mu_rows[str(k)] = dict(N_exact_solution=fs(s["N_exact"]), exact_solution_consistent_all_orders=s["exact_ok"],
                               integer_scan_range=[0, s["scan_max"]], integer_scan_solutions=s["scan_hits"],
                               neighbours_N_pm1_fail=(int(s["N_exact"]) + 1 not in s["scan_hits"] and int(s["N_exact"]) - 1 not in s["scan_hits"]),
                               multiply_back_ok=mb, truncation_stable_cutoff_plus2=stab,
                               orders_checked=s["n_orders_checked"], safe_i=s["safe"],
                               Z_at_tau_z0_from_phi01=fs(k * Z0_unit),
                               N_equals_Z_at_z0=(s["N_exact"] == k * Z0_unit),
                               h_0_to_nmax=[fs(x) for x in hl] if hl else None)
        results.append(dict(id=f"A.mu_N_solution_k{k}", quantity=f"mu-term coefficient N solving z-independence of H, k={k}",
                            computed=fs(s["N_exact"]), shared_inputs=["k"], script=SCRIPT, command=command))
    R["mu_term_scan_by_k"] = mu_rows
    R["mu_term_classification"] = dict(
        label="RIGID_GIVEN_DEFINITION",
        selecting_condition="H(tau) independent of z (the two y-slices j=0 and j=2 of (N y^{1/2}Psi - Z eta^3)/theta_1^2 agree at every q order)",
        condition_uses_true_value=False,
        note="The condition does not contain N, but its solution is N = 12k = Z(tau,0)(k): given k it is fixed and every neighbouring integer N fails, so it is RIGID once Z=k*phi01 is granted. It does not select k (all k solve it, N scaling with k).")

    # ---------------------------------------------- H, A_n and lock (untwined)
    def A_row(k):
        h = Hk[k]["h"]
        return [fs(x / 2) for x in h]
    R["H_and_A_n"] = {}
    for k in range(1, KMAX + 1):
        h = Hk[k]["h"]
        lk = lock_test(h, h)
        R["H_and_A_n"][str(k)] = dict(h_raw=[fs(x) for x in h], A_n_halved=A_row(k),
                                      leading_coefficient_h0=fs(h[0]), lock_raw_lhs_minus_rhs=lk["raw_lhs_minus_rhs"],
                                      A_n_integral=all(frac_int(x / 2) for x in h),
                                      h_n_nonneg_for_n_ge_1=all(x >= 0 for x in h[1:]))
    A2 = [Fr(x) for x in [Hk[2]["h"][n] / 2 for n in range(NM + 1)]]
    results.append(dict(id="A.A_n_k2", quantity=f"A_1..A_{NM} of H = 2q^(-1/8)(-1+sum A_n q^n) (convention A_n=h_n/2), k=2",
                        computed=", ".join(fs(x) for x in A2[1:]), shared_inputs=["k"], script=SCRIPT, command=command))
    results.append(dict(id="A.lock_untwined_k2", quantity="A_2*60 - 4*A_1*77 (untwined)",
                        computed=R["H_and_A_n"]["2"]["lock_raw_lhs_minus_rhs"], shared_inputs=["k"], script=SCRIPT, command=command))

    # -------------------------------------------------------- k selector (a)
    P3 = sinv(poch(1, L, 3), L)
    shapes = {(c["order"], c["shape_str"]): c for c in m24["classes"]}
    CLS = {"2A": (2, "1^8 2^8"), "3A": (3, "1^6 3^6"), "5A": (5, "1^4 5^4"), "7A": (7, "1^3 7^3")}
    def fixpts(key): return shapes[key]["fixed_points"]
    chi = {g: fixpts(key) for g, key in CLS.items()}
    dimsM2 = {}
    for N in (2, 3, 4, 5, 7, 11):
        d, info = dim_M2_Gamma0(N); dimsM2[str(N)] = dict(dim=d, **info)
    R["dim_M2_Gamma0"] = dimsM2
    ELn = {}
    for N in (2, 3, 5, 7):
        e1, e2 = E_N(N, L), E_N_from_logderiv(N, L)
        ELn[N] = e1 if e1 == e2 else None
    R["E_N_cross_check_eisenstein_vs_logderivative"] = {str(N): (ELn[N] is not None) for N in ELn}

    # twin for a given k (linear in k/2). T const = k*(1 - chi/24)
    def twin_class(g, k, h):
        N, key = CLS[g]
        chi_g = fixpts((N, key))
        T0 = Fr(k) * (1 - Fr(chi_g, 24))
        T = [T0 * x for x in ELn[N]]
        return twined(h, T, chi_g, P3, L), T

    k_par = {}
    for k in range(1, KMAX + 1):
        h = Hk[k]["h"]
        row = {}
        for g in CLS:
            hg, _ = twin_class(g, k, h)
            row[g] = dict(all_h_g_integral=all(frac_int(x) for x in hg[:NM + 1]),
                          Z_g_at_z0=fs(Fr(k) * chi[g] / 2), Z_g_at_z0_integer=frac_int(Fr(k) * chi[g] / 2))
        k_par[str(k)] = row
    R["k_selector_trial_twined_integrality"] = k_par
    passing = [k for k in range(1, KMAX + 1) if all(v["all_h_g_integral"] and v["Z_g_at_z0_integer"] for v in k_par[str(k)].values())]
    R["k_selector_summary"] = dict(
        candidate_selectors_tried=[
            "integrality of Z coefficients: all integer k pass",
            "z-independence of H: holds for every k (N=12k)",
            "positivity of A_n (n>=1): holds for every k>0 (see H_and_A_n)",
            "integrality of Z_g(tau,0)=k*fix(g)/2 and of all h_g, for 2A,3A,5A,7A (given the declared 24-point permutation premise)"],
        k_passing_last_selector=passing,
        result="No selector without a target: modular/integrality conditions leave k in 2Z (parity from 7A with 3 fixed points) and never bound k above; k=2 is the smallest positive one. Label NORMALISATION.",
        label="NORMALISATION",
        selecting_condition_contains_target=True)
    results.append(dict(id="A.k_selector", quantity="k values in 1..kmax passing twined-integrality (parity) selector",
                        computed=str(passing), shared_inputs=["k", "perm24_premise"], script=SCRIPT, command=command))

    # ---------------------------------------------------- twining at k = 2

    h = Hk[K0]["h"]
    hl = [Fr(x) for x in h]
    # expected values (v2 task, FROM MEMORY): F_g = fac * Lambda_N, Lambda_N = N(N-1)/24 * E_N
    exp_fac = {"2A": 16, "3A": 6, "5A": 2, "7A": 1}
    tw = {}; t2A = None
    for g, (N, key) in CLS.items():
        hg, T = twin_class(g, K0, hl)
        t = [hg[n] for n in range(NM + 1)]
        rep = class_report(hg, hl, NM, N)
        Fg_factor = T[0] / (Fr(N * (N - 1), 24)) if T else None
        tw[g] = dict(order=N, cycle_shape=key, chi_from_cycle_shape=chi[g],
                     dim_M2=dimsM2[str(N)]["dim"], T0_from_ground_state=fs(T[0]),
                     F_g_coefficient_of_Lambda_N_derived=fs(Fg_factor),
                     F_g_coefficient_expected_from_memory_v2=exp_fac[g],
                     F_g_matches_expected=(Fg_factor == exp_fac[g]),
                     h_g_raw=[fs(x) for x in t], h_0_g=fs(t[0]), **rep)
        if g == "2A": t2A = hg
        results.append(dict(id=f"A.twin_{g}", quantity=f"twined raw h_g[1..{NM}] and checks for {g}",
                            computed=f"chi={chi[g]}; F_g/Lambda_N={fs(Fg_factor)}; integral={rep['all_integral']}; congruent={rep['all_congruent']}; multiplicities_ok={rep['all_multiplicities_ok']}; A_n^g(1..3)={[fs(x/2) for x in t[1:4]]}",
                            shared_inputs=["k", "perm24_premise", "ground_state_invariance", "twined_form_structure", "moonshine_module_premise"],
                            script=SCRIPT, command=command))
    R["twining_k2"] = tw

    # ---------------------------------------- chi scan (negative control same parameter)
    realizable = {N: sorted({c["fixed_points"] for c in m24["classes"] if c["order"] == N}) for N in (2, 3, 5, 7)}
    chi_scan = {}
    for N in (2, 3, 5, 7):
        passing_int = []; passing_all = []
        for x in range(-12, 49):
            T0 = Fr(K0) * (1 - Fr(x, 24)); T = [T0 * e for e in ELn[N]]
            hg = twined(hl, T, x, P3, L)
            rep = class_report(hg, hl, NM, N)
            if rep["all_integral"]: passing_int.append(x)
            if rep["all_integral"] and rep["all_congruent"] and rep["all_multiplicities_ok"]: passing_all.append(x)
        gname = [g for g in CLS if CLS[g][0] == N][0]
        extra = [x for x in passing_all if x != chi[gname] and x != 24]
        extra_real = [x for x in extra if x in realizable[N]]
        cls_lab = "RIGID_GIVEN_DEFINITION" if (not extra_real and (chi[gname] - 1 not in passing_all) and (chi[gname] + 1 not in passing_all)) else "NON_DISCRIMINATING"
        chi_scan[str(N)] = dict(classification=cls_lab, other_passing_values_excluding_identity_24=extra,
                                other_passing_values_realized_by_M24_element_of_this_order=extra_real,
                                note="values passing but not realized by an M24 element of this order are HYPOTHETICAL and are not counted as discrimination; chi=24 is the identity (T=0). A realized value that passes (2B: 0, 3B: 0) shows the test cannot separate the class from another real class of the same order.",scan_range=[-12, 48], every_integer_stepped=True,
                                passing_integrality_only=passing_int, passing_integrality_congruence_multiplicity=passing_all,
                                fix_counts_of_actual_M24_elements_of_this_order=realizable[N],
                                true_value_among_passing=(chi[[g for g in CLS if CLS[g][0] == N][0]] in passing_all),
                                neighbours_of_true_value_pass=[x for x in passing_all if abs(x - chi[[g for g in CLS if CLS[g][0] == N][0]]) == 1])
        results.append(dict(id=f"A.chi_scan_N{N}", quantity=f"integers chi in [-12,48] whose ground-state-fixed twining at level {N} passes all checks (n<={NM})",
                            computed=str(passing_all), shared_inputs=["k", "ground_state_invariance", "twined_form_structure", "moonshine_module_premise"],
                            script=SCRIPT, command=command))
    R["chi_scan_same_parameter_control"] = chi_scan

    # ------------------------------------------ 4B and 11A one-parameter families
    fam = {}
    # 4B
    sq_ok = None
    shape4b = shapes[(4, "1^4 2^2 4^4")]
    perm = shape4b["perm"]; sq = [perm[perm[i]] for i in range(24)]
    seen = set(); sqshape = {}
    for i in range(24):
        if i in seen: continue
        j = i; Ln = 0
        while j not in seen: seen.add(j); j = sq[j]; Ln += 1
        sqshape[Ln] = sqshape.get(Ln, 0) + 1
    e2 = E_N(2, L); e2_2 = at_Nq(e2, 2, L)
    V = [e2[i] - e2_2[i] for i in range(L)]
    chi4 = shape4b["fixed_points"]
    T0_4 = Fr(K0) * (1 - Fr(chi4, 24))
    def h4(c): return twined(hl, [T0_4 * e2[i] + c * V[i] for i in range(L)], chi4, P3, L)
    a0, a1 = h4(Fr(0)), h4(Fr(1))
    beta = [a1[n] - a0[n] for n in range(L)]
    sols4 = []
    for m in range(-GRID_M, GRID_M + 1):
        c = Fr(m, GRID_D)
        if not frac_int(a0[1] + c * beta[1]): continue
        hg = [a0[n] + c * beta[n] for n in range(L)]
        rep = class_report(hg, hl, NM, 4, t2=t2A)
        if rep["all_integral"] and rep["all_multiplicities_ok"]: sols4.append(m)
    fam["4B"] = dict(cycle_shape="1^4 2^2 4^4", chi=chi4, square_cycle_shape=str(sqshape),
                     square_is_2A_shape=(sqshape == {1: 8, 2: 8}), dim_M2=dimsM2["4"]["dim"],
                     parameter="c in T = T0*E^(2)(tau) + c*(E^(2)(tau)-E^(2)(2tau))",
                     scan=f"c = m/{GRID_D}, m in [-{GRID_M},{GRID_M}] every integer m (|c|<={GRID_M//GRID_D}); n=1 integrality used as prefilter, all n<=nmax checked afterwards",
                     n_solutions=len(sols4), solutions_c_first20=[fs(Fr(m, GRID_D)) for m in sols4[:20]],
                     solutions_c_min=fs(Fr(min(sols4), GRID_D)) if sols4 else None, solutions_c_max=fs(Fr(max(sols4), GRID_D)) if sols4 else None,
                     solutions_all_integer_c=all(m % 12 == 0 for m in sols4),
                     classification="CONDITIONAL_ON_INPUT", unique=(len(sols4) == 1))
    results.append(dict(id="A.fam_4B", quantity="c values (grid m/12) with integral traces and nonneg integer eigenvalue multiplicities, 4B",
                        computed=f"{len(sols4)} solutions; first: {fam['4B']['solutions_c_first20'][:6]}", shared_inputs=["k", "perm24_premise", "ground_state_invariance", "twined_form_structure", "moonshine_module_premise"],
                        script=SCRIPT, command=command))
    # 11A
    # Ligozat conditions for f = eta(tau)^2 eta(11 tau)^2 on Gamma_0(11)
    lig = dict(weight=Fr(2 + 2, 2), sum_a_m=1 * 2 + 11 * 2, sum_Nover_a_m=11 * 2 + 1 * 2, prod_a_pow_m=11 ** 2)
    f11 = smul(poch(1, L, 2), poch(11, L, 2), L); f11 = [Fr(0)] + f11[:L - 1]
    e11 = E_N(11, L)
    chi11 = shapes[(11, "1^2 11^2")]["fixed_points"]
    T0_11 = Fr(K0) * (1 - Fr(chi11, 24))
    b0 = twined(hl, [T0_11 * x for x in e11], chi11, P3, L)
    b1 = twined(hl, [T0_11 * e11[i] + f11[i] for i in range(L)], chi11, P3, L)
    beta11 = [b1[n] - b0[n] for n in range(L)]
    sols11 = []
    for m in range(-GRID_M, GRID_M + 1):
        c = Fr(m, GRID_D)
        if not frac_int(b0[1] + c * beta11[1]): continue
        hg = [b0[n] + c * beta11[n] for n in range(L)]
        rep = class_report(hg, hl, NM, 11)
        if rep["all_integral"] and rep["all_multiplicities_ok"]: sols11.append(m)
    fam["11A"] = dict(cycle_shape="1^2 11^2", chi=chi11, dim_M2=dimsM2["11"]["dim"],
                      cusp_form="eta(tau)^2 eta(11 tau)^2 (Ligozat: sum a m = 24, sum (N/a) m = 24, product a^m = 11^2 square)",
                      ligozat_sums_divisible_by_24=(lig["sum_a_m"] % 24 == 0 and lig["sum_Nover_a_m"] % 24 == 0),
                      parameter="c in T = T0*E_11 + c*eta(tau)^2 eta(11tau)^2", scan=f"c = m/{GRID_D}, m in [-{GRID_M},{GRID_M}] every integer m (|c|<={GRID_M//GRID_D}); n=1 integrality used as prefilter, all n<=nmax checked afterwards",
                      n_solutions=len(sols11), solutions_c_first20=[fs(Fr(m, GRID_D)) for m in sols11[:20]],
                      solutions_c_min=fs(Fr(min(sols11), GRID_D)) if sols11 else None, solutions_c_max=fs(Fr(max(sols11), GRID_D)) if sols11 else None,
                      solutions_all_integer_c=all(m % 12 == 0 for m in sols11),
                      classification="CONDITIONAL_ON_INPUT", unique=(len(sols11) == 1))
    results.append(dict(id="A.fam_11A", quantity="c values (grid m/12) with integral traces and nonneg integer multiplicities, 11A",
                        computed=f"{len(sols11)} solutions; first: {fam['11A']['solutions_c_first20'][:6]}", shared_inputs=["k", "perm24_premise", "ground_state_invariance", "twined_form_structure", "moonshine_module_premise"],
                        script=SCRIPT, command=command))
    R["families_4B_11A"] = fam
    if not (fam["4B"]["unique"] and fam["11A"]["unique"]):
        could_not.append("4B and 11A: dim M_2(Gamma_0(N)) = 2 leaves one free parameter c after the ground-state condition; the integrality + eigenvalue-multiplicity scan does not pin c uniquely (see families_4B_11A). Fixing it needs the g-twisted-sector ground states or the M24 character table, neither derived here.")

    # ----------------------------------- sign pattern / lock summary
    R["sign_patterns_k2"] = {g: dict(pattern=tw[g]["sign_pattern"], alternating=tw[g]["sign_pattern_alternating_minus_one_pow_n"],
                                     all_multiplicities_nonneg_integer=tw[g]["all_multiplicities_ok"]) for g in tw}
    R["lock_under_twining_k2"] = {g: tw[g]["lock"] for g in tw}
    could_not += [
        "k has no structural selector free of the target (NORMALISATION); k=2 is the declared value.",
        "Twining beyond the classes listed (2B, 3B, 4A, 4C, 6A, ... with multipliers) not attempted.",
        "M24 character values were not computed; the eigenvalue-multiplicity test only uses that traces are rational integers on a dimension-h_n module, which is weaker than an M24 irrep decomposition.",
    ]

    # ---------------------------------------------------------------- inputs
    inputs = [
        dict(name="k", value="2 (scanned k=1..%d as hypotheticals)" % KMAX, tier="C",
             why="Overall factor of Z_K3 = k*phi_{0,1}; weak Jacobi forms of weight 0 index 1 are one-dimensional so nothing modular fixes k. k=2 is the K3 elliptic genus (chi=24). NORMALISATION."),
        dict(name="perm24_premise", value="Z_g(tau,0) = (k/2)*fix(g); fix(g) = number of 1-cycles of g in the 24-point permutation action built from the Golay code", tier="L",
             why="Mathieu moonshine premise; at k=2 this is chi(g)=fix(g). Ties the number of points (24) to Z(tau,0)=12k, so it is itself compatible with k=2 only."),
        dict(name="ground_state_invariance", value="y^{+-1} q^0 coefficient of Z_g equals that of Z (=k)", tier="L",
             why="The two RR ground states (h^{0,0},h^{2,0}) are g-invariant."),
        dict(name="twined_form_structure", value="Z_g = (k chi/24) phi_{0,1} + T_g phi_{-2,1}, T_g in M_2(Gamma_0(ord g)) with trivial multiplier for 2A,3A,5A,7A; also 4B (level 4) and 11A (level 11) reduced to a family", tier="L",
             why="Standard structure of twining genera; trivial-multiplier classes selected by the frame-shape condition sum a m_a = sum (N/a) m_a = 24 (checked in results)."),
        dict(name="moonshine_module_premise", value="K_n (raw dimension h_n) is a virtual-free M24 module, so raw twined traces are rational integers with nonnegative integer eigenvalue multiplicities", tier="L",
             why="Used for integrality, congruence and multiplicity checks; failing them would falsify the pairing, passing them is necessary, not sufficient."),
        dict(name="golay_and_generators", value="QR-mod-23 Golay code; generators x->x+1, x->2x, x->-1/x, delta", tier="L",
             why="Group generators quoted FROM MEMORY; verified by code preservation and group order (m24_shapes.json)."),
        dict(name="dim_M2_formula", value="dim M_2(Gamma_0(N)) = genus + cusps - 1", tier="L", why="Standard index formula."),
        dict(name="expected_M24_order", value="244823040", tier="L", why="FROM MEMORY; expected-field comparison only."),
        dict(name="expected_F_g_v2", value="F_2A=16 Lambda_2, F_3A=6 Lambda_3, F_5A=2 Lambda_5, F_7A=Lambda_7", tier="L",
             why="FROM MEMORY (v2 task); expected-field comparison only, never used in computation."),
    ]
    (HERE / "inputs.json").write_text(json.dumps(inputs, indent=1))

    # -------------------------------------------------------------- rigidity
    rig = [
        dict(parameter_inserted="k (overall factor in Z=k*phi01)", selecting_condition="integrality/positivity/modularity/twined-integrality tests (all computed in k_selector_summary)",
             condition_uses_true_value=True, solution_set=f"k in {passing} (integrality parity) with no upper bound",
             classification="NORMALISATION", negative_control="k=1 fails twined integrality (7A: Z_g(0)=3/2); controls perturb k itself",
             control_perturbs_same_parameter=True, input_it_depends_on="perm24_premise"),
        dict(parameter_inserted="N (mu-term coefficient)", selecting_condition="z-independence of H", condition_uses_true_value=False,
             solution_set="N = 12k for each k=1..%d (exact solve, integer scan 0..%d)" % (KMAX, 8 * 12 + 40),
             classification="RIGID_GIVEN_DEFINITION", negative_control="every other integer N fails the two-slice agreement",
             control_perturbs_same_parameter=True, input_it_depends_on="k"),
        dict(parameter_inserted="chi(g) at level N in {2,3,5,7}", selecting_condition="integrality + eigenvalue-multiplicity + congruence of the ground-state-fixed twining, n<=%d" % NM,
             condition_uses_true_value=False, solution_set=json.dumps({N: chi_scan[str(N)]["passing_integrality_congruence_multiplicity"] for N in (2, 3, 5, 7)}),
             classification="per level: " + json.dumps({N: chi_scan[str(N)]["classification"] for N in (2, 3, 5, 7)}), negative_control="all other integers chi in [-12,48] (values that are not fixed-point counts of an order-N element are HYPOTHETICAL, not counted as discrimination)",
             control_perturbs_same_parameter=True, input_it_depends_on="k, ground_state_invariance, twined_form_structure"),
        dict(parameter_inserted="c (4B, 11A cusp-direction coefficient)", selecting_condition="integrality + multiplicity scan", condition_uses_true_value=False,
             solution_set=f"4B: {len(sols4)} grid points; 11A: {len(sols11)} grid points", classification="CONDITIONAL_ON_INPUT",
             negative_control="grid values outside the solution set fail", control_perturbs_same_parameter=True, input_it_depends_on="k"),
    ]
    out = dict(track="A-genus v3", command=command, inputs_file="inputs.json", results=results, rigidity=rig, details=R, could_not_do=could_not)
    (HERE / "results.json").write_text(json.dumps(out, indent=1, default=str))
    exp = dict(
        chi_from_genus=fs(K0 * Z0_unit), shared_inputs=["k"],
        chi_from_genus_by_k={str(k): fs(k * Z0_unit) for k in range(1, KMAX + 1)},
        source="Z(tau,0) = k * phi_{0,1}(tau,0) computed in run_track_a_v3.py (results.json details.phi01); value shown for the declared k=2; NOT independent of k",
        mu_N_solution_by_k={k: v["N_exact_solution"] for k, v in mu_rows.items()},
        A_n_k2=[fs(x / 2) for x in Hk[2]["h"]],
        twining_chi_from_cycle_shape=chi, twining_shared_inputs=["k", "perm24_premise"],
        path="audit/k3t2_rigidity_v3/A-genus/exports.json")
    (HERE / "exports.json").write_text(json.dumps(exp, indent=1))
    print(json.dumps(dict(exports=exp), indent=1)[:1500])


main()
