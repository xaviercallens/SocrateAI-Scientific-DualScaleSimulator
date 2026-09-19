"""
Part (3): m=1 "immortal" polar subtraction: G_2/A - N*A_{2,1}  =?=  3*E4*A - M*Hhat.

=== v1 defect #1 fixed here (general-m A_{2,m}, s<=-1 branch) ===
A_{2,m}(tau,z) := sum_{s in Z} q^{m s^2 + s} y^{2ms+1} / (1 - q^s y)^2, expanded in the
region |q| < |y| < 1 (the convention stated for A_{2,m}; matches the region used
throughout this track for all y-series).

Derivation (done here, not copied): for s>=1, 1/(1-q^s y)^2 = sum_{k>=0} (k+1) q^{sk} y^k
directly (|q^s y|<1 in this region), so the term is q^{m s^2+s+sk} y^{2ms+1+k}.
For s=0 the prefactor collapses to y/(1-y)^2 = sum_{k>=0}(k+1) y^{k+1}, independent of m.
For s<=-1, write s=-t (t>=1). Then 1-q^{-t}y = -q^{-t}y(1-q^t/y), so
    (1-q^{-t}y)^{-2} = q^{2t} y^{-2} (1-q^t/y)^{-2} = sum_{k>=0} (k+1) q^{2t+tk} y^{-2-k}
                      = sum_{k>=0} (k+1) q^{-2s-sk} y^{-2-k}     (t=-s)
which converges in this region since |q^t/y| = |q|^t/|y| < |q|/|y| < 1 for t>=1 given
|q|<|y|. Multiplying by the term's own prefactor q^{m s^2+s} y^{2ms+1} gives
    q^{m s^2+s} y^{2ms+1} * q^{-2s-sk} y^{-2-k} = q^{m s^2 - s - s k} y^{2ms - 1 - k}.
So for s<=-1:  n = m*s^2 - s - s*k,  l = 2*m*s - 1 - k   (coefficient (k+1)).
v1's bug: it used n = -2s-sk (only the (1-q^sy)^{-2}-expansion part), DROPPING the
q^{m s^2+s} prefactor's contribution to n entirely. For m=1 this coincides with the
correct n=s^2-s-sk only at s=-1 (both give n=2+k); they first diverge at s=-2 where v1
puts weight at n=4+2k (correct: n=6+2k). REGRESSION CHECK below shows this precisely:
before the fix, A_{2,1}[(4,-5)] = 4 (spurious s=-2,k=0 term at n=4); after the fix it is 3,
and the previously "unresolved" remainder of exactly -324 at (n,l)=(4,-5) (D=-9), reported
honestly as unresolved in v1's REPORT, disappears (remainder becomes 0 -- checked below).

=== v1 defect #2: undocumented args -- fixed: every QCHK/YCAP/QMAX stated + exact command ===
=== v1 defect #3: coordinate-wise scans -- fixed: true joint (N,M) 2D grid below ===
=== v1 defect #4: D>=-4 window that hid the bug -- REMOVED. Structural fact (proved from
already-established index-1 vanishing of A's D<=-2 coefficients and H(D)=0 for D<0, not
assumed): 3*E4*A and M*Hhat both vanish identically at every D<=-2, for ANY M, so the
remainder G_2/A - N*A_{2,1} must ALSO vanish at every D<=-2 for the correct N -- a genuine,
target-free (M-free) structural selecting condition on N, checked over the FULL window.

Run: python part3_polar.py
Reads theta_forms_cache.json (`python theta_forms.py 40` in this directory) and
hurwitz_class_numbers_results.json (`python hurwitz_class_numbers.py` in this directory).
Writes part3_polar_results.json.
"""
import json
from fractions import Fraction as Fr
from series import mul, add, scal

with open("theta_forms_cache.json") as f:
    cache = json.load(f)

def pf(s):
    if s is None:
        return None
    if "/" in s:
        a, b = s.split("/")
        return Fr(int(a), int(b))
    return Fr(int(s))

def load2d(key):
    out = {}
    for k, v in cache[key].items():
        n, l = k.split(",")
        out[(int(n), int(l))] = pf(v)
    return out

def load1d(key):
    return {int(k): pf(v) for k, v in cache[key].items()}

A = load2d("A_series")
B = load2d("B_series")
E4 = load1d("E4_series")

QCHK = 8    # truncation order for the polar-subtraction identity check
YCAP = 60   # y-window cap for the invA / A_{2,m} construction
QMAX_theta = cache["QMAX"]
assert QMAX_theta >= 3 * QCHK, f"theta cache QMAX={QMAX_theta} too small for QCHK={QCHK} margin (need >=3x, learned from part2's QMAX16-vs-24 finding)"

def slice_n(d2d, n):
    return {l: v for (nn, l), v in d2d.items() if nn == n}

def conv_y(s1, s2, ycap):
    out = {}
    for l1, v1 in s1.items():
        for l2, v2 in s2.items():
            l = l1 + l2
            if abs(l) > ycap:
                continue
            out[l] = out.get(l, Fr(0)) + v1 * v2
    return {k: v for k, v in out.items() if v != 0}

def to_2d(seq_by_n):
    out = {}
    for n, s in seq_by_n.items():
        for l, v in s.items():
            out[(n, l)] = v
    return out

# ---- Step 1: invA = 1/A (unchanged from v1: not flagged as defective)
A_by_n = {n: slice_n(A, n) for n in range(0, QCHK + 1)}
invA_by_n = {}
R0 = {m: Fr(m) for m in range(1, YCAP + 1)}  # A0^{-1} = sum m y^m
invA_by_n[0] = R0
A0 = A_by_n[0]
chk = conv_y(A0, R0, YCAP)
a0_r0_ok_low = all(chk.get(l, Fr(0)) == (Fr(1) if l == 0 else Fr(0)) for l in range(-5, 6))
for n in range(1, QCHK + 1):
    acc = {}
    for i in range(1, n + 1):
        term = conv_y(A_by_n[i], invA_by_n[n - i], YCAP)
        for l, v in term.items():
            acc[l] = acc.get(l, Fr(0)) + v
    invA_by_n[n] = {l: -v for l, v in conv_y(R0, acc, YCAP).items() if v != 0}
invA = to_2d(invA_by_n)

# ---- Step 2: G2/A = 9*B^2/(4A) + 3*E4*A/4 (algebraic identity from part2's G2 formula;
# re-derived here directly via invA rather than re-importing part2's G2, so this script is
# self-contained and independently runnable per the reproducibility requirement)
B2 = mul(B, B, QCHK, None)
NineB2 = scal(9, B2)
NineB2_over_4A = {}
for n in range(0, QCHK + 1):
    acc = {}
    for i in range(0, n + 1):
        term = conv_y(slice_n(NineB2, i), invA_by_n.get(n - i, {}), YCAP)
        for l, v in term.items():
            acc[l] = acc.get(l, Fr(0)) + v
    for l, v in acc.items():
        if v != 0:
            NineB2_over_4A[(n, l)] = NineB2_over_4A.get((n, l), Fr(0)) + v / 4
E4_2d = {(n, 0): v for n, v in E4.items() if n <= QCHK}
threeE4A_over4 = scal(Fr(3, 4), mul(E4_2d, A, QCHK, None))
G2_over_A = add(NineB2_over_4A, threeE4A_over4)

# ---- Step 3: A_{2,m}(qchk, ycap), general m -- FIXED s<=-1 branch (see docstring).
# Bound s directly by "m*s^2+s <= qchk" (s>=1) / "m*t^2+t <= qchk" (t=-s>=1), NOT by a
# guessed smax with a skip-if-too-large guard (an earlier version of this function used a
# fixed smax=int((qchk/m)**0.5)+3 combined with a `continue` on base_n>qchk, which for
# m=1,QCHK=8 happens to cover everything needed but is not a safe general-m bound and
# silently drops valid s for larger m -- fixed here to loop exactly while the bound holds).
def build_A2m(m, qchk, ycap):
    out = {}
    # s >= 1: n = m*s^2+s+s*k, l = 2ms+1+k
    s = 1
    while m * s * s + s <= qchk:
        base_n = m * s * s + s
        kmax = (qchk - base_n) // s
        for k in range(0, kmax + 1):
            n = base_n + s * k
            l = 2 * m * s + 1 + k
            if abs(l) <= ycap:
                out[(n, l)] = out.get((n, l), Fr(0)) + (k + 1)
        s += 1
    # s = 0: y/(1-y)^2 = sum_{l>=1} l*y^l, independent of m (base_n=0 always <= qchk)
    for l in range(1, ycap + 1):
        out[(0, l)] = out.get((0, l), Fr(0)) + l
    # s <= -1 (FIXED branch): n = m*s^2-s-s*k = m*t^2+t+t*k (t=-s>=1), l=2ms-1-k=-2mt-1-k
    t = 1
    while m * t * t + t <= qchk:
        base_n = m * t * t + t
        kmax = (qchk - base_n) // t
        for k in range(0, kmax + 1):
            n = base_n + t * k
            l = -2 * m * t - 1 - k
            if abs(l) <= ycap:
                out[(n, l)] = out.get((n, l), Fr(0)) + (k + 1)
        t += 1
    return {k: v for k, v in out.items() if v != 0}

A21 = build_A2m(1, QCHK, YCAP)

# sanity: for n>=1, the s>=1 and s<=-1 branches contribute symmetric (k+1) weight to
# (n,l) and (n,-l) respectively (same base_n formula under s<->t=-s, l<->-l), so A21
# should be symmetric under l->-l for n>=1 (checked explicitly; the n=0/s=0 branch alone
# is NOT symmetric, since y/(1-y)^2 only has positive powers of y -- confirmed below, not
# assumed, and excluded from the check for that documented reason).
_A21_symmetric_n_geq_1 = all(A21.get((n, l), Fr(0)) == A21.get((n, -l), Fr(0))
                              for (n, l) in list(A21) if n >= 1)

# ---- REGRESSION CHECK: the v1 bug, precisely (see docstring). Reconstruct v1's (buggy)
# s<=-1 branch verbatim for comparison, at the exact point the bug first bites (s=-2,k=0).
def build_A21_v1_buggy(qchk, ycap):
    out = {}
    smax = int(qchk ** 0.5) + 3
    for s in range(1, smax + 1):
        kmax = (qchk - s * s - s) // s if s > 0 else 0
        for k in range(0, max(kmax, -1) + 1):
            n = s * s + s + s * k
            if n > qchk or n < 0:
                continue
            l = 2 * s + 1 + k
            if abs(l) > ycap:
                continue
            out[(n, l)] = out.get((n, l), Fr(0)) + (k + 1)
    for m_ in range(1, ycap + 1):
        out[(0, m_)] = out.get((0, m_), Fr(0)) + m_
    for s in range(-smax, 0):
        for k in range(0, qchk + 1):
            n = -2 * s - s * k  # v1's bug: missing the q^{s^2+s} prefactor contribution
            if n > qchk or n < 0:
                continue
            l0 = 2 * s + 1
            l = l0 - 2 - k
            if abs(l) > ycap:
                continue
            out[(n, l)] = out.get((n, l), Fr(0)) + (k + 1)
    return {k: v for k, v in out.items() if v != 0}

A21_v1_buggy = build_A21_v1_buggy(QCHK, YCAP)
regression_check = {
    "A21_fixed_at_(4,-5)": str(A21.get((4, -5), Fr(0))),
    "A21_v1_buggy_at_(4,-5)": str(A21_v1_buggy.get((4, -5), Fr(0))),
    "expected": "fixed=3, v1_buggy=4 (v1 double-counted a spurious s=-2,k=0 term at n=4 "
                "instead of the correct n=6)",
}

# ---- Step 4: Hhat from (extended, whole-table-validated) Hurwitz numbers
with open("hurwitz_class_numbers_results.json") as f:
    Hdata_full = json.load(f)
Hdata = Hdata_full["H_table"]

def Hval(D):
    if D < 0:
        return Fr(0)
    if D == 0:
        return Fr(-1, 12)
    if D % 4 not in (0, 3):
        return Fr(0)
    if str(D) in Hdata:
        return pf(Hdata[str(D)]["H"])
    raise ValueError(f"H({D}) not in tabulated range (extend hurwitz_class_numbers.py DMAX)")

Hhat = {}
for n in range(0, QCHK + 1):
    lmax = int((4 * n + 1) ** 0.5) + 2 if n > 0 else 1
    for l in range(-lmax, lmax + 1):
        D = 4 * n - l * l
        if D > 40:
            continue  # outside tabulated Hurwitz range; excluded from window below anyway
        h = Hval(D)
        if h != 0:
            Hhat[(n, l)] = h

# ---- Step 5: remainder and JOINT (N,M) scan, full window (NO D>=-4 restriction)
def remainder(N, M):
    Npart = scal(N, A21)
    rem = {}
    for k in set(G2_over_A) | set(Npart):
        v = G2_over_A.get(k, Fr(0)) - Npart.get(k, Fr(0))
        if v != 0:
            rem[k] = v
    target = add(scal(3, mul(E4_2d, A, QCHK, None)), scal(-M, Hhat))
    return rem, target

def compare(N, M, qchk, drop_top=2, keep_l=20, keep_D_leq=None):
    rem, target = remainder(N, M)
    keys = set(rem) | set(target)
    keys = {k for k in keys if k[0] <= qchk - drop_top and abs(k[1]) <= keep_l}
    if keep_D_leq is not None:
        keys = {k for k in keys if (4 * k[0] - k[1] * k[1]) <= keep_D_leq}
    mism = {k: (rem.get(k, Fr(0)), target.get(k, Fr(0))) for k in keys
            if rem.get(k, Fr(0)) != target.get(k, Fr(0))}
    return mism, len(keys)

N_window = list(range(310, 340))
M_window = list(range(620, 680, 2))
joint_grid = {}
solutions = []
for N in N_window:
    for M in M_window:
        mism, nkeys = compare(N, M, QCHK)
        if len(mism) == 0:
            solutions.append((N, M))
            joint_grid[f"({N},{M})"] = {"matches": True, "num_keys_checked": nkeys}

# ---- structural D<=-2 test (target-free: 3E4A-M*Hhat is IDENTICALLY zero there for any M;
# checked explicitly it's zero in our own data too, not just asserted) -- run for each
# candidate N found above, over the FULL window (this is what v1's D>=-4 filter hid)
def d_leq_m2_test(N, qchk):
    rem, _ = remainder(N, 0)  # target is 0 at D<=-2 regardless of M; use M=0 for the target arg
    bad = {k: v for k, v in rem.items() if (4 * k[0] - k[1] * k[1]) <= -2 and k[0] <= qchk - 2 and abs(k[1]) <= 20}
    # also confirm target really is zero there in our own tables (not assumed)
    target_check_M0 = add(scal(3, mul(E4_2d, A, qchk, None)), scal(0, Hhat))
    target_nonzero_at_Dleqm2 = {k: v for k, v in target_check_M0.items()
                                 if (4 * k[0] - k[1] * k[1]) <= -2 and v != 0}
    return bad, target_nonzero_at_Dleqm2

d_leq_m2_results = {}
for N in sorted(set(N for N, M in solutions)):
    bad, tnz = d_leq_m2_test(N, QCHK)
    d_leq_m2_results[N] = {"remainder_nonzero_at_D_leq_-2": {str(k): str(v) for k, v in bad.items()},
                            "target_3E4A_nonzero_at_D_leq_-2 (should be empty, M-independent)":
                                {str(k): str(v) for k, v in tnz.items()},
                            "structural_D_leq_-2_test_passes": (len(bad) == 0)}

# ---- q^0 tail-slope structural check for N (second, independent selector; advisor-derived)
# At q^0: 1/A = sum_{m>=1} m*y^m (R0 above), so G2/A|_{q0} = (9B^2/(4A)+3E4A/4)|_{q0} has a
# linear-in-l tail for large l; A_{2,1}|_{q0} = sum_{l>=1} l*y^l (same R0, since the s=0
# branch of A_{2,m} is m-independent). N is fixed structurally (not from ANY target) by
# requiring the l-> l+1 tail slope of G2_over_A|_{q0} - N*A21|_{q0} to vanish (bounded
# support) rather than grow linearly. Computed from OUR OWN series, not hardcoded.
G2A_q0 = {l: v for (n, l), v in G2_over_A.items() if n == 0}
A21_q0 = {l: v for (n, l), v in A21.items() if n == 0}
tail_ls = [l for l in range(10, 30) if l in G2A_q0]
slopes_G2A = {l: (G2A_q0.get(l + 1, Fr(0)) - G2A_q0.get(l, Fr(0))) for l in tail_ls[:-1]}
slopes_A21 = {l: (A21_q0.get(l + 1, Fr(0)) - A21_q0.get(l, Fr(0))) for l in tail_ls[:-1]}
# both slopes should individually be constant (=: sG2A, sA21) for l large enough (bounded
# support of G2_over_A - N*A21 means slopes match once we pick N = sG2A/sA21)
slope_G2A_vals = set(slopes_G2A.values())
slope_A21_vals = set(slopes_A21.values())
N_from_slope = None
if len(slope_G2A_vals) == 1 and len(slope_A21_vals) == 1:
    sg = next(iter(slope_G2A_vals))
    sa = next(iter(slope_A21_vals))
    if sa != 0:
        N_from_slope = sg / sa

tail_slope_check = {
    "slopes_G2_over_A_q0 (l -> l+1 diffs, should be constant)": {str(k): str(v) for k, v in slopes_G2A.items()},
    "slopes_A21_q0": {str(k): str(v) for k, v in slopes_A21.items()},
    "N_from_tail_slope_ratio": str(N_from_slope),
    "note": "Structural, target-free determination of N from the q^0-slice linear tail "
            "alone; compared below against the (N,M)-grid solution set.",
}

# ---- negative control (a): M=0 (class numbers not used) must fail everywhere
if solutions:
    N0 = solutions[0][0]
    mism_M0, nkeys_M0 = compare(N0, 0, QCHK)
    M0_control = {"N_used": N0, "num_mismatches_with_M=0": len(mism_M0),
                  "num_keys_checked": nkeys_M0, "fails_as_expected": len(mism_M0) > 0}
else:
    M0_control = {"note": "no (N,M) solution found in the grid; M=0 control skipped"}

# ---- negative control (b): joint (delta_c3, delta_c4) perturbation of c(D), rebuilding A,
# B, invA, A21, G2/A, Hhat all consistently, then checking whether ANY (N,M) near the
# unperturbed solution still works. Kept at a smaller QCHK for speed (25 joint runs).
def rebuild_and_check(delta3, delta4, qchk, ycap, N_try, M_try):
    cD_local = {}
    for k, v in cache["cD_table"].items():
        cD_local[int(k)] = pf(v)
    if delta3:
        cD_local[3] = cD_local.get(3, Fr(0)) + delta3
    if delta4:
        cD_local[4] = cD_local.get(4, Fr(0)) + delta4
    B_pert = {}
    for (n, l), v in B.items():
        D = 4 * n - l * l
        B_pert[(n, l)] = cD_local[D] / 2 if D in cD_local else v
    A_local = A  # A's low-D coefficients (D<=... ) are unaffected by perturbing c(3),c(4)
    # (A = phi_{-2,1} is independent of B's c(D) table by construction; only B, and hence
    # anything built FROM B, changes under this perturbation)
    A_by_n_l = {n: slice_n(A_local, n) for n in range(0, qchk + 1)}
    invA_l = {0: R0}
    for n in range(1, qchk + 1):
        acc = {}
        for i in range(1, n + 1):
            term = conv_y(A_by_n_l[i], invA_l[n - i], ycap)
            for l, v in term.items():
                acc[l] = acc.get(l, Fr(0)) + v
        invA_l[n] = {l: -v for l, v in conv_y(R0, acc, ycap).items() if v != 0}
    B2p = mul(B_pert, B_pert, qchk, None)
    NineB2p_over_4A = {}
    for n in range(0, qchk + 1):
        acc = {}
        for i in range(0, n + 1):
            term = conv_y(slice_n(scal(9, B2p), i), invA_l.get(n - i, {}), ycap)
            for l, v in term.items():
                acc[l] = acc.get(l, Fr(0)) + v
        for l, v in acc.items():
            if v != 0:
                NineB2p_over_4A[(n, l)] = NineB2p_over_4A.get((n, l), Fr(0)) + v / 4
    G2A_p = add(NineB2p_over_4A, threeE4A_over4)
    Npart = scal(N_try, A21)
    rem = {}
    for k in set(G2A_p) | set(Npart):
        v = G2A_p.get(k, Fr(0)) - Npart.get(k, Fr(0))
        if v != 0:
            rem[k] = v
    target = add(scal(3, mul(E4_2d, A_local, qchk, None)), scal(-M_try, Hhat))
    keys = set(rem) | set(target)
    keys = {k for k in keys if k[0] <= qchk - 2 and abs(k[1]) <= 15}
    mism = sum(1 for k in keys if rem.get(k, Fr(0)) != target.get(k, Fr(0)))
    return mism

QCHK_PERT = QCHK  # MUST match QCHK: an earlier version used QCHK_PERT=5 while
# threeE4A_over4 (used inside rebuild_and_check) was the QCHK=8 module-level value -- a
# truncation-order mismatch that silently made the perturbation test meaningless (caught
# by the advisor review, not self-discovered). Fixed by using the SAME qchk throughout, so
# threeE4A_over4, A21, and Hhat (all correctly built once at QCHK -- they do NOT depend on
# c(3),c(4), only B does, so they are legitimately reused unperturbed) are on the same
# truncation order as the perturbed invA_l/G2A_p built inside the function.
N_center, M_center = (solutions[0] if solutions else (None, None))
# sanity: the zero-perturbation point must reproduce the unperturbed (N,M) match exactly
sanity_zero_pert_mismatches = (
    rebuild_and_check(Fr(0), Fr(0), QCHK_PERT, YCAP, N_center, M_center)
    if N_center is not None else None
)
joint_pert_results = {}
survivors = []
if N_center is not None:
    for d3 in range(-2, 3):
        for d4 in range(-2, 3):
            mism = rebuild_and_check(Fr(d3), Fr(d4), QCHK_PERT, YCAP, N_center, M_center)
            key = f"({d3},{d4})"
            joint_pert_results[key] = {"num_mismatches": mism, "breaks": mism > 0}
            if mism == 0:
                survivors.append(key)

result = {
    "QCHK": QCHK, "YCAP": YCAP, "QMAX_theta_cache": QMAX_theta,
    "margin_note": "QMAX=40 is 5x QCHK=8; see part2's QMAX16-vs-24 finding and "
                   "margin_stability_check.py for the empirical margin-adequacy evidence "
                   "this track relies on (same theta cache, same construction).",
    "invA0_times_A0_correct_near_low_l": a0_r0_ok_low,
    "regression_check_v1_bug": regression_check,
    "N_M_joint_grid": {
        "N_window": [N_window[0], N_window[-1]], "M_window": [M_window[0], M_window[-1], "step2"],
        "num_pairs_scanned": len(N_window) * len(M_window),
        "solutions_found (N,M) matching on drop_top2_keep_l20_no_D_filter": solutions,
    },
    "structural_D_leq_-2_test": d_leq_m2_results,
    "tail_slope_structural_N_check": tail_slope_check,
    "negative_control_M0": M0_control,
    "negative_control_joint_perturb_c3_c4": {
        "QCHK_used": QCHK_PERT,
        "note": "perturbing c(3),c(4) moves B (hence G2_over_A) only -- A, A_{2,1}(=A21) "
                "and Hhat do not depend on B's Fourier coefficients and are correctly "
                "reused unperturbed at the SAME QCHK; this control therefore tests "
                "whether the B-side of the identity is rigid, given the (independently "
                "constructed) A-side and Hhat fixed. sanity_zero_perturbation_reproduces_unperturbed_match "
                "confirms (0,0) recovers exactly the unperturbed result.",
        "sanity_zero_perturbation_mismatches (must be 0)": sanity_zero_pert_mismatches,
        "N_M_tested": [N_center, M_center],
        "grid": joint_pert_results,
        "survivors (expect only (0,0))": survivors,
    },
    "A21_symmetric_under_l_to_-l_for_n_geq_1": _A21_symmetric_n_geq_1,
}

with open("part3_polar_results.json", "w") as f:
    json.dump(result, f, indent=1)
print(json.dumps({
    "regression_check": regression_check,
    "solutions_found": solutions,
    "structural_D_leq_-2_all_pass": all(v["structural_D_leq_-2_test_passes"] for v in d_leq_m2_results.values()) if d_leq_m2_results else None,
    "N_from_tail_slope": str(N_from_slope),
    "M0_control_fails_as_expected": M0_control.get("fails_as_expected"),
    "joint_pert_survivors": survivors,
    "sanity_zero_pert_mismatches": sanity_zero_pert_mismatches,
    "A21_symmetric_n_geq_1": _A21_symmetric_n_geq_1,
}, indent=1))
