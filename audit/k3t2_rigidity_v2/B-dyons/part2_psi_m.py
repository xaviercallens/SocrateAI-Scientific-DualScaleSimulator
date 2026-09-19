"""
Part (2): psi_m via Delta*psi_m = G_{m+1}/A, for m = -1, 0, 1, 2.

DMZ tier-L ansatz (Dabholkar-Murthy-Zagier, arXiv:1208.4074v2, eq. (5.16), p.38 --
located by fetching the arXiv PDF and grepping for "5.16", then VISUALLY CONFIRMED by
rendering page 38 to a PNG (pdftoppm) and reading it, because plain pdftotext garbled the
m=2 line's exponents (it printed "48E4 AB2 + 10E6 A2B", which the rendered image shows is
WRONG -- the true typeset formula has no exponent on B in the second term and no B at all
in the third term; this was caught precisely BECAUSE the m=2 identity failed to close
under the pdftotext reading while the m=1 and m=3 lines -- independently cross-checked,
m=3 against the paper's own worked eq.(5.17) example -- closed exactly, so a mismatch on
m=2 alone flagged the transcription error rather than a real discrepancy):

    Delta psi_{-1} = A^{-1}
    Delta psi_0    = 2 A^{-1} B
    4 Delta psi_1  = 9 A^{-1} B^2 + 3 E4 A
    27 Delta psi_2 = 50 A^{-1} B^3 + 48 E4 A B + 10 E6 A^2

These four formulas are used ONLY as a literature-cited target to compare against G_k
computed HERE, independently, from the DMVV product (1.1)/(5.12)-(5.13) truncated in p,
expanded to order p^2 (for m=1) and p^3 (for m=2). They are tier-L (a cited identity, not
re-derived from first principles in this script) but the G_k series they are compared
against ARE derived here (not copied from the paper).

G_k construction (derived, not assumed):
  sum_k G_k p^k = prod_{r>=1} prod_{s>=0,t} (1 - p^r q^s y^t)^{-c(4rs-t^2)}
  Group factors by r: log(prod_r F_r(p)) = sum_r sum_{k>=1} P_k^{(r)} p^{rk} / k,
  where P_k^{(r)}(q,y) := sum_{s,t} c(4rs-t^2) q^{ks} y^{kt}  (the "A_r" power-sum series,
  evaluated at (q^k,y^k)). This is the standard prod(1-x_i)^{-c_i} = exp(sum c_i x_i^k p^k/k)
  identity, exact order by order in p. G_k is then recovered from L(p):=sum p^n L_n via the
  standard exp-of-power-series recurrence n*G_n = sum_{j=1}^n j*L_j*G_{n-j}, G_0=1 -- computed
  generically in code below (NOT hand-derived per k, to avoid algebra-transcription error).
  G_1 = 2B is recovered automatically (P_1^{(1)} = sum c(4s-t^2)q^sy^t = 2B, by definition of
  how c(D) was extracted from 2B -- flagged, as in v1, as carrying no independent content for
  the m=0 check).

Run: python part2_psi_m.py
Reads theta_forms_cache.json (written by `python theta_forms.py 16`) and
hurwitz_class_numbers_results.json in this directory (only for provenance echo; not used
in the m<=2 checks below, which need no class numbers).
Writes part2_psi_m_results.json.
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

QMAX = cache["QMAX"]
A = load2d("A_series")
B = load2d("B_series")
E4 = load1d("E4_series")
E6 = load1d("E6_series")
Delta = load1d("Delta_series")
cD = {int(k): pf(v) for k, v in cache["cD_table"].items()}

QCHK = 10  # truncation order for the identity checks. IMPORTANT MARGIN NOTE (found empirically,
# not assumed): a first attempt with theta cache QMAX=16 (only 2x QCHK=8) gave a SPURIOUS
# mismatch in the m=2 identity for n=6,7,8 while n=0..5 matched exactly with the unique
# forced integer coefficients (50,48,10) (over-determined solve, checked). Rebuilding the
# theta cache at QMAX=24 (`python theta_forms.py 24`) made the SAME QCHK=8 check pass at
# ALL n=0..8 -- i.e. QMAX=16 was an insufficient truncation margin for this cubic-degree
# identity, not a formula error. We use QMAX=40 here (4x QCHK=10) and separately verify
# stability against an independent QMAX=60 rebuild in margin_check_QMAX60/ (run
# `python margin_check_QMAX60/theta_forms.py 60` inside that dir first, then
# `python margin_stability_check.py` in this dir) -- see that script's output for the
# explicit QMAX=40-vs-60 agreement check, not assumed.

def to_qtuple(d1):
    return {(k, 0): v for k, v in d1.items()}

E4_2d = to_qtuple(E4)
E6_2d = to_qtuple(E6)

def A_r(r, qchk):
    """P_1^{(r)} = sum_{s>=0,t} c(4rs-t^2) q^s y^t, truncated s<=qchk."""
    out = {}
    for s in range(0, qchk + 1):
        tmax = int((4 * r * s + 1) ** 0.5) + 2
        for t in range(-tmax, tmax + 1):
            D = 4 * r * s - t * t
            if D in cD:
                out[(s, t)] = out.get((s, t), Fr(0)) + cD[D]
    return {k: v for k, v in out.items() if v != 0}

def subst_qk_yk(series, k, qchk):
    """Substitute q->q^k, y->y^k (i.e. evaluate the SAME series at (q^k,y^k))."""
    out = {}
    for (n, l), v in series.items():
        if k * n > qchk:
            continue
        out[(k * n, k * l)] = out.get((k * n, k * l), Fr(0)) + v
    return out

KMAX_P = 3  # need G_1, G_2, G_3

# ---- build L_n := coefficient of p^n in log(prod_r F_r(p)), for n=1..KMAX_P
L = {n: {} for n in range(1, KMAX_P + 1)}
for r in range(1, KMAX_P + 1):
    Ar = A_r(r, QCHK)
    for k in range(1, KMAX_P // r + 1):
        pk = r * k
        term = scal(Fr(1, k), subst_qk_yk(Ar, k, QCHK))
        L[pk] = add(L[pk], term)

# ---- exp-of-power-series recurrence: n*G_n = sum_{j=1}^n j*L_j*G_{n-j}, G_0 = 1
G = {0: {(0, 0): Fr(1)}}
for n in range(1, KMAX_P + 1):
    acc = {}
    for j in range(1, n + 1):
        Lj = L.get(j, {})
        Gnj = G[n - j]
        term = scal(j, mul(Lj, Gnj, QCHK, None))
        acc = add(acc, term)
    G[n] = scal(Fr(1, n), acc)

G1, G2, G3 = G[1], G[2], G[3]

# sanity: G_1 should equal 2B exactly (definitional, as flagged)
twoB = scal(2, B)
G1_eq_2B = (G1 == {k: v for k, v in twoB.items() if k[0] <= QCHK})

# All A^{-1} terms below are handled by multiplying the whole identity through by A first
# (this series engine has no division); m=1: 4*Delta*psi_1 = 9*B^2/A + 3*E4*A  <=>  (mult by A)  4*G_2 = 9*B^2 + 3*E4*A^2
B2 = mul(B, B, QCHK, None)
A2 = mul(A, A, QCHK, None)
lhs1 = scal(4, G2)
rhs1 = add(scal(9, B2), scal(3, mul(E4_2d, A2, QCHK, None)))
lhs1_t = {k: v for k, v in lhs1.items() if k[0] <= QCHK}
rhs1_t = {k: v for k, v in rhs1.items() if k[0] <= QCHK}
mism1 = {str(k): [str(lhs1_t.get(k, Fr(0))), str(rhs1_t.get(k, Fr(0)))]
         for k in set(lhs1_t) | set(rhs1_t) if lhs1_t.get(k, Fr(0)) != rhs1_t.get(k, Fr(0))}
m1_holds = (len(mism1) == 0)

# m=2: 27*Delta*psi_2 = 50*B^3/A + 48*E4*A*B + 10*E6*A^2  <=>  (mult by A)
#      27*G_3 = 50*B^3 + 48*E4*A^2*B + 10*E6*A^3
B3 = mul(B2, B, QCHK, None)
A3 = mul(A2, A, QCHK, None)
lhs2 = scal(27, G3)
rhs2 = add(scal(50, B3),
           scal(48, mul(E4_2d, mul(A2, B, QCHK, None), QCHK, None)),
           scal(10, mul(E6_2d, A3, QCHK, None)))
lhs2_t = {k: v for k, v in lhs2.items() if k[0] <= QCHK}
rhs2_t = {k: v for k, v in rhs2.items() if k[0] <= QCHK}
mism2 = {str(k): [str(lhs2_t.get(k, Fr(0))), str(rhs2_t.get(k, Fr(0)))]
         for k in set(lhs2_t) | set(rhs2_t) if lhs2_t.get(k, Fr(0)) != rhs2_t.get(k, Fr(0))}
m2_holds = (len(mism2) == 0)

result = {
    "QCHK": QCHK, "QMAX_theta_cache": QMAX,
    "dmz_5_16_source": "arXiv:1208.4074v2 (Dabholkar-Murthy-Zagier), eq. (5.16), p.36 -- fetched from https://people.mpim-bonn.mpg.de/zagier/files/arxiv/1208.4074/1208.4074v2.pdf and located by text search for '5.16'. Quoted formulas used as tier-L ansatz only.",
    "G1_equals_2B_exactly": G1_eq_2B,
    "m_minus1": {
        "statement": "Delta*psi_{-1} = A^{-1} (G_0=1 since Sym^0(K3)=point)",
        "status": "definitional (G_0 by construction of the DMVV product); no further identity to test",
    },
    "m_0": {
        "statement": "Delta*psi_0 = 2*A^{-1}*B  <=>  G_1 = 2B",
        "status": "holds identically BY CONSTRUCTION: c(D) is defined as the Fourier coefficients "
                  "of 2B, and G_1 = P_1^{(1)} = sum c(4s-t^2)q^sy^t is 2B term-for-term. Not an "
                  "independent test (flagged, as in v1).",
    },
    "m_1": {
        "statement": "4*Delta*psi_1 = 9*A^{-1}*B^2 + 3*E4*A  <=>  4*G_2 = 9*B^2 + 3*E4*A^2 (A cleared)",
        "holds_to_order_q^QCHK": m1_holds,
        "mismatches": mism1,
    },
    "m_2": {
        "statement": "27*Delta*psi_2 = 50*A^{-1}*B^3 + 48*E4*A*B + 10*E6*A^2  <=>  "
                     "27*G_3 = 50*B^3 + 48*E4*A^2*B + 10*E6*A^3 (A cleared)",
        "G3_derivation": "G_3 built from the exp-of-power-series recurrence above, from A_r(1,QCHK), "
                          "A_r(2,QCHK), A_r(3,QCHK) (r-graded pieces of the DMVV product), NOT copied "
                          "from the paper.",
        "holds_to_order_q^QCHK": m2_holds,
        "mismatches": mism2,
        "sample_G3_terms_n_leq_3": {str(k): str(v) for k, v in sorted(G3.items()) if k[0] <= 3},
    },
}

# ---- negative control: JOINT perturbation of (c(3), c(4)) by (delta3,delta4) in {-2..2}^2,
# re-run BOTH the m=1 and m=2 identity checks under each joint perturbation.
def rebuild_with_perturbed_cD(delta3, delta4, qchk):
    cD_pert = dict(cD)
    if delta3:
        cD_pert[3] = cD_pert.get(3, Fr(0)) + delta3
    if delta4:
        cD_pert[4] = cD_pert.get(4, Fr(0)) + delta4
    B_pert = {}
    for (n, l), v in B.items():
        D = 4 * n - l * l
        B_pert[(n, l)] = cD_pert[D] / 2 if D in cD_pert else v

    def Ar_pert(r, qchk_):
        out = {}
        for s in range(0, qchk_ + 1):
            tmax = int((4 * r * s + 1) ** 0.5) + 2
            for t in range(-tmax, tmax + 1):
                D = 4 * r * s - t * t
                if D in cD_pert:
                    out[(s, t)] = out.get((s, t), Fr(0)) + cD_pert[D]
        return {k: v for k, v in out.items() if v != 0}

    Lp = {n: {} for n in range(1, KMAX_P + 1)}
    for r in range(1, KMAX_P + 1):
        Arp = Ar_pert(r, qchk)
        for k in range(1, KMAX_P // r + 1):
            pk = r * k
            term = scal(Fr(1, k), subst_qk_yk(Arp, k, qchk))
            Lp[pk] = add(Lp[pk], term)
    Gp = {0: {(0, 0): Fr(1)}}
    for n in range(1, KMAX_P + 1):
        acc = {}
        for j in range(1, n + 1):
            term = scal(j, mul(Lp.get(j, {}), Gp[n - j], qchk, None))
            acc = add(acc, term)
        Gp[n] = scal(Fr(1, n), acc)
    G2p, G3p = Gp[2], Gp[3]
    B2p = mul(B_pert, B_pert, qchk, None)
    A2q = mul(A, A, qchk, None)
    lhs1p = scal(4, G2p)
    rhs1p = add(scal(9, B2p), scal(3, mul(E4_2d, A2q, qchk, None)))
    diffs1 = sum(1 for k in set(lhs1p) | set(rhs1p)
                 if lhs1p.get(k, Fr(0)) != rhs1p.get(k, Fr(0)))
    B3p = mul(B2p, B_pert, qchk, None)
    A3q = mul(A2q, A, qchk, None)
    lhs2p = scal(27, G3p)
    rhs2p = add(scal(50, B3p), scal(48, mul(E4_2d, mul(A2q, B_pert, qchk, None), qchk, None)),
                scal(10, mul(E6_2d, A3q, qchk, None)))
    diffs2 = sum(1 for k in set(lhs2p) | set(rhs2p)
                 if lhs2p.get(k, Fr(0)) != rhs2p.get(k, Fr(0)))
    return diffs1, diffs2

QCHK_PERT = 5  # smaller window to keep the 5x5 joint perturbation grid fast
joint_grid = {}
any_survives_m1 = []
any_survives_m2 = []
for d3 in range(-2, 3):
    for d4 in range(-2, 3):
        diffs1, diffs2 = rebuild_with_perturbed_cD(Fr(d3), Fr(d4), QCHK_PERT)
        key = f"({d3},{d4})"
        joint_grid[key] = {"m1_num_mismatches": diffs1, "m1_breaks": diffs1 > 0,
                            "m2_num_mismatches": diffs2, "m2_breaks": diffs2 > 0}
        if diffs1 == 0:
            any_survives_m1.append(key)
        if diffs2 == 0:
            any_survives_m2.append(key)

result["negative_control_joint_perturb_c3_c4"] = {
    "QCHK_used_for_this_grid": QCHK_PERT,
    "grid": "delta3,delta4 in {-2,-1,0,1,2} x {-2,-1,0,1,2}, 25 joint combinations",
    "results": joint_grid,
    "combinations_where_m1_identity_survives (expect only (0,0))": any_survives_m1,
    "combinations_where_m2_identity_survives (expect only (0,0))": any_survives_m2,
    "purpose": "Demonstrates the m=1 and m=2 DMZ (5.16) identities are genuine, "
               "non-tautological consistency checks: the ONLY joint (delta3,delta4) that "
               "preserves either identity is (0,0), i.e. no 1- or multi-unit deformation of "
               "c(3) or c(4), singly or jointly, survives.",
}

with open("part2_psi_m_results.json", "w") as f:
    json.dump(result, f, indent=1)
print(json.dumps({"G1_eq_2B": G1_eq_2B, "m1_holds": m1_holds, "m2_holds": m2_holds,
                   "m1_survivors": any_survives_m1, "m2_survivors": any_survives_m2}, indent=1))
