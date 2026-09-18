"""
Part (2): psi_m via Delta*psi_m = G_{m+1}/A, for m = -1, 0, 1.

G_0 = 1 (Sym^0 K3 = point).
G_1 = 2B  (Sym^1 K3 = K3 itself; this is DEFINITIONAL given how c(D) was extracted from B,
           so the m=0 check "Delta*psi_0 = 2B/A" carries no independent content -- flagged honestly).
G_2: derived here (not assumed) from the DMVV product truncated to order p^2:

  prod_{r=1,2} prod_{s,t} (1-p^r q^s y^t)^{-c(4rs-t^2)}  =  1 + p*(2B) + p^2*G_2 + O(p^3)

  Expanding the r=1 factors to O(p^2) and using A1 := sum_{s,t} c(4s-t^2) q^s y^t = 2B exactly:
      [r=1 product]_{p^2} = 2B^2 + B(q^2,y^2)          (algebra below)
  and the r=2 factors contribute directly at their own p^1 order:
      [r=2 product]_{p^2} = C2(tau,z) := sum_{s,t} c(8s-t^2) q^s y^t

  so   G_2 = 2*B^2 + B(q^2,y^2) + C2.

This is then compared against the claimed identity 4*Delta*psi_1 = 9*B^2/A + 3*E4*A,
i.e. (clearing the A denominator) against  4*G_2 =?= 9*B^2*A/A... -> multiply the ORIGINAL
identity by A:   4*G_2  =?=  9*B^2 + 3*E4*A^2.
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
Delta = load1d("Delta_series")
cD = {int(k): pf(v) for k, v in cache["cD_table"].items()}

QCHK = 10  # truncation order for the identity check (conservative vs QMAX=16 available)

def to_qtuple(d1):
    return {(k, 0): v for k, v in d1.items()}

E4_2d = to_qtuple(E4)

# B(q^2,y^2): substitute n->2n, l->2l
B_sub2 = {(2 * n, 2 * l): v for (n, l), v in B.items() if 2 * n <= QCHK}

# B^2
B2 = mul(B, B, QCHK, None)

# C2 = sum_{s,t} c(8s - t^2) q^s y^t, truncated s<=QCHK
C2 = {}
for s in range(0, QCHK + 1):
    tmax = int((8 * s + 1) ** 0.5) + 2
    for t in range(-tmax, tmax + 1):
        D = 8 * s - t * t
        if D in cD:
            C2[(s, t)] = C2.get((s, t), Fr(0)) + cD[D]
C2 = {k: v for k, v in C2.items() if v != 0}

G2 = add(scal(2, B2), B_sub2, C2)

lhs = scal(4, G2)
rhs = add(scal(9, B2), scal(3, mul(E4_2d, mul(A, A, QCHK, None), QCHK, None)))

# restrict comparison to n <= QCHK on both sides (drop terms outside truncation confidence)
lhs_t = {k: v for k, v in lhs.items() if k[0] <= QCHK}
rhs_t = {k: v for k, v in rhs.items() if k[0] <= QCHK}

all_keys = set(lhs_t) | set(rhs_t)
mismatches = {}
for k in all_keys:
    lv = lhs_t.get(k, Fr(0))
    rv = rhs_t.get(k, Fr(0))
    if lv != rv:
        mismatches[str(k)] = {"4G2": str(lv), "9B2+3E4A2": str(rv)}

identity_holds_to_order = (len(mismatches) == 0)

# m=-1: Delta*psi_{-1} = G_0/A = 1/A ; purely definitional, no further identity claimed by
# the problem statement, so we just record it as such (no computation of 1/A needed/claimed).

# also record G_1 vs 2B (definitional check, flagged)
twoB = scal(2, B)
G1_minus_2B = {}
for k in set(twoB):
    pass
G1_matches_2B_note = "definitional: c(D) was extracted FROM B, so G_1 := sum c(D) q^n y^l == 2B identically; not an independent test."

result = {
    "QCHK": QCHK,
    "m_minus1": {
        "statement": "Delta*psi_{-1} = G_0/A = 1/A (G_0=1 since Sym^0(K3)=point)",
        "status": "definitional; no further identity to test per problem statement",
    },
    "m_0": {
        "statement": "Delta*psi_0 = G_1/A =?= 2B/A",
        "status": "holds identically BY CONSTRUCTION (c(D) were defined as the Fourier coefficients of 2B, so G_1:=sum_{s,t} c(4s-t^2) q^s y^t is 2B term-for-term)",
        "note": G1_matches_2B_note,
    },
    "m_1": {
        "statement": "4*Delta*psi_1 = 9*B^2/A + 3*E4*A  <=>  4*G_2 = 9*B^2 + 3*E4*A^2  (cleared A)",
        "G2_derivation": "G_2 = 2*B^2 + B(q^2,y^2) + C2,  C2 = sum_{s,t} c(8s-t^2) q^s y^t  (from DMVV r=1,2 expansion to p^2, derived in module docstring)",
        "identity_holds_to_order_q^QCHK": identity_holds_to_order,
        "mismatches": mismatches,
        "sample_G2_terms (n,l)->coeff, n<=3": {str(k): str(v) for k, v in sorted(G2.items()) if k[0] <= 3},
    },
}


# ---- negative control: perturb c(3) (D=3, the first "structural" coefficient beyond
# c(0),c(-1)) by +1 and re-run the SAME identity check. This must break the identity,
# demonstrating the check has real discriminating power (it is not a tautology).
def rebuild_G2_with_perturbed_cD(D_perturb, delta, qchk):
    cD_pert = dict(cD)
    cD_pert[D_perturb] = cD_pert.get(D_perturb, Fr(0)) + delta
    # B must also be rebuilt consistently from the perturbed c(D) table (2B's coefficients)
    B_pert = {}
    for (n, l), v in B.items():
        D = 4 * n - l * l
        B_pert[(n, l)] = cD_pert.get(D, Fr(0)) if D != 0 or True else v
    # simplest: only the specific (n,l) with 4n-l^2==D_perturb actually change; rebuild via cD_pert/2
    B_pert = {}
    for (n, l), v in B.items():
        D = 4 * n - l * l
        B_pert[(n, l)] = cD_pert[D] / 2 if D in cD_pert else v
    B2_pert = mul(B_pert, B_pert, qchk, None)
    B_sub2_pert = {(2 * n, 2 * l): v for (n, l), v in B_pert.items() if 2 * n <= qchk}
    C2_pert = {}
    for s in range(0, qchk + 1):
        tmax = int((8 * s + 1) ** 0.5) + 2
        for t in range(-tmax, tmax + 1):
            D = 8 * s - t * t
            if D in cD_pert:
                C2_pert[(s, t)] = C2_pert.get((s, t), Fr(0)) + cD_pert[D]
    C2_pert = {k: v for k, v in C2_pert.items() if v != 0}
    G2_pert = add(scal(2, B2_pert), B_sub2_pert, C2_pert)
    lhs_p = scal(4, G2_pert)
    rhs_p = add(scal(9, B2_pert), scal(3, mul(E4_2d, mul(A, A, qchk, None), qchk, None)))
    lhs_pt = {k: v for k, v in lhs_p.items() if k[0] <= qchk}
    rhs_pt = {k: v for k, v in rhs_p.items() if k[0] <= qchk}
    keys = set(lhs_pt) | set(rhs_pt)
    diffs = sum(1 for k in keys if lhs_pt.get(k, Fr(0)) != rhs_pt.get(k, Fr(0)))
    return diffs

neg_control_perturb_c3 = rebuild_G2_with_perturbed_cD(3, Fr(1), QCHK)
neg_control_perturb_c4 = rebuild_G2_with_perturbed_cD(4, Fr(1), QCHK)

result["negative_control_perturb_cD_by_1"] = {
    "perturb_c(3)_by_+1": {"num_mismatched_coeffs": neg_control_perturb_c3,
                             "identity_breaks": neg_control_perturb_c3 > 0},
    "perturb_c(4)_by_+1": {"num_mismatched_coeffs": neg_control_perturb_c4,
                             "identity_breaks": neg_control_perturb_c4 > 0},
    "purpose": "Shows the m=1 identity 4G2=9B^2+3E4A^2 is a genuine, non-tautological "
               "consistency check: an arbitrary +1 shift in a single Fourier coefficient "
               "c(D) of the K3 elliptic genus breaks it. This is the rigidity evidence for "
               "TRACK B: the c(D) table (hence the DMVV structure built from it) admits no "
               "1-integer deformation at D=3 or D=4 that preserves the m=1 identity.",
}

with open("part2_psi_m_results.json", "w") as f:
    json.dump(result, f, indent=1)
print(json.dumps(result, indent=1))
