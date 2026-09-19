"""
Part (1): G_k at z=0 (y=1), q^0, for k=0..5 must equal e(Hilb^k(K3)), compared against
the coefficients of prod_n (1-p^n)^{-24} (Goettsche's formula for e(S)=24).

DMVV: sum_k G_k p^k = prod_{r>=1,s>=0,t in Z} (1 - p^r q^s y^t)^{-c(4rs-t^2)}.

Restriction to q^0, y=1 (derived, not assumed):
  Expanding (1-p^r q^s y^t)^{-c} = sum_j C(c+j-1,j) p^{rj} q^{sj} y^{tj}.
  A term contributes to q^0 only if sj=0, i.e. j=0 (trivial) or s=0.
  So the q^0 part of the whole product = prod_{r>=1} prod_{t} (1 - p^r y^t)^{-c(-t^2)}.
  c(-t^2) is nonzero (from the theta computation) only for t=0 (c(0)) and t=+-1 (c(-1)),
  since c(D)=0 for D<=-2 for a weak Jacobi form of this weight/index (checked below).
  At y=1: prod_r (1-p^r)^{-c(0)} (1-p^r)^{-c(-1)} (1-p^r)^{-c(-1)} = prod_r (1-p^r)^{-(c(0)+2c(-1))}.

This script computes that restricted product EXACTLY from the c(D) table (theta_forms_cache.json),
for k=0..5, and separately computes prod_n(1-p^n)^{-24} coefficients, and compares.
It also runs the kappa rigidity scan: exponent -kappa*c(...), kappa in a grid around 1.
"""
import json
from fractions import Fraction as Fr
from math import comb

with open("theta_forms_cache.json") as f:
    cache = json.load(f)

def parse_frac(s):
    if s is None:
        return None
    if "/" in s:
        a, b = s.split("/")
        return Fr(int(a), int(b))
    return Fr(int(s))

cD = {int(k): parse_frac(v) for k, v in cache["cD_table"].items()}
c0 = cD[0]
cm1 = cD[-1]

# Check no other D<=-2 contributes (should be absent/zero from the table; the table only
# has D=4n-l^2 for n>=0 with the (n,l) pairs that actually arose from the truncated theta
# series, so D=-4,-5,... simply never appear -- confirm explicitly none present with D<=-2.
d_leq_m2 = {D: v for D, v in cD.items() if D <= -2}

def euler_numbers_from_c0_cm1(kappa, kmax):
    """[p^k] of prod_r (1-p^r)^{-kappa*(c0+2*cm1)}  (the q^0,y=1 restriction of DMVV with
    exponent multiplier kappa), for k=0..kmax. kappa*(c0+2cm1) need not be an integer for
    general kappa; we still expand formally with binomial coefficients over Q."""
    exponent = kappa * (c0 + 2 * cm1)  # a Fraction, possibly non-integer for kappa!=1
    # prod_{r=1}^{kmax} (1-p^r)^{-exponent}, coefficients of p^0..p^kmax
    s = {0: Fr(1)}
    for r in range(1, kmax + 1):
        # (1-p^r)^{-exponent} = sum_j C(-exponent, j) (-1)^j p^{r j} ... use generalized binomial
        factor = {}
        jmax = kmax // r
        for j in range(0, jmax + 1):
            # C(-exponent, j) * (-1)^j = C(exponent+j-1, j)  [generalized binomial, exponent may be non-integer Fraction]
            c = Fr(1)
            for i in range(j):
                c *= (exponent + i)
            if j > 0:
                fact = 1
                for i in range(1, j + 1):
                    fact *= i
                c = c / fact
            factor[r * j] = c
        new_s = {}
        for e1, v1 in s.items():
            for e2, v2 in factor.items():
                e = e1 + e2
                if e > kmax:
                    continue
                new_s[e] = new_s.get(e, Fr(0)) + v1 * v2
        s = new_s
    return s

def gottsche_prod_1_minus_pn_neg24(kmax):
    s = {0: Fr(1)}
    for n in range(1, kmax + 1):
        factor = {}
        jmax = kmax // n
        for j in range(0, jmax + 1):
            factor[n * j] = Fr(comb(24 + j - 1, j))
        new_s = {}
        for e1, v1 in s.items():
            for e2, v2 in factor.items():
                e = e1 + e2
                if e > kmax:
                    continue
                new_s[e] = new_s.get(e, Fr(0)) + v1 * v2
        s = new_s
    return s

KMAX = 5
euler_kappa1 = euler_numbers_from_c0_cm1(Fr(1), KMAX)
gottsche = gottsche_prod_1_minus_pn_neg24(KMAX)

match_k0_5 = all(euler_kappa1[k] == gottsche[k] for k in range(KMAX + 1))

# rigidity scan on kappa (rational grid around 1)
kappa_grid = [Fr(n, 12) for n in range(0, 37)]  # 0, 1/12, ..., 3 in steps of 1/12
kappa_solutions = []
for kap in kappa_grid:
    eul = euler_numbers_from_c0_cm1(kap, KMAX)
    ok = all(eul.get(k) == gottsche[k] for k in range(KMAX + 1))
    if ok:
        kappa_solutions.append(str(kap))

# negative control: kappa = 1 + 1/24 and kappa = 1 - 1/24 must fail
neg_controls = {}
for kap in [Fr(1) + Fr(1, 24), Fr(1) - Fr(1, 24), Fr(23, 24), Fr(25, 24)]:
    eul = euler_numbers_from_c0_cm1(kap, KMAX)
    ok = all(eul.get(k) == gottsche[k] for k in range(KMAX + 1))
    neg_controls[str(kap)] = {"matches_all_k_0_5": ok,
                               "k1": str(eul.get(1)), "gottsche_k1": str(gottsche[1])}

result = {
    "c0": str(c0), "cm1": str(cm1),
    "c0_plus_2cm1": str(c0 + 2 * cm1),
    "D_leq_-2_present_in_table (should be empty)": {str(k): str(v) for k, v in d_leq_m2.items()},
    "euler_numbers_k0_5_from_DMVV_q0_y1_kappa1": {str(k): str(v) for k, v in sorted(euler_kappa1.items())},
    "euler_numbers_k0_5_from_gottsche_prod_1-pn_neg24": {str(k): str(v) for k, v in sorted(gottsche.items())},
    "match_for_k_0_to_5": match_k0_5,
    "kappa_rigidity_scan": {
        "grid": "n/12 for n=0..36 (0 to 3)",
        "solution_set (kappa matching ALL k=0..5)": kappa_solutions,
    },
    "negative_control_near_kappa_1": neg_controls,
    "honesty_note": (
        "This q^0,y=1 restriction of DMVV collapses ANALYTICALLY (derived above, not assumed) "
        "to prod_r (1-p^r)^{-kappa*(c0+2*cm1)}, i.e. a single linear equation kappa*(c0+2*cm1)=24 "
        "in kappa. It is tier B exact arithmetic but constrains only the ONE combination "
        "c(0)+2*c(-1), not the full DMVV structure. Do not oversell this as testing 'the whole "
        "product'; Part 2's G_2 identity (4*Delta*psi_1 = 9B^2/A + 3E4*A) is the structural test "
        "that uses c(D) for D>4 as well."
    ),
}

with open("part1_euler_numbers_results.json", "w") as f:
    json.dump(result, f, indent=1)
print(json.dumps(result, indent=1))
