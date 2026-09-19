"""
Hurwitz class numbers H(D), D >= 0, D = 0 or 3 (mod 4), by exact reduced-form counting.

METHOD (independent of any literature table): H(D) for D>0 is the number of SL(2,Z)-
equivalence classes of positive-definite integral binary quadratic forms (a,b,c) of
discriminant b^2-4ac = -D, summed over ALL discriminants -D*f^2 for f=1,2,... that also
have disc = -D (i.e. imprimitive forms included -- this is exactly the Hurwitz, as
opposed to Gauss/class-number-h, convention). A form is reduced iff
    -a < b <= a <= c,   and additionally b >= 0 whenever a==c or a==b.
Weight 1/2 is given to the (unique, when present) form equivalent to a*(x^2+y^2)
(b=0, a=c); weight 1/3 to the form equivalent to a*(x^2+xy+y^2) (a=b=c); weight 1
otherwise. H(0) := -1/12 (the standard convention for the D=0 case of Zagier's
non-holomorphic Eisenstein series; not obtained from the counting method, stated as
given -- exactly as in v1).

This is UNCHANGED from v1's method (v1's counting code had no defect identified by the
task). What is NEW here:
  (a) extended from D<=15 to D<=40 (needed for Track B's m=1 immortal check window);
  (b) a WHOLE-TABLE structural negative control: the classical Kronecker-Hurwitz class
      number relation (tier-L ansatz, cited, NOT built into the counting code above --
      the counting code above knows nothing about this relation)
          sum_{s in Z, s^2<=4n} H(4n-s^2) = 2*sigma_1(n) - sum_{d|n} min(d, n/d)
      is checked against the H(D) table PRODUCED BY THE COUNTING METHOD ABOVE, for
      n=1..10. This is a genuine test of the counting/weighting code (a weight error at
      a=c,b=0 or a=b=c would break it) that individual spot-checks would not catch,
      since it sums many H(D) values into one integer identity.

Run: python hurwitz_class_numbers.py
Writes hurwitz_class_numbers_results.json.
"""
import json
from fractions import Fraction as Fr

def H(D):
    if D == 0:
        return Fr(-1, 12), []
    if D % 4 not in (0, 3):
        return Fr(0), []
    total = Fr(0)
    forms = []
    amax = int((D / 3) ** 0.5) + 2
    for a in range(1, amax + 1):
        for b in range(-a + 1, a + 1):
            num = b * b + D
            if num % (4 * a) != 0:
                continue
            c = num // (4 * a)
            if c < a:
                continue
            if c == a and b < 0:
                continue
            if b == 0 and a == c:
                w = Fr(1, 2)
            elif a == b and a == c:
                w = Fr(1, 3)
            else:
                w = Fr(1)
            total += w
            forms.append((a, b, c, str(w)))
    return total, forms

DMAX = 40
result = {}
Htab = {0: Fr(-1, 12)}
for D in range(0, DMAX + 1):
    if D != 0 and D % 4 not in (0, 3):
        continue
    val, forms = H(D)
    Htab[D] = val
    result[str(D)] = {"H": str(val), "reduced_forms": forms}

# ---- whole-table negative control: Kronecker-Hurwitz class number relation
# sum_{s in Z, s^2<=4n} H(4n-s^2) = 2*sigma_1(n) - sum_{d|n} min(d,n/d)
def sigma1(n):
    return sum(d for d in range(1, n + 1) if n % d == 0)

def min_sum(n):
    return sum(min(d, n // d) for d in range(1, n + 1) if n % d == 0)

def Hlook(D):
    if D < 0:
        return Fr(0)
    return Htab.get(D, None)

khr_checks = {}
NMAX = 10
khr_ok = True
for n in range(1, NMAX + 1):
    smax = int((4 * n) ** 0.5) + 1
    total = Fr(0)
    missing = []
    for s in range(-smax, smax + 1):
        D = 4 * n - s * s
        if D < 0:
            continue
        h = Hlook(D)
        if h is None:
            missing.append(D)
            continue
        total += h
    rhs = 2 * sigma1(n) - min_sum(n)
    ok = (not missing) and (total == rhs)
    khr_ok = khr_ok and ok
    khr_checks[n] = {
        "lhs_sum_H(4n-s^2)": str(total), "rhs_2sigma1(n)-sum_min(d,n/d)": str(rhs),
        "sigma1(n)": sigma1(n), "sum_min(d,n_over_d)": min_sum(n),
        "missing_H_values_needed": missing, "matches": ok,
    }

# negative control on the relation itself: perturb H(3) by +1/6 (an amount NOT a multiple
# of the weights that appear, so it cannot cancel) and confirm n=1's relation breaks
Htab_pert = dict(Htab)
Htab_pert[3] = Htab_pert[3] + Fr(1, 6)
def Hlook_pert(D):
    if D < 0:
        return Fr(0)
    return Htab_pert.get(D, None)
n = 1
smax = int((4 * n) ** 0.5) + 1
total_pert = sum((Hlook_pert(4 * n - s * s) or Fr(0)) for s in range(-smax, smax + 1) if 4 * n - s * s >= 0)
rhs1 = 2 * sigma1(1) - min_sum(1)
relation_breaks_under_perturbation = (total_pert != rhs1)

out = {
    "DMAX": DMAX,
    "H_table": result,
    "kronecker_hurwitz_relation": {
        "statement_tier_L_cited_not_assumed_in_counting_code": (
            "sum_{s in Z, s^2<=4n} H(4n-s^2) = 2*sigma_1(n) - sum_{d|n} min(d,n/d); "
            "classical (Kronecker 1860 / Hurwitz); used here only as an EXTERNAL "
            "structural check on the H(D) table produced by the independent reduced-form "
            "counting above, for n=1..10."
        ),
        "checks_n_1_to_10": {str(n): khr_checks[n] for n in range(1, NMAX + 1)},
        "all_match": khr_ok,
        "negative_control_perturb_H3_by_1_6": {
            "n=1_lhs_after_perturbation": str(total_pert), "n=1_rhs": str(rhs1),
            "relation_breaks": relation_breaks_under_perturbation,
        },
    },
}
with open("hurwitz_class_numbers_results.json", "w") as f:
    json.dump(out, f, indent=1)
print(json.dumps({"DMAX": DMAX, "khr_all_match": khr_ok,
                   "khr_breaks_under_perturbation": relation_breaks_under_perturbation}, indent=1))
