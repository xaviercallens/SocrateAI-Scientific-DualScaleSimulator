"""
Track C, rigidity parameter (b): T-duality radius map  R -> lambda * alpha' / R.

lambda is treated as a free parameter and scanned over a mechanically
generated set of rationals (not hand-picked around the "known" answer).
Two SEPARATE conditions are checked symbolically/exactly (sympy, R, n, w,
alpha' kept as formal symbols -- no floats, no sampled points):

  (I) INVOLUTION: f(f(R)) == R identically, where f(R) = lambda*alpha'/R.

  (II) SPECTRUM PRESERVATION: the closed-string mass formula
       m^2(R; n, w) = (n/R)^2 + (w R / alpha')^2
       is preserved under R -> R' = lambda*alpha'/R together with the
       momentum/winding exchange (n,w) -> (w,n), i.e.
       m^2(R'; w, n) == m^2(R; n, w) identically in n, w, R, alpha'
       (checked by exact symbolic simplification of the difference to 0,
       not by numeric sampling).

  (III) PHYSICAL POSITIVITY: R > 0 must map to R' = lambda*alpha'/R > 0 for
        all R>0, alpha'>0 -- i.e. lambda > 0 -- since a radius is a physical
        length.

The reported solution set is the intersection of (I), (II), (III).
Negative control: lambda values that fail (II) or (III) are listed
explicitly from the same scan (not filtered out silently).
"""
import json
import sys
from fractions import Fraction

import sympy as sp

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality")

R, n, w, ap = sp.symbols('R n w alpha_prime', positive=True)
# n, w are momentum/winding integers; treat as formal positive symbols for
# the algebraic identity check (the identity must hold for the symbols, which
# implies it holds for every integer instantiation -- stronger than sampling).
# lambda kept as a plain rational substituted per scan point.


def mass2(Rval, nval, wval):
    return (nval / Rval) ** 2 + (wval * Rval / ap) ** 2


def check_lambda(lam):
    lam = sp.nsimplify(lam)
    f = lam * ap / R

    # (I) involution: f(f(R)) == R
    f_of_f = f.subs(R, f)
    involution_diff = sp.simplify(f_of_f - R)
    involution_holds = (involution_diff == 0)

    # (II) spectrum preservation under R->f(R), (n,w)->(w,n)
    Rprime = f
    m2_R = mass2(R, n, w)
    m2_Rprime_swapped = mass2(Rprime, w, n)
    spectrum_diff = sp.simplify(sp.expand(m2_Rprime_swapped - m2_R))
    spectrum_holds = (spectrum_diff == 0)

    # (III) positivity: lambda > 0 (algebraic condition on the parameter itself;
    # exact rational comparison, no floats)
    positivity_holds = bool(lam > 0)

    all_hold = involution_holds and spectrum_holds and positivity_holds

    return {
        "lambda": str(lam),
        "involution_holds": bool(involution_holds),
        "involution_residual": str(involution_diff),
        "spectrum_preservation_holds": bool(spectrum_holds),
        "spectrum_residual": str(spectrum_diff),
        "positivity_holds": positivity_holds,
        "all_conditions_hold": bool(all_hold),
    }


# mechanically generated scan set: all lambda = p/q for p in -4..4, q in 1..2,
# p!=0 (lambda=0 degenerate/undefined map), deduplicated -- includes the
# "true" value 1 and its neighbors/negatives without being hand-curated to
# only include plausible candidates.
scan_vals = set()
for p in range(-4, 5):
    if p == 0:
        continue
    for q in (1, 2):
        scan_vals.add(Fraction(p, q))
scan_vals = sorted(scan_vals)

scan_results = [check_lambda(sp.Rational(v.numerator, v.denominator)) for v in scan_vals]

solution_set = [r["lambda"] for r in scan_results if r["all_conditions_hold"]]
failing = [r for r in scan_results if not r["all_conditions_hold"]]

# Solution set under the task's LITERAL stated conditions -- "an involution
# preserving the spectrum" -- which is conditions (I) and (II) ONLY.
# Positivity (III) is an ADDED physical condition (a radius must be a
# positive length), not one of the two conditions the task names.
solution_set_I_and_II_only = [r["lambda"] for r in scan_results
                               if r["involution_holds"] and r["spectrum_preservation_holds"]]
task_literal_conditions_fix_lambda_uniquely = (len(solution_set_I_and_II_only) == 1)

results = {
    "track": "C",
    "item": "rigidity_b_tduality_radius_map",
    "rigidity_test": {
        "parameter_inserted": "lambda in R -> lambda*alpha'/R",
        "scanned": f"{len(scan_vals)} rationals p/q, p in -4..4 (p!=0), q in {{1,2}}: {[str(v) for v in scan_vals]}",
        "conditions_checked": [
            "(I) involution: f(f(R)) == R identically (symbolic)",
            "(II) spectrum preservation: m^2 formula invariant under R->f(R), (n,w)->(w,n), identically (symbolic)",
            "(III) positivity: lambda > 0 (R>0 must map to R'>0)",
        ],
        "solution_set_under_task_literal_conditions_I_and_II_only": solution_set_I_and_II_only,
        "task_literal_conditions_fix_lambda_uniquely": bool(task_literal_conditions_fix_lambda_uniquely),
        "solution_set_with_added_physical_positivity_condition_III": solution_set,
        "negative_control_failing_candidates_under_full_condition_set": [
            {"lambda": r["lambda"], "involution": r["involution_holds"],
             "spectrum": r["spectrum_preservation_holds"], "positivity": r["positivity_holds"]}
            for r in failing
        ],
        "zero_free_parameters_supported_by_task_literal_conditions_alone": bool(
            task_literal_conditions_fix_lambda_uniquely
        ),
        "zero_free_parameters_supported_with_added_positivity_condition": bool(
            len(solution_set) == 1 and len(failing) >= 1
        ),
        "honest_summary": (
            f"The task's two literal conditions (involution + spectrum preservation) leave "
            f"{len(solution_set_I_and_II_only)} solutions: {solution_set_I_and_II_only} "
            "(lambda^2=1, a discrete Z2 sign ambiguity -- NOT a single point). Condition (III), "
            "radius positivity, is an ADDED physical input beyond what the task names; only with "
            "it does the solution set collapse to the single point lambda=+1."
        ),
    },
    "full_scan": scan_results,
    "note": "Condition (I) holds algebraically for EVERY nonzero lambda (f(f(R))=lambda*alpha'/(lambda*alpha'/R)=R "
            "identically) -- it does not by itself constrain lambda. Condition (II), the closed-string spectrum "
            "identity, holds exactly for lambda=+1 and lambda=-1 (both make (n/R)^2+(wR/ap)^2 symmetric), i.e. the "
            "physics only fixes lambda^2=1. Condition (III), radius positivity, then breaks the +-1 degeneracy and "
            "leaves the single physical solution lambda=+1 (the standard R -> alpha'/R self-T-duality map).",
    "expected": {
        "lambda": 1,
        "source": "standard bosonic-string T-duality self-dual map R -> alpha'/R (memory, unverified prior to computation)",
    },
}

out_path = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality/07_rigidity_b_result.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(json.dumps(results, indent=2))
