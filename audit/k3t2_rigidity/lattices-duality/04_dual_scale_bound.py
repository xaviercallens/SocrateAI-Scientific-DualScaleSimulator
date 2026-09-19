"""
Track C, item (4): dual-scale bound  f(G) = tr(G) + tr(G^{-1}) - 2d  for
random SPD matrices G, d = 1..6, 2000 seeded samples each.

Exact-arithmetic part: for each sampled SPD G (built from an integer matrix
M as G = M^T M + c I, so entries are exact integers/rationals), compute
tr(G) + tr(G^{-1}) - 2d EXACTLY via sympy Rational, and confirm it is >= 0
for every sample (this is an exact statement per sample, since it reduces to
comparing two rationals).

Rigidity / "equality iff G = I" part (necessarily float, since it's a
perturbation *scan*, not a single exact evaluation): scan G = I + eps * P for
a symmetric perturbation P and eps in a grid including 0, and show f(G) > 0
for every eps != 0 in the scan and f(I) = 0 exactly -- i.e. I is an isolated
minimum in the scanned directions.

Proof sketch (reported, not just numerically checked): G is SPD so it has an
orthonormal eigenbasis with eigenvalues x_1..x_d > 0, and tr G + tr G^{-1} =
sum_i (x_i + 1/x_i). For x > 0, x + 1/x - 2 = (x-1)^2 / x >= 0, with equality
iff x = 1. So f(G) = sum_i (x_i - 1)^2 / x_i >= 0, equality iff every
eigenvalue is 1, iff G = I (since G is diagonalizable with eigenvalues all
1 implies G = I for a symmetric matrix). This script verifies the algebraic
identity (x-1)^2/x = x + 1/x - 2 exactly (sympy) and the SPD-sample and
perturbation-scan claims numerically/exactly as described above.
"""
import json
import random
import sys

import sympy as sp
import numpy as np

sys.path.insert(0, "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality")

N_SAMPLES = 2000


def random_spd_exact(d, rng, entry_range=3, ridge=None):
    if ridge is None:
        ridge = d + 1
    M = sp.Matrix(d, d, lambda i, j: rng.randint(-entry_range, entry_range))
    G = M.T * M + ridge * sp.eye(d)
    return G


def exact_bound_check(d, rng, n_samples):
    min_f = None
    all_nonneg = True
    n_checked = 0
    for _ in range(n_samples):
        G = random_spd_exact(d, rng)
        Ginv = G.inv()
        f = sp.trace(G) + sp.trace(Ginv) - 2 * d
        f = sp.nsimplify(f)
        if f < 0:
            all_nonneg = False
        if (min_f is None) or (f < min_f):
            min_f = f
        n_checked += 1
    return {
        "n_samples": n_checked,
        "all_f_nonnegative": all_nonneg,
        "min_f_observed_exact": str(min_f),
    }


def identity_algebraic_check():
    x = sp.symbols('x', positive=True)
    lhs = x + 1 / x - 2
    rhs = (x - 1) ** 2 / x
    identity_holds = sp.simplify(lhs - rhs) == 0
    return bool(identity_holds)


def perturbation_scan(d, rng, n_directions=20, eps_grid=None):
    """Scan G = I + eps*P for random symmetric P, eps over a grid including 0
    and small nonzero values (both signs). Verify f(I)=0 exactly and f(G) > 0
    for every eps != 0 small enough that G stays SPD (float, as noted)."""
    if eps_grid is None:
        eps_grid = [-0.2, -0.1, -0.05, -0.01, 0.0, 0.01, 0.05, 0.1, 0.2]
    I_d = np.eye(d)
    f_at_identity = float(np.trace(I_d) + np.trace(np.linalg.inv(I_d)) - 2 * d)
    violations = []
    min_nonzero_f = None
    checked = 0
    for _ in range(n_directions):
        P = rng.normal(size=(d, d))
        P = (P + P.T) / 2.0
        for eps in eps_grid:
            G = I_d + eps * P
            eigvals = np.linalg.eigvalsh(G)
            if np.any(eigvals <= 1e-9):
                continue  # not SPD at this eps; skip (out of domain)
            f = float(np.trace(G) + np.trace(np.linalg.inv(G)) - 2 * d)
            checked += 1
            if eps == 0.0:
                if abs(f) > 1e-8:
                    violations.append({"eps": eps, "f": f, "reason": "f(I) should be exactly 0"})
            else:
                if min_nonzero_f is None or f < min_nonzero_f:
                    min_nonzero_f = f
                if f <= 1e-10:
                    violations.append({"eps": eps, "f": f, "reason": "f(G)<=0 for eps!=0"})
    return {
        "f_at_identity_exact_formula": f_at_identity,
        "n_perturbations_checked": checked,
        "min_f_for_eps_neq_0": min_nonzero_f,
        "violations": violations,
        "equality_iff_G_eq_I_supported": (len(violations) == 0) and (min_nonzero_f is not None) and (min_nonzero_f > 0),
    }


results = {"track": "C", "item": "4_dual_scale_bound",
           "method": "Exact sympy Rational for random-SPD-sample bound (tr G + tr G^-1 - 2d >= 0, "
                     "each sample checked exactly); algebraic identity x+1/x-2=(x-1)^2/x checked "
                     "exactly via sympy; equality-iff-G=I demonstrated by a float perturbation "
                     "scan around G=I (explicitly marked non-exact -- eigenvalue signs / trace "
                     "bound are the only floats in this run, per ground rules).",
           "algebraic_identity_x_plus_1_over_x_minus_2_eq_square_over_x": identity_algebraic_check(),
           "by_dimension": {}}

master_rng_seed = 20260918
for d in range(1, 7):
    rng_exact = random.Random(master_rng_seed + 100 * d)
    exact_res = exact_bound_check(d, rng_exact, N_SAMPLES)
    exact_res["exact"] = True

    rng_pert = np.random.default_rng(master_rng_seed + 100 * d + 1)
    pert_res = perturbation_scan(d, rng_pert)
    pert_res["exact"] = False  # float scan, per ground rules

    results["by_dimension"][str(d)] = {
        "random_spd_sample_bound": exact_res,
        "perturbation_scan_equality_check": pert_res,
    }

out_path = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/lattices-duality/04_dual_scale_bound_result.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)

print(json.dumps(results, indent=2))
