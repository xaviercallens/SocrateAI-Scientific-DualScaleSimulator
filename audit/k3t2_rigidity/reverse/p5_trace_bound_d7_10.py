"""
Reverse pass, prediction P5.

lattices-duality/04_dual_scale_bound.py established (exact sympy Rational,
2000 SPD samples per dimension) that tr(G) + tr(G^-1) - 2d >= 0 for d = 1..6,
with the algebraic identity x + 1/x - 2 = (x-1)^2/x verified symbolically
(so the bound is d-independent by construction: it is a per-eigenvalue
statement summed over d eigenvalues). PREDICTION: the same bound holds for
d = 7, 8, 9, 10 (the dimensions relevant to heterotic-on-K3xT2 / Gamma^{6,22}
constructions with extra spectator directions), tested the same way -- exact
Rational arithmetic, no floats except the (already-flagged) perturbation scan.
"""
import json
import random
import sys

import sympy as sp
import numpy as np

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
    min_f_matrix = None
    for _ in range(n_samples):
        G = random_spd_exact(d, rng)
        Ginv = G.inv()
        f = sp.trace(G) + sp.trace(Ginv) - 2 * d
        f = sp.nsimplify(f)
        if f < 0:
            all_nonneg = False
        if (min_f is None) or (f < min_f):
            min_f = f
            min_f_matrix = G
        n_checked += 1
    return {
        "n_samples": n_checked,
        "all_f_nonnegative": all_nonneg,
        "min_f_observed_exact": str(min_f),
    }


def perturbation_scan(d, rng, n_directions=20, eps_grid=None):
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
                continue
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


x = sp.symbols('x', positive=True)
identity_holds = bool(sp.simplify((x + 1 / x - 2) - (x - 1) ** 2 / x) == 0)

results = {"track": "reverse-pass extension of lattices-duality/04_dual_scale_bound.py",
           "prediction": "tr(G)+tr(G^-1)-2d >= 0 holds for d=7..10 (Lean-adjacent bound tr G + tr G^-1 >= 2d, "
                          "d<=6 previously checked)",
           "algebraic_identity_x_plus_1_over_x_minus_2_eq_square_over_x": identity_holds,
           "by_dimension": {}}

master_rng_seed = 20260918
for d in range(7, 11):
    rng_exact = random.Random(master_rng_seed + 100 * d)
    exact_res = exact_bound_check(d, rng_exact, N_SAMPLES)
    exact_res["exact"] = True

    rng_pert = np.random.default_rng(master_rng_seed + 100 * d + 1)
    pert_res = perturbation_scan(d, rng_pert)
    pert_res["exact"] = False

    results["by_dimension"][str(d)] = {
        "random_spd_sample_bound": exact_res,
        "perturbation_scan_equality_check": pert_res,
    }

out_path = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-k3t2/audit/k3t2_rigidity/reverse/p5_trace_bound_results.json"
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)
print(json.dumps(results, indent=2)[:3000])
