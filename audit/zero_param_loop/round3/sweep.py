#!/usr/bin/env python3
"""
ROUND 3 forward loop, step 1: sweep the CURRENTLY FREE parameters
(mu_sym, c4_pta_product) over wide log ranges. Everything else is held
at its ledger (accepted round-1/round-2) value:
  - a_pot, b_pot: the quintessence sector is DELETED (round-2 accepted,
    "delete_unobservable"); dark energy is a zero-parameter frozen flat
    LCDM(Om=0.31115, H0=67.661 km/s/Mpc from omegaLambda/hubbleRadius_m,
    DarkEnergyScale.lean:84 / SelfDualCutoff.lean:100) via
    round2/reverse/reduced_model.py, reused here unmodified.
  - lambda_sym: fixed at 1.0 (round-2 accepted "delete_unobservable" via
    exact psi-form algebra; the code-level psi-rewrite in
    workshopcosmo.py has STILL not landed this round either -- see
    ledger note below -- so this stays behaviourally-equivalent-but-not
    literally-deleted, tier B not A).
  - pta_suppression / c4_c0_ratio: only the product c4_pta_product is
    swept; the split denominator is fixed at pta_suppression=0.005
    (round-1 Check C: any positive split gives byte-identical gamma_theta,
    re-asserted below with a fresh split-invariance check).

Sampler: scipy.stats.qmc.Sobol(d=2, scramble=True, seed=42), m=9 -> 512
points (power-of-2, keeps Sobol's balance property), plus the nominal
point at idx 0 (idx 1..512 are the Sobol points).

Ranges (wide, >=3 decades each, log-uniform):
  mu_sym          in [1e-2, 1e2]   (4 decades) -- WIDER than round2's
                  reverse-loop stability-restricted [0.2,3.0]; the known
                  BVP relaxation-fallback failure mode (silent NaN with
                  workshopcosmo success=True, first caught round 1) is
                  expected to bite in the upper end -- filtered via the
                  numerically_stable field, fraction reported honestly.
  c4_pta_product  in [1e-5, 1e3]   (8 decades, same range round1/round2
                  forward already used for this product)

Parallel: multiprocessing.Pool(6).

Writes: sweep.csv, sweep_meta.json
"""
import csv
import json
import math
import multiprocessing as mp
import os
import sys
import time

import numpy as np
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
ROUND2_REVERSE_DIR = os.path.abspath(os.path.join(HERE, "..", "round2", "reverse"))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
sys.path.insert(0, ROUND2_REVERSE_DIR)

SEED = 42
M_SOBOL = 9  # 2**9 = 512
MU_LOG10_BOUNDS = (-2.0, 2.0)   # mu_sym in [1e-2, 1e2], 4 decades
C4_LOG10_BOUNDS = (-5.0, 3.0)  # c4_pta_product in [1e-5, 1e3], 8 decades
PTA_SUPPRESSION_SPLIT_DENOM = 0.005  # arbitrary positive split, matches DEFAULT_PARAMS

FIELDNAMES = [
    "idx", "mu_sym", "c4_pta_product", "c4_c0_ratio_split", "pta_suppression_split",
    "screening_suppression_factor", "phi_center_ratio", "numerically_stable",
    "log10_ssf", "log10_pcr",
] + [f"gamma_theta_{i}" for i in range(15)] + [
    "pta_max_deviation_from_hd", "wall_time_seconds", "error",
]


def _eval_one(args):
    idx, mu_sym, c4_pta_product = args
    import reduced_model as rm  # re-imported per worker (fork-safe, cheap)
    t0 = time.time()
    row = {"idx": idx, "mu_sym": mu_sym, "c4_pta_product": c4_pta_product, "error": ""}
    try:
        res = rm.evaluate_reduced_point({"mu_sym": mu_sym, "c4_pta_product": c4_pta_product})
        sc, pta = res["screening"], res["pta"]
        ssf = sc["screening_suppression_factor"]
        pcr = sc["phi_center_ratio"]
        row["screening_suppression_factor"] = ssf
        row["phi_center_ratio"] = pcr
        row["numerically_stable"] = bool(sc["numerically_stable"])
        row["log10_ssf"] = math.log10(ssf) if ssf > 0 and np.isfinite(ssf) else float("nan")
        row["log10_pcr"] = math.log10(pcr) if pcr > 0 and np.isfinite(pcr) else float("nan")
        gtheta = pta["gamma_theta"]
        for i in range(15):
            row[f"gamma_theta_{i}"] = gtheta[i]
        row["pta_max_deviation_from_hd"] = pta["max_deviation_from_hd"]
        row["c4_c0_ratio_split"] = pta["c4_c0_ratio"]
        row["pta_suppression_split"] = pta["pta_suppression"]
    except Exception as exc:
        row["error"] = str(exc)
        row["screening_suppression_factor"] = float("nan")
        row["phi_center_ratio"] = float("nan")
        row["numerically_stable"] = False
        row["log10_ssf"] = float("nan")
        row["log10_pcr"] = float("nan")
        for i in range(15):
            row[f"gamma_theta_{i}"] = float("nan")
        row["pta_max_deviation_from_hd"] = float("nan")
        row["c4_c0_ratio_split"] = float("nan")
        row["pta_suppression_split"] = float("nan")
    row["wall_time_seconds"] = time.time() - t0
    return idx, row


def split_invariance_check():
    """Re-assert round-1's Check C inside round 3: two different positive
    splits of the SAME c4_pta_product must give byte-identical gamma_theta,
    and a DIFFERENT product must move it (positive control)."""
    import reduced_model as rm
    p_nominal = rm.pls.compute_pta_observable(0.005, 16.07)          # product 0.08035
    p_altsplit = rm.pls.compute_pta_observable(0.0005, 160.7)        # same product, different split
    p_diffprod = rm.pls.compute_pta_observable(0.005, 32.14)         # product 0.1607, DIFFERENT
    g0 = np.array(p_nominal["gamma_theta"])
    g1 = np.array(p_altsplit["gamma_theta"])
    g2 = np.array(p_diffprod["gamma_theta"])
    max_diff_same_product = float(np.max(np.abs(g0 - g1)))
    max_diff_diff_product = float(np.max(np.abs(g0 - g2)))
    return {
        "same_product_different_split_max_diff": max_diff_same_product,
        "different_product_max_diff": max_diff_diff_product,
        "split_invariance_holds": bool(max_diff_same_product < 1e-9),
        "positive_control_product_change_moves_output": bool(max_diff_diff_product > 1e-3),
    }


def run():
    points = [(0, 0.08035 ** 0 * 1.0, 0.08035)]  # idx 0: nominal (mu_sym=1.0, product=0.08035)
    points[0] = (0, 1.0, 0.08035)

    sampler = qmc.Sobol(d=2, scramble=True, seed=SEED)
    u = sampler.random_base2(m=M_SOBOL)  # 512 points, preserves Sobol balance
    mu = 10 ** (MU_LOG10_BOUNDS[0] + u[:, 0] * (MU_LOG10_BOUNDS[1] - MU_LOG10_BOUNDS[0]))
    c4 = 10 ** (C4_LOG10_BOUNDS[0] + u[:, 1] * (C4_LOG10_BOUNDS[1] - C4_LOG10_BOUNDS[0]))
    for i in range(len(u)):
        points.append((i + 1, float(mu[i]), float(c4[i])))

    t0 = time.time()
    with mp.Pool(6) as pool:
        results = pool.map(_eval_one, points, chunksize=8)
    elapsed = time.time() - t0

    results.sort(key=lambda r: r[0])  # sort by idx: worker return order is not guaranteed
    rows = [r for _, r in results]

    out_csv = os.path.join(HERE, "sweep.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    n_total = len(rows)
    n_stable = sum(1 for r in rows if r["numerically_stable"])
    stable_rows = [r for r in rows if r["numerically_stable"]]
    mu_stable = [r["mu_sym"] for r in stable_rows]
    mu_unstable = [r["mu_sym"] for r in rows if not r["numerically_stable"]]

    split_check = split_invariance_check()

    meta = {
        "generated": "2026-09-18",
        "sampler": f"scipy.stats.qmc.Sobol(d=2, scramble=True, seed={SEED}).random_base2(m={M_SOBOL})",
        "n_points_total": n_total,
        "n_points_sobol": len(u),
        "n_numerically_stable": n_stable,
        "n_numerically_unstable": n_total - n_stable,
        "fraction_stable": n_stable / n_total,
        "wall_time_seconds_pool6": elapsed,
        "ranges": {
            "mu_sym_log10_bounds": MU_LOG10_BOUNDS,
            "mu_sym_bounds": [10 ** MU_LOG10_BOUNDS[0], 10 ** MU_LOG10_BOUNDS[1]],
            "c4_pta_product_log10_bounds": C4_LOG10_BOUNDS,
            "c4_pta_product_bounds": [10 ** C4_LOG10_BOUNDS[0], 10 ** C4_LOG10_BOUNDS[1]],
        },
        "mu_sym_stable_range_actually_explored": [
            float(np.min(mu_stable)) if mu_stable else None,
            float(np.max(mu_stable)) if mu_stable else None,
        ],
        "mu_sym_unstable_range": [
            float(np.min(mu_unstable)) if mu_unstable else None,
            float(np.max(mu_unstable)) if mu_unstable else None,
        ],
        "instability_mechanism": (
            "known symmetron BVP relaxation-fallback silent divergence to "
            "NaN/Inf at large mu_sym (workshopcosmo.run_symmetron_screening_"
            "simulation reports success=True regardless; caught via the "
            "numerically_stable field, first identified round 1/round 2)."
        ),
        "split_invariance_check_round3": split_check,
        "params_held_fixed_at_ledger_values": {
            "a_pot_b_pot_sector": "DELETED, replaced by frozen flat LCDM (round2/reverse/reduced_model.py, zero free params)",
            "lambda_sym": 1.0,
            "pta_suppression_split_denominator": PTA_SUPPRESSION_SPLIT_DENOM,
        },
        "lambda_sym_condition_status": (
            "ROUND-2 CONDITIONAL NOT MET: round2/ledger.json's convention_only_"
            "not_counted entry required the psi-form rewrite to 'land next "
            "round or this reverts to convention_only'. That code change "
            "(deleting lambda_sym from workshopcosmo.py's ODE and from "
            "param_loop_sim.py's DEFAULT_PARAMS/signature) was NOT made this "
            "round either -- ground rules restrict this worktree to audit/ "
            "experiment files, and no separate code-change round was run. "
            "Per the ground rules this reverts lambda_sym's removal to "
            "'fixed by convention', 0.5-evidence, reported separately, NOT "
            "counted in the free-parameter total. See round3 ledger update."
        ),
    }
    with open(os.path.join(HERE, "sweep_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Wrote {out_csv}: {n_total} rows ({n_stable} stable, {n_total-n_stable} unstable) in {elapsed:.1f}s (Pool(6))")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    run()
