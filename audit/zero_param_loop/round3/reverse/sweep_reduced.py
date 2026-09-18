#!/usr/bin/env python3
"""
ROUND 3 REVERSE loop, step 1: build the REDUCED model and re-sweep it
independently (fresh seed 43, distinct from the forward round's seed 42).

APPLYING THE FORWARD ROUND'S PROPOSED REDUCTIONS
--------------------------------------------------
Per the task brief's own "FORWARD RESULT" and this repo's
../../ledger.json (ladder round-3 entry, "how": "none -- no count
reduction this round"), NEITHER of round 3's two proposed items was
actually accepted as a tier A/B/L reduction:

  - mu_sym: "no value is proposed for mu_sym... mu_sym itself has no
    length calibration anywhere in workshopcosmo.py; no bridge formula
    exists" (tier C, NOT derived). Not removable under ground rule (a).
  - c4_pta_product: mechanism "convention" -- drop the PTA sector from
    RETAINED OBSERVABLES (no real angular Hellings-Downs dataset exists).
    This is a scope change (stop testing it), not rule (b) (the parameter
    still affects gamma_theta, it is just untested) or rule (a)/(c). Per
    the ground rules this does NOT count as a reduction.
  - lambda_sym: already excluded from "the 2" (mu_sym, c4_pta_product);
    its convention-only status is unchanged this round.

So the REDUCED model for this reverse pass is, honestly, IDENTICAL to the
round-3 FULL model: reduced_params = [mu_sym, c4_pta_product]. This
script re-verifies that null result independently by resampling the same
2D parameter space with a different Sobol seed (43 vs the forward round's
42) and, additionally, applying the S3 "scope, not derivation" convention
explicitly: the PTA columns are still COMPUTED (so we can check the
convention doesn't silently change any number) but are marked
`retained_for_chi2=False` in the meta file, documenting the scope
decision without pretending it removed a parameter.

Reuses round2/reverse/reduced_model.py unmodified (same frozen-LCDM dark
energy block, same screening/PTA wiring) -- no code in workshopcosmo.py
or param_loop_sim.py is touched, per ground rules.

Writes: sweep_reduced.csv, sweep_reduced_meta.json
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
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
ROUND2_REVERSE_DIR = os.path.abspath(os.path.join(HERE, "..", "..", "round2", "reverse"))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
sys.path.insert(0, ROUND2_REVERSE_DIR)

SEED = 43
M_SOBOL = 9  # 2**9 = 512, same sample size as forward round for a fair comparison
MU_LOG10_BOUNDS = (-2.0, 2.0)
C4_LOG10_BOUNDS = (-5.0, 3.0)
PTA_SUPPRESSION_SPLIT_DENOM = 0.005

FIELDNAMES = [
    "idx", "mu_sym", "c4_pta_product", "c4_c0_ratio_split", "pta_suppression_split",
    "screening_suppression_factor", "phi_center_ratio", "numerically_stable",
    "log10_ssf", "log10_pcr",
] + [f"gamma_theta_{i}" for i in range(15)] + [
    "pta_max_deviation_from_hd", "wall_time_seconds", "error",
]

REDUCED_PARAMS_STILL_FREE = ["mu_sym", "c4_pta_product"]  # unchanged: no reduction accepted


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


def run():
    points = [(0, 1.0, 0.08035)]  # idx 0: nominal
    sampler = qmc.Sobol(d=2, scramble=True, seed=SEED)
    u = sampler.random_base2(m=M_SOBOL)
    mu = 10 ** (MU_LOG10_BOUNDS[0] + u[:, 0] * (MU_LOG10_BOUNDS[1] - MU_LOG10_BOUNDS[0]))
    c4 = 10 ** (C4_LOG10_BOUNDS[0] + u[:, 1] * (C4_LOG10_BOUNDS[1] - C4_LOG10_BOUNDS[0]))
    for i in range(len(u)):
        points.append((i + 1, float(mu[i]), float(c4[i])))

    t0 = time.time()
    with mp.Pool(6) as pool:
        results = pool.map(_eval_one, points, chunksize=8)
    elapsed = time.time() - t0

    results.sort(key=lambda r: r[0])
    rows = [r for _, r in results]

    out_csv = os.path.join(HERE, "sweep_reduced.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDNAMES)
        w.writeheader()
        for row in rows:
            w.writerow(row)

    n_total = len(rows)
    n_stable = sum(1 for r in rows if r["numerically_stable"])
    stable_rows = [r for r in rows if r["numerically_stable"]]
    mu_stable = [r["mu_sym"] for r in stable_rows]

    meta = {
        "generated": "2026-09-18",
        "purpose": "REVERSE loop round 3, step 1: independent re-sweep (seed 43) of the REDUCED model, which is numerically identical to the round-3 forward FULL model because neither proposed reduction (mu_sym length bridge, c4_pta_product PTA-scope) was accepted per ../../ledger.json ladder round-3 entry.",
        "reduced_params_still_free": REDUCED_PARAMS_STILL_FREE,
        "reduction_count_change": "2 -> 2 (no change; this is the honest headline of round 3's reverse loop too)",
        "sampler": f"scipy.stats.qmc.Sobol(d=2, scramble=True, seed={SEED}).random_base2(m={M_SOBOL})",
        "n_points_total": n_total,
        "n_numerically_stable": n_stable,
        "fraction_stable": n_stable / n_total,
        "wall_time_seconds_pool6": elapsed,
        "mu_sym_stable_range_actually_explored": [
            float(np.min(mu_stable)) if mu_stable else None,
            float(np.max(mu_stable)) if mu_stable else None,
        ],
        "s3_scope_convention_applied": (
            "PTA gamma_theta columns are still computed and written (so the "
            "convention decision is auditable) but pta_max_deviation_from_hd "
            "is NOT used in chi2_aic_bic.py's fit (retained_for_chi2=False) "
            "per the S3 'scope, not derivation' convention -- this does not "
            "change the free-parameter count."
        ),
        "cross_seed_check_vs_forward_round3": "see chi2_aic_bic.py for the chi2 cross-check; TDA cross-seed comparison in tda_reduced.py.",
    }
    with open(os.path.join(HERE, "sweep_reduced_meta.json"), "w") as f:
        json.dump(meta, f, indent=2)

    print(f"Wrote {out_csv}: {n_total} rows ({n_stable} stable) in {elapsed:.2f}s")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    run()
