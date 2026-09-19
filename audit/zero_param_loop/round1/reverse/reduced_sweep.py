#!/usr/bin/env python3
"""
REVERSE LOOP: Latin-hypercube sweep of the 5 REDUCED free parameters
(a_pot, b_pot, mu_sym, lambda_sym, c4_pta_product), seed 43 (task explicitly
asks for "a fresh sweep of the reduced model (seed 43)" -- this OVERRIDES
the global ground-rule default of seed 42 for this specific fresh sweep).

Adapted directly from ../param_sweep.py: same DESI-z grid extension, same
per-point columns, same timeout/error handling. Only the sampled parameter
set changed (5-dim instead of 6-dim) and the model call goes through
reduced_model.evaluate_reduced_point instead of pls.evaluate_point.

c4_pta_product bounds: round 1 sampled pta_suppression in log10[-4,0] and
c4_c0_ratio in log10[-1,3] independently; their SUM in log10-space (i.e.
their PRODUCT) has support log10[-5,3] but a convolved (triangular-in-log)
density, not uniform. Here we sample c4_pta_product log-UNIFORMLY over
log10[-5,3] directly -- a genuinely different marginal density, flagged in
tda_prediction.json as an expected source of nonzero (but topology-neutral)
bottleneck distance vs round 1's model cloud.

Writes: reduced_sweep.csv
"""
import json
import multiprocessing as mp
import os
import signal
import sys
import time

import numpy as np
from scipy.stats import qmc

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

SEED = 43  # task override for this fresh reduced-model sweep
N_POINTS = 520  # margin above round-1's 316-stable-out-of-420 stability rate
N_WORKERS = 6
PER_POINT_TIMEOUT_S = 45

PARAM_LOG10_BOUNDS = {
    "a_pot": (-2.0, 2.0),
    "b_pot": (-4.0, 0.0),
    "mu_sym": (-2.0, 2.0),
    "lambda_sym": (-2.0, 2.0),
    "c4_pta_product": (-5.0, 3.0),
}
PARAM_NAMES = list(PARAM_LOG10_BOUNDS.keys())

DESI_ROWS = [
    (0.295, 7.92512927, "DV_over_rs"),
    (0.510, 13.62003080, "DM_over_rs"),
    (0.510, 20.98334647, "DH_over_rs"),
    (0.706, 16.84645313, "DM_over_rs"),
    (0.706, 20.07872919, "DH_over_rs"),
    (0.930, 21.70841761, "DM_over_rs"),
    (0.930, 17.87612922, "DH_over_rs"),
    (1.317, 27.78720817, "DM_over_rs"),
    (1.317, 13.82372285, "DH_over_rs"),
    (1.491, 26.07217182, "DV_over_rs"),
    (2.330, 39.70838281, "DM_over_rs"),
    (2.330, 8.52256583, "DH_over_rs"),
]
DESI_Z_UNIQUE = sorted(set(z for z, _, _ in DESI_ROWS))

_rm = None


def _worker_init():
    global _rm
    import param_loop_sim as pls
    pls.Z_GRID = np.linspace(0.0, 2.45, 80)
    pls.Z_POINTS = list(DESI_Z_UNIQUE)
    import reduced_model as rm
    _rm = rm
    signal.signal(signal.SIGALRM, _timeout_handler)


class _TimeoutError(Exception):
    pass


def _timeout_handler(signum, frame):
    raise _TimeoutError("reduced point exceeded timeout")


def _eval_one(row):
    idx, params = row
    global _rm
    signal.alarm(PER_POINT_TIMEOUT_S)
    t0 = time.time()
    try:
        r = _rm.evaluate_reduced_point(params)
        elapsed = time.time() - t0
        out = {"idx": idx, "error": None, "elapsed": elapsed}
        out.update(params)
        de = r["dark_energy"]
        sc = r["screening"]
        pta = r["pta"]
        out["success"] = bool(de["success"] and sc["success"])
        out["numerically_stable"] = bool(sc["numerically_stable"])
        out["a0_today_convention"] = de["a0_today_convention"]
        out["w0_cpl_latetime"] = de["w0_cpl_latetime"]
        out["wa_cpl_latetime"] = de["wa_cpl_latetime"]
        out["legacy_w0_fit"] = de["legacy_w0_fit"]
        out["legacy_wa_fit"] = de["legacy_wa_fit"]
        for z in DESI_Z_UNIQUE:
            j = DESI_Z_UNIQUE.index(z)
            out[f"DM_H0_z{z}"] = de["D_M_times_H0_points"][j]
            out[f"HzH0_z{z}"] = de["H_of_z_over_H0_points"][j]
        out["z_grid_json"] = json.dumps(de["z_grid"])
        out["DM_H0_grid_json"] = json.dumps(de["D_M_times_H0_grid"])
        out["screening_suppression_factor"] = sc["screening_suppression_factor"]
        out["phi_center_ratio"] = sc["phi_center_ratio"]
        out["is_screened"] = bool(sc["is_screened"])
        out["max_deviation_from_hd"] = pta["max_deviation_from_hd"]
        out["c4_pta_product"] = pta["c4_pta_product"]
        return out
    except _TimeoutError:
        out = {"idx": idx, "error": "timeout", "elapsed": time.time() - t0}
        out.update(params)
        return out
    except Exception as exc:  # pragma: no cover
        out = {"idx": idx, "error": f"{type(exc).__name__}: {exc}", "elapsed": time.time() - t0}
        out.update(params)
        return out
    finally:
        signal.alarm(0)


def sample_lhs(n, seed):
    sampler = qmc.LatinHypercube(d=len(PARAM_NAMES), seed=seed)
    u = sampler.random(n)
    rows = []
    for i in range(n):
        params = {}
        for j, name in enumerate(PARAM_NAMES):
            lo, hi = PARAM_LOG10_BOUNDS[name]
            log10_val = lo + u[i, j] * (hi - lo)
            params[name] = float(10.0 ** log10_val)
        rows.append((i, params))
    return rows


def run_sweep():
    rows = sample_lhs(N_POINTS, SEED)
    import reduced_model as rm_nominal
    rows = [(-1, dict(rm_nominal.REDUCED_DEFAULT_PARAMS))] + rows

    t0 = time.time()
    with mp.Pool(processes=N_WORKERS, initializer=_worker_init, maxtasksperchild=8) as pool:
        results = pool.map(_eval_one, rows, chunksize=1)
    elapsed = time.time() - t0
    print(f"[reduced_sweep] {len(rows)} points in {elapsed:.1f}s using {N_WORKERS} workers")

    import pandas as pd
    df = pd.DataFrame(results)
    n_error = int(df["error"].notna().sum()) if "error" in df else 0
    n_stable = int(df.get("numerically_stable", False).fillna(False).sum()) if "numerically_stable" in df else 0
    print(f"[reduced_sweep] errors/timeouts: {n_error}/{len(df)}; numerically_stable: {n_stable}/{len(df)}")

    csv_path = os.path.join(OUT_DIR, "reduced_sweep.csv")
    df.to_csv(csv_path, index=False)
    print(f"[reduced_sweep] wrote {csv_path}")
    return df, csv_path, elapsed, n_error, n_stable


if __name__ == "__main__":
    run_sweep()
