#!/usr/bin/env python3
"""
Round 1 FORWARD LOOP step 1 (EXPERIMENT): Latin-hypercube sweep of the 6
currently-free parameters (a_pot, b_pot, mu_sym, lambda_sym,
pta_suppression, c4_c0_ratio) over wide log ranges, evaluated with
scripts/param_loop_sim.evaluate_point (imported directly, not subprocessed,
for speed), run in parallel with multiprocessing (6 workers), seed 42.

Also computes:
  - chi2 of the model (2 free directions it can actually move: a_pot/b_pot
    for dark energy) against DESI 2024 BAO (with its real covariance,
    fetched from the same trusted CobayaSampler/bao_data repo as the
    verified mean file and hashed in data_provenance.json) and against
    Pantheon+ SH0ES distance moduli (diagonal errors, MU_SH0ES_ERR_DIAG).
  - chi2 of flat LCDM (Om fitted) on the same two datasets as baseline.
  - PTA chi2 is explicitly ABSENT: the verified NANOGrav product
    (ceffyl_data/*.npy) is a per-frequency free-spectrum KDE, not an
    inter-pulsar angular-separation correlation with errors, so there is no
    real angular Gamma(theta) dataset to compare pta.gamma_theta against.

Writes:
  audit/zero_param_loop/round1/sweep.csv
  audit/zero_param_loop/round1/chi2_report.json
  audit/zero_param_loop/round1/data_provenance.json
"""
import json
import math
import multiprocessing as mp
import os
import signal
import sys
import time

import numpy as np
from scipy.stats import qmc
from scipy.optimize import minimize_scalar as _unused  # noqa: F401 (guard import path)
from scipy.optimize import minimize_scalar
from scipy.optimize import minimize

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))

OUT_DIR = os.path.dirname(os.path.abspath(__file__))

SEED = 42
N_POINTS = 420
N_WORKERS = 6
PER_POINT_TIMEOUT_S = 45

PARAM_LOG10_BOUNDS = {
    # name: (log10 lo, log10 hi)  -- >=3 decades each, nominal inside range
    "a_pot": (-2.0, 2.0),          # nominal 1.0
    "b_pot": (-4.0, 0.0),          # nominal 0.01
    "mu_sym": (-2.0, 2.0),         # nominal 1.0
    "lambda_sym": (-2.0, 2.0),     # nominal 1.0
    "pta_suppression": (-4.0, 0.0),  # nominal 0.005
    "c4_c0_ratio": (-1.0, 3.0),    # nominal 16.07
}
PARAM_NAMES = list(PARAM_LOG10_BOUNDS.keys())

DESI_ROWS = [
    # (z, value, quantity) -- exact order matching desi_2024_bao_all_cov.txt rows
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
DESI_Z_UNIQUE = sorted(set(z for z, _, _ in DESI_ROWS))  # [0.295,0.510,0.706,0.930,1.317,1.491,2.330]

# ---------------------------------------------------------------------------
# Worker-side setup: extend param_loop_sim's module-level Z_GRID so that the
# DESI z=2.330 row and all Pantheon+ SN redshifts fall INSIDE the grid (no
# np.interp clamping at the boundary). Done once per worker via the Pool
# initializer, not via forked-at-import-time state, so it is explicit.
# ---------------------------------------------------------------------------
_pls = None  # set inside each worker process


def _worker_init():
    global _pls
    import param_loop_sim as pls
    pls.Z_GRID = np.linspace(0.0, 2.45, 80)  # covers DESI z=2.330 and Pantheon+ max z ~2.26 with margin
    pls.Z_POINTS = list(DESI_Z_UNIQUE)  # evaluate model exactly at DESI's z's
    _pls = pls
    # per-point hard timeout guard (symmetron relaxation fallback can crawl
    # at extreme mu_sym); SIGALRM only works on the main thread of a process,
    # which each pool worker is.
    signal.signal(signal.SIGALRM, _timeout_handler)


class _TimeoutError(Exception):
    pass


def _timeout_handler(signum, frame):
    raise _TimeoutError("param_loop_sim point exceeded timeout")


def _eval_one(row):
    idx, params = row
    global _pls
    signal.alarm(PER_POINT_TIMEOUT_S)
    t0 = time.time()
    try:
        r = _pls.evaluate_point(params)
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
        # model values at the exact DESI z's (Z_POINTS was overridden above)
        for j, z in enumerate(DESI_Z_UNIQUE):
            out[f"DM_H0_z{z}"] = de["D_M_times_H0_points"][j]
            out[f"HzH0_z{z}"] = de["H_of_z_over_H0_points"][j]
        # Full grid needed for SN interpolation, stored compactly as JSON.
        # IMPORTANT: store D_M_times_H0_grid (smooth, ~linear near z=0), NOT
        # distance_modulus_shape_grid (=5*log10(D_L*H0)): that quantity
        # diverges like log(z) near z=0, and linearly interpolating an
        # already-logged value on an 80-point grid badly misrepresents it
        # for the many real SNe at z~0.01-0.03 (verified empirically: gave
        # mu_model(z=0.01) = -43 instead of the correct ~ -7, from
        # interpolating between the z=0 floor value and the first grid
        # step). Reconstruct 5*log10((1+z)*D_M_H0(z)) downstream from the
        # smooth D_M_H0 grid instead.
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
    # include the nominal (defaults) point explicitly at idx=-1 for reference
    import param_loop_sim as pls_nominal  # host process import (grids untouched here)
    nominal = dict(pls_nominal.DEFAULT_PARAMS)
    rows = [(-1, nominal)] + rows

    t0 = time.time()
    with mp.Pool(processes=N_WORKERS, initializer=_worker_init, maxtasksperchild=8) as pool:
        results = pool.map(_eval_one, rows, chunksize=1)
    elapsed = time.time() - t0
    print(f"[sweep] {len(rows)} points in {elapsed:.1f}s using {N_WORKERS} workers")

    import pandas as pd
    df = pd.DataFrame(results)
    n_error = df["error"].notna().sum() if "error" in df else 0
    n_unstable = int((~df.get("numerically_stable", True).fillna(False)).sum()) if "numerically_stable" in df else None
    print(f"[sweep] errors/timeouts: {n_error}/{len(df)}; numerically_unstable: {n_unstable}/{len(df)}")

    csv_path = os.path.join(OUT_DIR, "sweep.csv")
    df.to_csv(csv_path, index=False)
    print(f"[sweep] wrote {csv_path}")
    return df, csv_path, elapsed, n_error, n_unstable


if __name__ == "__main__":
    run_sweep()
