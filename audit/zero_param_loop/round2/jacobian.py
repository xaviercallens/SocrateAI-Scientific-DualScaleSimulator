#!/usr/bin/env python3
"""
ROUND 2 forward-loop step 2 (RESPONSE RANK): numerical Jacobian of the
standardized observable vector w.r.t. ln(params), at this round's sweep
best-fit point (sweep_chi2_report.json) and 5 random numerically-stable
points from sweep_chi2.csv, seed 42. Central differences, step 1e-2 in
ln-space, over the 5 currently free parameters (a_pot, b_pot, mu_sym,
lambda_sym, c4_pta_product). SVD -> singular values -> effective_dimension
= count of singular values > 1e-3 * max singular value, at each point.

Adapted from round1/reverse/reduced_jacobian.py (same method, same
observable set, fresh points/seed for this round).

Writes: jacobian_report.json
"""
import json
import math
import multiprocessing as mp
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
ROUND1_REVERSE_DIR = os.path.abspath(os.path.join(HERE, "..", "round1", "reverse"))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
sys.path.insert(0, ROUND1_REVERSE_DIR)

import param_loop_sim as pls  # noqa: E402
pls.Z_GRID = np.linspace(0.0, 2.45, 80)

import reduced_model as rm  # noqa: E402

PARAM_NAMES = ["a_pot", "b_pot", "mu_sym", "lambda_sym", "c4_pta_product"]
H_STEP_LN = 1e-2

OBS_NAMES = ["w0_cpl", "wa_cpl", "H_z1", "H_z2.3", "DM_z1", "DM_z2.3",
             "log10_ssf", "log10_pcr", "pta_max_dev", "legacy_w0"]


def observable_vector(params):
    r = rm.evaluate_reduced_point(params)
    de, sc, pta = r["dark_energy"], r["screening"], r["pta"]
    z_grid = np.array(de["z_grid"])
    Hz = np.array(de["H_of_z_over_H0_grid"])
    DM = np.array(de["D_M_times_H0_grid"])
    H_z1 = float(np.interp(1.0, z_grid, Hz))
    H_z23 = float(np.interp(2.3, z_grid, Hz))
    DM_z1 = float(np.interp(1.0, z_grid, DM))
    DM_z23 = float(np.interp(2.3, z_grid, DM))
    ssf = sc["screening_suppression_factor"]
    pcr = sc["phi_center_ratio"]
    vec = np.array([
        de["w0_cpl_latetime"], de["wa_cpl_latetime"],
        H_z1, H_z23, DM_z1, DM_z23,
        math.log10(ssf) if ssf > 0 and np.isfinite(ssf) else float("nan"),
        math.log10(pcr) if pcr > 0 and np.isfinite(pcr) else float("nan"),
        pta["max_deviation_from_hd"],
        de["legacy_w0_fit"],
    ], dtype=float)
    stable = bool(sc["numerically_stable"]) and np.all(np.isfinite(vec))
    return vec, stable


def _pooled(params):
    try:
        return observable_vector(params)
    except Exception:
        return np.full(len(OBS_NAMES), np.nan), False


def jacobian_at(params, sweep_std, sweep_mean):
    base_vec, base_stable = observable_vector(params)
    if not base_stable:
        return None, base_vec, False
    n_obs = len(base_vec)
    J = np.zeros((n_obs, len(PARAM_NAMES)))
    for j, pname in enumerate(PARAM_NAMES):
        p_plus = dict(params); p_minus = dict(params)
        p_plus[pname] = params[pname] * math.exp(H_STEP_LN)
        p_minus[pname] = params[pname] * math.exp(-H_STEP_LN)
        v_plus, ok_plus = observable_vector(p_plus)
        v_minus, ok_minus = observable_vector(p_minus)
        if not (ok_plus and ok_minus):
            J[:, j] = np.nan
            continue
        v_plus_std = (v_plus - sweep_mean) / sweep_std
        v_minus_std = (v_minus - sweep_mean) / sweep_std
        J[:, j] = (v_plus_std - v_minus_std) / (2.0 * H_STEP_LN)
    return J, base_vec, True


def run():
    df = pd.read_csv(os.path.join(HERE, "sweep_chi2.csv"))
    stable = df[df["numerically_stable"] == True].copy()  # noqa: E712
    print(f"[round2 jacobian] stable sweep rows: {len(stable)}")

    with open(os.path.join(HERE, "sweep_chi2_report.json")) as f:
        chi2_report = json.load(f)
    bf = chi2_report["best_fit_row"]
    best_params = {p: float(bf[p]) for p in PARAM_NAMES}

    rng = np.random.RandomState(42)
    idxs = stable["idx"].to_numpy()
    random_idxs = rng.choice(idxs, size=min(5, len(idxs)), replace=False)

    param_dicts = [{p: float(row[p]) for p in PARAM_NAMES} for _, row in stable.iterrows()]
    print(f"[round2 jacobian] computing standardization stats over {len(param_dicts)} points (6 workers) ...")
    with mp.Pool(processes=6, maxtasksperchild=8) as pool:
        pooled = pool.map(_pooled, param_dicts, chunksize=1)
    raw_vecs = np.array([v for v, ok in pooled if ok])
    sweep_mean = np.nanmean(raw_vecs, axis=0)
    sweep_std = np.nanstd(raw_vecs, axis=0)
    zero_var_cols = [OBS_NAMES[i] for i in range(len(OBS_NAMES)) if sweep_std[i] < 1e-12]
    sweep_std[sweep_std < 1e-12] = 1.0  # guard: would otherwise give 0/0 -> NaN SVD garbage

    points_to_probe = [("best_fit", best_params)] + [
        ("random_%d" % idx, {p: float(stable[stable["idx"] == idx].iloc[0][p]) for p in PARAM_NAMES})
        for idx in random_idxs
    ]

    all_results = {}
    eff_dims = {}
    for label, params in points_to_probe:
        J, base_vec, ok = jacobian_at(params, sweep_std, sweep_mean)
        if not ok or J is None or not np.all(np.isfinite(J)):
            all_results[label] = {"params": params, "ok": False, "reason": "unstable or nonfinite Jacobian"}
            continue
        U, S, Vt = np.linalg.svd(J, full_matrices=False)
        thresh = 1e-3 * S[0]
        eff_dim = int(np.sum(S > thresh))
        eff_dims[label] = eff_dim
        null_dirs = [
            {"singular_value": float(S[k]),
             "loadings_on_ln_param": {PARAM_NAMES[i]: float(Vt[k, i]) for i in range(len(PARAM_NAMES))}}
            for k in range(len(S)) if S[k] <= thresh
        ]
        all_results[label] = {
            "params": params, "ok": True,
            "singular_values": S.tolist(),
            "effective_dimension_at_1e-3_of_max": eff_dim,
            "near_null_right_singular_vectors": null_dirs,
        }

    report = {
        "param_names": PARAM_NAMES,
        "observable_names": OBS_NAMES,
        "zero_variance_observable_columns_guarded": zero_var_cols,
        "step_ln_param": H_STEP_LN,
        "n_stable_sweep_points_used_for_standardization": int(len(raw_vecs)),
        "probed_points": all_results,
        "effective_dimension_summary": eff_dims,
        "effective_dimension_best_fit": eff_dims.get("best_fit"),
        "effective_dimension_spread_random": [v for k, v in eff_dims.items() if k != "best_fit"],
    }
    with open(os.path.join(HERE, "jacobian_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: (v if k != "probed_points" else {
        kk: {"effective_dimension_at_1e-3_of_max": vv.get("effective_dimension_at_1e-3_of_max"),
             "singular_values": vv.get("singular_values")}
        for kk, vv in v.items()}) for k, v in report.items()}, indent=2))
    return report


if __name__ == "__main__":
    run()
