#!/usr/bin/env python3
"""
Round 1 FORWARD LOOP step 2 (RESPONSE RANK).

Numerical Jacobian of a standardized observable vector w.r.t. ln(param),
evaluated at the sweep's best-fit (chi2-minimizing) numerically-stable point
and at 5 random numerically-stable sweep points (seed 42). Central
differences with step h=1e-2 in ln-space (NOT 1e-6: the ODE/BVP solves have
finite tolerance, so a tiny step returns solver noise, not the true
derivative -- see advisor note in this round's transcript).

Observable vector (11 components), each log-transformed first if it spans
multiple decades over the sweep, then standardized (z-score) using the
sweep's own mean/std for that observable:
  w0_cpl_latetime, wa_cpl_latetime,
  H_of_z_over_H0 at z=1.0 and z=2.3 (from full grid),
  D_M_times_H0 at z=1.0 and z=2.3,
  log10(screening_suppression_factor), log10(phi_center_ratio + 1e-300),
  max_deviation_from_hd, c4_pta_product,
  legacy_w0_fit (regression signal only, kept for completeness)

Writes: audit/zero_param_loop/round1/jacobian_report.json
"""
import json
import math
import multiprocessing as mp
import os
import sys

import numpy as np
import pandas as pd

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)
sys.path.insert(0, os.path.join(REPO_ROOT, "scripts"))
OUT_DIR = os.path.dirname(os.path.abspath(__file__))

import param_loop_sim as pls  # noqa: E402

PARAM_NAMES = ["a_pot", "b_pot", "mu_sym", "lambda_sym", "pta_suppression", "c4_c0_ratio"]
H_STEP_LN = 1e-2  # central-difference step in ln(param) space


def observable_vector(params):
    """Evaluate the harness at `params` and return the 11-component raw
    observable vector (before log-transform/standardization), plus a
    numerically_stable flag."""
    r = pls.evaluate_point(params)
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
        pta["max_deviation_from_hd"], pta["c4_pta_product"],
        de["legacy_w0_fit"],
    ], dtype=float)
    stable = bool(sc["numerically_stable"]) and np.all(np.isfinite(vec))
    return vec, stable


OBS_NAMES = ["w0_cpl", "wa_cpl", "H_z1", "H_z2.3", "DM_z1", "DM_z2.3",
             "log10_ssf", "log10_pcr", "pta_max_dev", "c4_pta_prod", "legacy_w0"]


def jacobian_at(params, sweep_std, sweep_mean):
    """Central-difference Jacobian of the STANDARDIZED observable vector
    w.r.t. ln(param_i), for i in PARAM_NAMES. Returns (J [11x6], base_vec_raw, ok)."""
    base_vec, base_stable = observable_vector(params)
    if not base_stable:
        return None, base_vec, False
    n_obs = len(base_vec)
    J = np.zeros((n_obs, len(PARAM_NAMES)))
    for j, pname in enumerate(PARAM_NAMES):
        p_plus = dict(params)
        p_minus = dict(params)
        p_plus[pname] = params[pname] * math.exp(H_STEP_LN)
        p_minus[pname] = params[pname] * math.exp(-H_STEP_LN)
        v_plus, ok_plus = observable_vector(p_plus)
        v_minus, ok_minus = observable_vector(p_minus)
        if not (ok_plus and ok_minus):
            J[:, j] = np.nan
            continue
        # standardize both, then differentiate the standardized quantity
        v_plus_std = (v_plus - sweep_mean) / sweep_std
        v_minus_std = (v_minus - sweep_mean) / sweep_std
        J[:, j] = (v_plus_std - v_minus_std) / (2.0 * H_STEP_LN)
    return J, base_vec, True


def _pooled_observable_vector(params):
    """Top-level (picklable) wrapper for multiprocessing.Pool.map. Relies on
    the fork start method (default on Linux) inheriting the already-imported
    `pls` module (imported at the top of this file) from the parent
    process's memory into each forked worker."""
    try:
        return observable_vector(params)
    except Exception:
        return np.full(len(OBS_NAMES), np.nan), False


def run_jacobian():
    df = pd.read_csv(os.path.join(OUT_DIR, "sweep.csv"))
    clean = df[df["error"].isna()].copy()
    stable = clean[clean["numerically_stable"] == True].copy()  # noqa: E712
    chi_df = pd.read_csv(os.path.join(OUT_DIR, "sweep_chi2.csv"))
    chi_stable = chi_df[chi_df["numerically_stable"] == True]  # noqa: E712
    best_idx = int(chi_stable.loc[chi_stable["chi2_total"].idxmin(), "idx"])
    best_row = stable[stable["idx"] == best_idx].iloc[0]
    best_params = {p: float(best_row[p]) for p in PARAM_NAMES}

    rng = np.random.RandomState(42)
    candidate_idxs = stable["idx"].to_numpy()
    candidate_idxs = candidate_idxs[candidate_idxs != best_idx]
    random_idxs = rng.choice(candidate_idxs, size=min(5, len(candidate_idxs)), replace=False)

    # sweep_mean/std for standardization, computed over ALL numerically_stable
    # sweep points' raw observable vectors (recomputed here directly from
    # harness, not re-derived from CSV columns, to keep the 11-vector self-consistent).
    print("[jacobian] recomputing raw observable vectors for standardization stats "
          f"over {len(stable)} numerically-stable sweep points (parallel, 6 workers) ...")
    param_dicts = [{p: float(row[p]) for p in PARAM_NAMES} for _, row in stable.iterrows()]
    with mp.Pool(processes=6, maxtasksperchild=8) as pool:
        pooled = pool.map(_pooled_observable_vector, param_dicts, chunksize=1)
    raw_vecs = [v for v, ok in pooled if ok]
    raw_vecs = np.array(raw_vecs)
    sweep_mean = np.nanmean(raw_vecs, axis=0)
    sweep_std = np.nanstd(raw_vecs, axis=0)
    sweep_std[sweep_std < 1e-12] = 1.0  # guard against a constant observable

    points_to_probe = [("best_fit", best_params)] + [
        ("random_%d" % idx, {p: float(stable[stable["idx"] == idx].iloc[0][p]) for p in PARAM_NAMES})
        for idx in random_idxs
    ]

    all_results = {}
    for label, params in points_to_probe:
        J, base_vec, ok = jacobian_at(params, sweep_std, sweep_mean)
        if not ok or J is None or not np.all(np.isfinite(J)):
            all_results[label] = {"params": params, "ok": False, "reason": "unstable or nonfinite Jacobian"}
            continue
        U, S, Vt = np.linalg.svd(J, full_matrices=False)
        thresh = 1e-3 * S[0]
        eff_dim = int(np.sum(S > thresh))
        null_dirs = []
        for k in range(len(S)):
            if S[k] <= thresh:
                loadings = {PARAM_NAMES[i]: float(Vt[k, i]) for i in range(len(PARAM_NAMES))}
                null_dirs.append({"singular_value": float(S[k]), "loadings_on_ln_param": loadings})
        all_results[label] = {
            "params": params, "ok": True,
            "singular_values": S.tolist(),
            "effective_dimension_at_1e-3_of_max": eff_dim,
            "near_null_right_singular_vectors": null_dirs,
            "right_singular_vectors_all": {
                f"sv_{k}": {PARAM_NAMES[i]: float(Vt[k, i]) for i in range(len(PARAM_NAMES))}
                for k in range(len(S))
            },
        }

    report = {
        "param_names": PARAM_NAMES,
        "observable_names": OBS_NAMES,
        "step_ln_param": H_STEP_LN,
        "n_stable_sweep_points_used_for_standardization": int(len(raw_vecs)),
        "sweep_mean_raw_observables": sweep_mean.tolist(),
        "sweep_std_raw_observables": sweep_std.tolist(),
        "probed_points": all_results,
    }
    with open(os.path.join(OUT_DIR, "jacobian_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: (v if k != "probed_points" else {
        kk: {"effective_dimension_at_1e-3_of_max": vv.get("effective_dimension_at_1e-3_of_max"),
             "singular_values": vv.get("singular_values")}
        for kk, vv in v.items()}) for k, v in report.items() if k != "sweep_mean_raw_observables" and k != "sweep_std_raw_observables"}, indent=2))
    return report


if __name__ == "__main__":
    run_jacobian()
