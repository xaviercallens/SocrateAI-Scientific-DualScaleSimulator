#!/usr/bin/env python3
"""
ROUND 3 forward loop, step 2 (RESPONSE RANK): numerical Jacobian of the
standardized OBSERVABLE vector w.r.t. ln(mu_sym), ln(c4_pta_product), at
this round's sweep "best-fit" point and 5 random numerically-stable
points from sweep.csv, seed 42. Central differences, step 1e-2 in
ln-space.

"Best-fit" here means the nominal point (mu_sym=1.0, c4_pta_product=
0.08035): chi2.py showed chi2_total is EXACTLY CONSTANT across the whole
sweep (dark energy frozen, no screening/PTA dataset), so there is no
actual chi2-minimizing point to report -- using the nominal point is the
honest choice, stated explicitly rather than manufacturing a spurious
"best fit" from a flat likelihood surface.

Observable vector (11 components, dark-energy block DROPPED because it
is a zero-variance constant -- see round2/reverse/tda_reduced.py for the
same zero-variance trap and its fix, applied here too):
  log10(screening_suppression_factor), log10(phi_center_ratio),
  gamma_theta[0..14] (15 PTA angular bins), pta_max_deviation_from_hd
  = 18 components total.

Prediction (stated BEFORE running, per advisor review): mu_sym only
enters the screening block; c4_pta_product only enters PTA, and does so
LINEARLY (gamma_theta = hd_curve + c4_pta_product * l4_response, exact
per workshopcosmo.py:766) -- so all 15 PTA bins + max_deviation should be
EXACTLY COLLINEAR in that one direction. Expected Jacobian rank = 2 (one
direction per parameter, no internal degeneracy), i.e.
effective_dimension = 2 at every point.

Writes: jacobian_report.json
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROUND2_REVERSE_DIR = os.path.abspath(os.path.join(HERE, "..", "round2", "reverse"))
sys.path.insert(0, ROUND2_REVERSE_DIR)

import reduced_model as rm  # noqa: E402

PARAM_NAMES = ["mu_sym", "c4_pta_product"]
H_STEP_LN = 1e-2
N_PTA_BINS = 15


def observable_vector(params):
    r = rm.evaluate_reduced_point(params)
    sc, pta = r["screening"], r["pta"]
    ssf, pcr = sc["screening_suppression_factor"], sc["phi_center_ratio"]
    vec = np.array(
        [math.log10(ssf) if ssf > 0 and np.isfinite(ssf) else float("nan"),
         math.log10(pcr) if pcr > 0 and np.isfinite(pcr) else float("nan")]
        + list(pta["gamma_theta"])
        + [pta["max_deviation_from_hd"]],
        dtype=float,
    )
    stable = bool(sc["numerically_stable"]) and np.all(np.isfinite(vec))
    return vec, stable


OBS_NAMES = ["log10_ssf", "log10_pcr"] + [f"gamma_theta_{i}" for i in range(N_PTA_BINS)] + ["pta_max_dev"]


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
        J[:, j] = ((v_plus - sweep_mean) / sweep_std - (v_minus - sweep_mean) / sweep_std) / (2.0 * H_STEP_LN)
    return J, base_vec, True


def run():
    df = pd.read_csv(os.path.join(HERE, "sweep.csv"))
    stable = df[df["numerically_stable"] == True].copy()  # noqa: E712
    print(f"[round3 jacobian] stable sweep rows: {len(stable)} / {len(df)}")

    # standardization stats from the stable sweep (for the standardized
    # observable vector used by both Jacobian and, later, TDA)
    vecs, oks = [], []
    for _, row in stable.iterrows():
        v, ok = observable_vector({"mu_sym": float(row["mu_sym"]), "c4_pta_product": float(row["c4_pta_product"])})
        vecs.append(v); oks.append(ok)
    V = np.array(vecs)
    ok_mask = np.array(oks) & np.all(np.isfinite(V), axis=1)
    V = V[ok_mask]
    sweep_mean = V.mean(axis=0)
    sweep_std = V.std(axis=0)
    zero_var = sweep_std < 1e-14
    n_zero_var = int(np.sum(zero_var))
    sweep_std_safe = np.where(zero_var, 1.0, sweep_std)

    best_params = {"mu_sym": 1.0, "c4_pta_product": 0.08035}

    rng = np.random.RandomState(42)
    idxs = stable["idx"].to_numpy()
    random_idxs = rng.choice(idxs, size=min(5, len(idxs)), replace=False)

    def point_report(label, params):
        J, base_vec, ok = jacobian_at(params, sweep_std_safe, sweep_mean)
        if not ok or J is None:
            return {"label": label, "params": params, "ok": False}
        finite_cols = np.all(np.isfinite(J), axis=0)
        Jf = J[:, finite_cols] if np.any(~finite_cols) else J
        U, S, Vt = np.linalg.svd(Jf, full_matrices=False)
        thresh = 1e-3 * S[0] if len(S) and S[0] > 0 else 0.0
        eff_dim = int(np.sum(S > thresh)) if len(S) else 0
        near_null = []
        for k in range(Vt.shape[0]):
            if k >= len(S) or S[k] <= thresh:
                loadings = {PARAM_NAMES[j]: float(Vt[k, j]) for j in range(Vt.shape[1])}
                near_null.append({"singular_value": float(S[k]) if k < len(S) else 0.0, "right_singular_vector": loadings})
        return {
            "label": label, "params": params, "ok": True,
            "singular_values": S.tolist(),
            "effective_dimension_1e-3_of_max": eff_dim,
            "near_null_right_singular_vectors": near_null,
            "n_finite_param_columns": int(np.sum(finite_cols)),
        }

    reports = [point_report("nominal_flat_chi2_point", best_params)]
    for idx in random_idxs:
        row = stable[stable["idx"] == idx].iloc[0]
        p = {"mu_sym": float(row["mu_sym"]), "c4_pta_product": float(row["c4_pta_product"])}
        reports.append(point_report(f"random_idx_{idx}", p))

    eff_dims = [r["effective_dimension_1e-3_of_max"] for r in reports if r.get("ok")]

    out = {
        "generated": "2026-09-18",
        "method": "central difference, step=1e-2 in ln-param, standardized 18-component observable vector (log10_ssf, log10_pcr, gamma_theta x15, pta_max_dev), dark-energy block dropped (zero-variance by construction, frozen LCDM)",
        "n_zero_variance_observable_columns_before_dropping_de": n_zero_var,
        "n_stable_points_used_for_standardization": int(V.shape[0]),
        "param_names": PARAM_NAMES,
        "points": reports,
        "effective_dimension_summary": {"best_fit": eff_dims[0] if eff_dims else None, "random": eff_dims[1:]},
        "prediction_check": {
            "predicted_effective_dimension": 2,
            "observed": eff_dims,
            "matches_prediction": bool(all(d == 2 for d in eff_dims)) if eff_dims else False,
        },
    }
    with open(os.path.join(HERE, "jacobian_report.json"), "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(out["effective_dimension_summary"], indent=2))
    print(json.dumps(out["prediction_check"], indent=2))


if __name__ == "__main__":
    run()
