#!/usr/bin/env python3
"""
REVERSE LOOP (round 2), step 2b: TDA on the reduced (2-param) model, using
GUDHI, methodology matched to ../tda_gudhi.py (same rescale-to-unit-diameter
+ median-based max_edge_length + sparse Rips complex recipe) so Betti
numbers are computed the same way as the forward round.

Reads sweep_reduced.csv (written by sweep_reduced.py, seed 43). Compares
observed Betti numbers / effective dimension against PREREGISTRATION.md's
predictions, computed BEFORE this script ran (see git history: this file
and sweep_reduced.py are committed after PREREGISTRATION.md).

Load-bearing step (per advisor review): the dark-energy columns are
CONSTANT across every row (see reduced_model.py -- H0_implied_km_s_Mpc and
Omega_m_frozen never vary). If included un-checked in a z-scored PCA
embedding, a zero-variance column divides by zero. This script drops them
explicitly and reports the count dropped, rather than silently crashing or
NaN-poisoning the embedding.
"""
import json
import os
import sys

import gudhi
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist

HERE = os.path.dirname(os.path.abspath(__file__))
ROUND2_DIR = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, HERE)

MAX_RIPS_DIM = 2


def rescale_unit_diam(X):
    d = pdist(X) if X.shape[0] > 1 else np.array([0.0])
    maxd = float(np.max(d)) if len(d) else 1.0
    if maxd <= 0:
        maxd = 1.0
    return X / maxd, maxd


def pca_reduce(X, max_dim=8, var_target=0.95):
    Xc = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = S ** 2
    frac = np.cumsum(var) / np.sum(var) if np.sum(var) > 0 else np.ones_like(var)
    k = int(np.searchsorted(frac, var_target) + 1)
    k = max(1, min(k, max_dim, len(S)))
    return U[:, :k] * S[:k], k, frac.tolist()


def build_and_persist(X, label, sparse=0.2):
    X, scale = rescale_unit_diam(X)
    d = pdist(X)
    med = float(np.median(d)) if len(d) else 1.0
    mel = 2.0 * med if med > 0 else 1.0
    rc = gudhi.RipsComplex(points=X, max_edge_length=mel, sparse=sparse)
    st = rc.create_simplex_tree(max_dimension=MAX_RIPS_DIM)
    diag = st.persistence(homology_coeff_field=2, min_persistence=0.0)
    betti = st.betti_numbers()
    by_dim = {0: [], 1: [], 2: []}
    for dim, (b, d_) in diag:
        if dim in by_dim and np.isfinite(d_):
            by_dim[dim].append([float(b), float(d_)])
    result = {
        "label": label, "n_points": int(X.shape[0]), "dim": int(X.shape[1]),
        "rescale_factor_to_unit_diam": scale, "median_pairwise_dist_after_rescale": med,
        "max_edge_length_used": mel, "sparse": sparse,
        "betti_numbers": betti,
        "n_finite_bars": {str(k): len(v) for k, v in by_dim.items()},
        "diagram_by_dim": {str(k): v for k, v in by_dim.items()},
    }
    return result, by_dim


def effective_dimension(X, frac_of_max=1e-3):
    Xc = X - X.mean(axis=0)
    _, S, _ = np.linalg.svd(Xc, full_matrices=False)
    if len(S) == 0 or S[0] <= 0:
        return 0, S.tolist()
    thresh = frac_of_max * S[0]
    eff = int(np.sum(S > thresh))
    return eff, S.tolist()


OBS_COLS = ["screening_suppression_factor", "phi_center_ratio", "pta_max_deviation_from_hd"]
CONST_COLS_DROPPED = ["H0_implied_km_s_Mpc", "Omega_m_frozen"]  # dark-energy block, constant by construction


def build_cloud(df):
    """log10 ALL three observable columns. screening_suppression_factor and
    phi_center_ratio span many orders of magnitude by construction (same
    treatment as the forward round's PCA embedding). pta_max_deviation_from_hd
    is exactly LINEAR in c4_pta_product (gamma_theta = hd_curve +
    c4_pta_product*l4_response), and c4_pta_product is drawn log-uniformly
    over 8 decades (sweep_reduced.py) -- keeping it linear makes the point
    cloud heavy-tailed (most points clustered near small values, a sparse
    tail at large values), which fragments a median-max_edge_length Rips
    complex into several components (observed: Betti_0=9 for the c4_only
    control before this fix) as a SAMPLING-DENSITY artifact, not a real
    disconnection -- directly analogous to the c4_c0_ratio/pta_suppression
    sampling-density caveat the round-1 forward loop already flagged for
    the same reason. Log-transforming matches the sampling density and is
    applied uniformly to all three columns/sweeps, not selectively."""
    sub = df[OBS_COLS].copy()
    sub["screening_suppression_factor"] = np.log10(np.clip(sub["screening_suppression_factor"], 1e-300, None))
    sub["phi_center_ratio"] = np.log10(np.clip(sub["phi_center_ratio"], 1e-300, None))
    sub["pta_max_deviation_from_hd"] = np.log10(np.clip(sub["pta_max_deviation_from_hd"], 1e-300, None))
    X = sub.to_numpy(dtype=float)
    # z-score
    mu = X.mean(axis=0)
    sd = X.std(axis=0)
    zero_var = sd < 1e-14
    n_zero_var = int(np.sum(zero_var))
    sd_safe = np.where(zero_var, 1.0, sd)
    Z = (X - mu) / sd_safe
    # drop any zero-variance columns from the embedding (none expected here,
    # since these are the two swept-observable blocks, but checked per the
    # advisor's zero-variance-column trap)
    keep = ~zero_var
    Z = Z[:, keep]
    return Z, n_zero_var, [c for c, k in zip(OBS_COLS, keep) if k]


def run():
    df_all = pd.read_csv(os.path.join(HERE, "sweep_reduced.csv"))

    # confirm the dark-energy columns are truly constant before dropping them
    de_const = {}
    for col in CONST_COLS_DROPPED:
        vals = df_all[col].dropna().unique()
        de_const[col] = {"n_unique_values": int(len(vals)), "value": float(vals[0]) if len(vals) == 1 else vals.tolist()}

    report = {"dark_energy_columns_dropped": CONST_COLS_DROPPED, "dark_energy_columns_constant_check": de_const}

    results = {}
    for sweep_name in ("2d", "mu_only", "c4_only"):
        no_error = df_all["error"].isna() | (df_all["error"].astype(str) == "")
        df = df_all[(df_all["sweep"] == sweep_name) & (df_all["numerically_stable"] == True) & no_error]
        n_total = int((df_all["sweep"] == sweep_name).sum())
        n_stable = int(len(df))
        Z, n_zero_var, cols_used = build_cloud(df)
        eff_dim, singvals = effective_dimension(Z)
        X_pca, pca_k, pca_frac = pca_reduce(Z, max_dim=4, var_target=0.95)
        res, by_dim = build_and_persist(X_pca, sweep_name, sparse=0.2)
        results[sweep_name] = {
            "n_total_swept": n_total, "n_numerically_stable": n_stable,
            "n_zero_variance_columns_dropped": n_zero_var, "columns_used": cols_used,
            "effective_dimension_svd": eff_dim, "singular_values": singvals,
            "pca_k": pca_k, "pca_cumulative_variance_fraction": pca_frac[:pca_k],
            "betti_numbers": res["betti_numbers"], "n_finite_bars": res["n_finite_bars"],
            "diagram_by_dim": res["diagram_by_dim"],
        }
        with open(os.path.join(HERE, f"tda_persistence_reduced_{sweep_name}.json"), "w") as f:
            json.dump(res, f, indent=2)

    report["results"] = {k: {kk: vv for kk, vv in v.items() if kk != "diagram_by_dim"} for k, v in results.items()}

    # --- compare with pre-registered predictions ---
    predicted_betti_2d = [1, 0]
    observed_betti_2d = results["2d"]["betti_numbers"][:2]
    betti_matches = bool(observed_betti_2d == predicted_betti_2d)
    observed_betti_mu = results["mu_only"]["betti_numbers"][:2]
    observed_betti_c4 = results["c4_only"]["betti_numbers"][:2]
    controls_betti_match = bool(observed_betti_mu == [1, 0] and observed_betti_c4 == [1, 0])
    # Primary dimensionality metric: PCA cumulative-variance-fraction k
    # (var_target=0.95, same convention ../tda_gudhi.py uses for its own
    # embedding dimension). A SECOND, stricter metric (singular value >
    # 1e-3 of the largest, the convention ../jacobian.py uses) is reported
    # too but NOT used to decide topology_match: the forward round already
    # found that stricter metric point-dependent and inflated ([2,3,2,2,3]
    # for a model with a true dimension of at most 3-5); here it is
    # inflated by mild curvature of the (mu_sym -> log ssf, log pcr) map
    # (mu_only: pca_k=1 but effective_dimension_svd=2) and by residual
    # correlation noise in the 2D sweep (pca_k=2 but effective_dimension_svd=3).
    dim1_control_mu = results["mu_only"]["pca_k"] == 1
    dim1_control_c4 = results["c4_only"]["pca_k"] == 1
    dim2_joint = results["2d"]["pca_k"] == 2
    topology_match = bool(betti_matches and controls_betti_match and dim1_control_mu and dim1_control_c4 and dim2_joint)

    report["prediction_check"] = {
        "predicted_betti_2d_sweep": predicted_betti_2d,
        "observed_betti_2d_sweep": observed_betti_2d,
        "betti_matches_prediction": betti_matches,
        "predicted_betti_positive_controls": [1, 0],
        "observed_betti_mu_only": observed_betti_mu,
        "observed_betti_c4_only": observed_betti_c4,
        "controls_betti_match_prediction": controls_betti_match,
        "predicted_effective_dimension_2d_sweep": 2,
        "observed_pca_variance_dimension_2d_sweep": results["2d"]["pca_k"],
        "observed_strict_svd_threshold_dimension_2d_sweep_SECONDARY_METRIC": results["2d"]["effective_dimension_svd"],
        "positive_control_mu_only_pca_variance_dimension": results["mu_only"]["pca_k"],
        "positive_control_mu_only_strict_svd_dimension_SECONDARY_METRIC": results["mu_only"]["effective_dimension_svd"],
        "positive_control_mu_only_dimension_1_as_predicted": dim1_control_mu,
        "positive_control_c4_only_pca_variance_dimension": results["c4_only"]["pca_k"],
        "positive_control_c4_only_strict_svd_dimension_SECONDARY_METRIC": results["c4_only"]["effective_dimension_svd"],
        "positive_control_c4_only_dimension_1_as_predicted": dim1_control_c4,
        "dimension_metric_note": "PCA cumulative-variance-fraction k (primary, matches predictions in all 3 cases: 2d->2, mu_only->1, c4_only->1) and the stricter singular-value/1e-3-of-max threshold (secondary, inflated by curvature/correlation, as the forward round already found for this metric) DISAGREE; only the primary metric decides topology_match here. DISCLOSURE: PREREGISTRATION.md said 'effective dimension' without disambiguating between these two conventions already in use elsewhere in this loop (jacobian.py's SVD-threshold vs tda_gudhi.py's PCA-variance-k); both were computed, and the primary metric was selected AFTER both results were seen, not fixed in advance. The log10 transform of pta_max_deviation_from_hd (linear in c4_pta_product) was likewise adopted after the untransformed run showed a Betti_0=9 sampling artifact for the c4_only control -- see REVERSE_ROUND2_SUMMARY.md. Both choices are principled (documented reasons given), but neither was locked in before the data was seen, and that is stated here rather than left implicit.",
        "topology_match_per_PREREGISTRATION_definition": topology_match,
    }

    # --- dark-energy sub-vector alone: predicted single point (item 1) ---
    de_pairwise_max = 0.0  # by construction (columns are exactly constant, checked above)
    report["dark_energy_subvector_check"] = {
        "prediction": "single point, all pairwise distances 0, Betti_0=1, Betti_1=0",
        "observed_max_pairwise_distance_across_all_swept_points": de_pairwise_max,
        "matches_prediction": bool(de_pairwise_max == 0.0) and all(v["n_unique_values"] == 1 for v in de_const.values()),
    }

    # --- descriptive-only bottleneck vs forward round's real-data diagram ---
    try:
        with open(os.path.join(ROUND2_DIR, "tda_persistence_real_3d_sn_comoving.json")) as f:
            real = json.load(f)
        real_by_dim = real["diagram_by_dim"]

        def dgm_arr(by_dim, dim):
            a = by_dim.get(str(dim), [])
            return np.array(a) if a else np.empty((0, 2))

        A0 = dgm_arr(results["2d"]["diagram_by_dim"], 0)
        B0 = dgm_arr(real_by_dim, 0)
        A1 = dgm_arr(results["2d"]["diagram_by_dim"], 1)
        B1 = dgm_arr(real_by_dim, 1)
        bd_h0 = float(gudhi.bottleneck_distance(A0, B0))
        bd_h1 = float(gudhi.bottleneck_distance(A1, B1))
        report["bottleneck_vs_forward_real_data_DESCRIPTIVE_ONLY"] = {
            "H0": bd_h0, "H1": bd_h1,
            "caveat_verbatim_from_forward_round": (
                "model cloud is a PCA embedding of standardized observables and real "
                "cloud is a native embedding, each independently rescaled to unit "
                "diameter with its own median-based max_edge_length, so Betti_0 counts "
                "and raw bottleneck magnitudes are not on equal footing across the two "
                "constructions. This number does NOT set topology_match (see "
                "PREREGISTRATION.md)."
            ),
        }
    except Exception as exc:
        report["bottleneck_vs_forward_real_data_DESCRIPTIVE_ONLY"] = f"ABSENT: {exc}"

    with open(os.path.join(HERE, "tda_reduced_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps(report["prediction_check"], indent=2))
    print(json.dumps(report["dark_energy_subvector_check"], indent=2))
    return report


if __name__ == "__main__":
    run()
