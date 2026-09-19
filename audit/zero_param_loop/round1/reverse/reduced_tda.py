#!/usr/bin/env python3
"""
REVERSE LOOP step 2 (TDA, real INRIA GUDHI `import gudhi`): persistent
homology of the REDUCED 5-param model's observable-image point cloud, seed
43, compared against the pre-registered prediction in tda_prediction.json
(written before this script was run) and against round 1's real-data
diagrams and model diagram.

Fixes applied per pre-registration / advisor review of round 1:
  - Model cloud subsampled to EXACTLY N=316 (round 1's N) before building
    the Rips complex, so H1-bar-count / bottleneck comparisons are not an
    N artifact.
  - The null control is VARIANCE-MATCHED to the model cloud's own per-axis
    (post-PCA) standard deviations, not a uniform unit box -- round 1
    flagged the box-null ordering (null_vs_real_sn_H1 > model_vs_real_sn_H1)
    as an unresolved possible bounding-box artifact; this control is what
    tests that.

Writes: reduced_tda_persistence_model.json/.png,
        reduced_tda_persistence_real_sn.json/.png,
        reduced_tda_persistence_null_variance_matched.json/.png,
        reduced_tda_persistence_null_uniform_box.json/.png (round-1-style,
          kept for direct comparability),
        reduced_tda_persistence_circle.json/.png,
        reduced_tda_report.json
"""
import json
import math
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gudhi

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "..", ".."))
ROUND1_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, REPO_ROOT)

SEED = 43
rng = np.random.RandomState(SEED)
MAX_RIPS_DIM = 2
N_TARGET = 316  # match round-1's model cloud N exactly


def rescale_unit_diam(X):
    from scipy.spatial.distance import pdist
    d = pdist(X)
    maxd = float(np.max(d)) if len(d) else 1.0
    if maxd <= 0:
        maxd = 1.0
    return X / maxd, maxd


def pca_reduce(X, max_dim=8, var_target=0.95):
    Xc = X - X.mean(axis=0)
    U, S, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = S ** 2
    frac = np.cumsum(var) / np.sum(var)
    k = int(np.searchsorted(frac, var_target) + 1)
    k = max(1, min(k, max_dim, len(S)))
    return U[:, :k] * S[:k], k, frac.tolist()


def build_and_persist(X, label, sparse=0.2):
    from scipy.spatial.distance import pdist
    X, scale = rescale_unit_diam(X)
    d = pdist(X)
    med = float(np.median(d))
    mel = 2.0 * med
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
    fig, ax = plt.subplots(figsize=(5, 5))
    colors = {0: "tab:blue", 1: "tab:orange", 2: "tab:green"}
    for dim in (0, 1, 2):
        pts = by_dim[dim]
        if pts:
            arr = np.array(pts)
            ax.scatter(arr[:, 0], arr[:, 1], s=10, color=colors[dim], label=f"H{dim} (n={len(pts)})")
    lim = max([p[1] for v in by_dim.values() for p in v] + [1e-6]) * 1.1
    ax.plot([0, lim], [0, lim], "k--", lw=0.5)
    ax.set_xlabel("birth"); ax.set_ylabel("death"); ax.set_title(f"Reduced-model persistence diagram: {label}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, f"reduced_tda_persistence_{label}.png"), dpi=120)
    plt.close(fig)
    with open(os.path.join(HERE, f"reduced_tda_persistence_{label}.json"), "w") as f:
        json.dump(result, f, indent=2)
    return result, by_dim


# ---------------------------------------------------------------------------
# Model cloud from reduced_sweep.csv (5-param model), SAME 20-dim raw
# observable layout as round-1's build_model_cloud (4 non-DESI + 7*2 DESI +
# 2 log-screening = 20), except column 4 is log10(c4_pta_product) directly
# (no product-of-two-params reconstruction needed -- it's the parameter).
# ---------------------------------------------------------------------------
def build_model_cloud():
    df = pd.read_csv(os.path.join(HERE, "reduced_sweep.csv"))
    clean = df[df["error"].isna()]
    stable = clean[clean["numerically_stable"] == True].copy()  # noqa: E712
    desi_z = [0.295, 0.510, 0.706, 0.930, 1.317, 1.491, 2.330]
    cols = ["w0_cpl_latetime", "wa_cpl_latetime", "max_deviation_from_hd", "c4_pta_product"]
    X = stable[cols].to_numpy(dtype=float)
    for z in desi_z:
        X = np.column_stack([X, stable[f"DM_H0_z{z}"].to_numpy(dtype=float),
                              stable[f"HzH0_z{z}"].to_numpy(dtype=float)])
    ssf = stable["screening_suppression_factor"].to_numpy(dtype=float)
    pcr = stable["phi_center_ratio"].to_numpy(dtype=float)
    log_ssf = np.log10(np.clip(ssf, 1e-300, None))
    log_pcr = np.log10(np.clip(pcr, 1e-300, None))
    log_c4pta = np.log10(np.clip(X[:, 3], 1e-300, None))
    X[:, 3] = log_c4pta
    X = np.column_stack([X, log_ssf, log_pcr])
    n_before_subsample = X.shape[0]
    if X.shape[0] > N_TARGET:
        idx = rng.choice(X.shape[0], size=N_TARGET, replace=False)
        X = X[idx]
    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd < 1e-12] = 1.0
    Xz = (X - mu) / sd
    return Xz, n_before_subsample


def build_real_sn_cloud(n_target):
    """Reuses this round's OWN independently-recomputed LCDM (Om,offset)
    from reduced_chi2_report.json (verified consistent with round 1's to
    abs_diff=0.0), not round 1's file."""
    path = os.path.join(REPO_ROOT, "data", "real", "dark_energy", "pantheon_plus_sh0es.dat")
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df = df[(df["zHD"] > 0.01) & (df["zHD"] <= 2.4) & np.isfinite(df["MU_SH0ES"])
            & np.isfinite(df["c"]) & np.isfinite(df["x1"])].copy()
    with open(os.path.join(HERE, "reduced_chi2_report.json")) as f:
        chi2_rep = json.load(f)
    Om = chi2_rep["chi2_lcdm_independent_recompute"]["sn"]["Om_fit"]
    offset = chi2_rep["chi2_lcdm_independent_recompute"]["sn"]["offset_fit"]
    zg = np.linspace(0.0, 2.45, 4000)
    Ez = np.sqrt(Om * (1.0 + zg) ** 3 + (1.0 - Om))
    inv_E = 1.0 / Ez
    DC = np.concatenate(([0.0], np.cumsum(0.5 * (inv_E[1:] + inv_E[:-1]) * np.diff(zg))))
    DC_q = np.interp(df["zHD"].to_numpy(), zg, DC)
    DL_q = (1.0 + df["zHD"].to_numpy()) * DC_q
    mu_lcdm = 5.0 * np.log10(np.clip(DL_q, 1e-12, None)) + offset
    resid = df["MU_SH0ES"].to_numpy() - mu_lcdm
    X = np.column_stack([df["zHD"].to_numpy(), df["c"].to_numpy(), df["x1"].to_numpy(), resid])
    if len(X) > n_target:
        idx = rng.choice(len(X), size=n_target, replace=False)
        X = X[idx]
    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd < 1e-12] = 1.0
    return (X - mu) / sd, Om, offset


def build_circle_cloud(n, noise_sigma=0.05):
    theta = rng.uniform(0, 2 * math.pi, size=n)
    r = 1.0 + rng.normal(0, noise_sigma, size=n)
    return np.column_stack([r * np.cos(theta), r * np.sin(theta)])


def build_null_uniform_box(n, dim, box_halfwidth=1.0):
    """Round-1-style null, kept for direct comparability."""
    return rng.uniform(-box_halfwidth, box_halfwidth, size=(n, dim))


def build_null_variance_matched(n, model_cloud_post_pca):
    """NEW this round (advisor-recommended fix to round-1's flagged, unresolved
    caveat): per-axis-variance-matched null. model_cloud_post_pca is the
    PCA-reduced (NOT yet standardized-to-unit-variance) model cloud, whose
    k columns have decreasing variance (PCA orders by singular value). A
    Uniform(-a,a) distribution has variance a^2/3, so a_i = sqrt(3*var_i)
    reproduces each axis's actual variance exactly, unlike a fixed-halfwidth
    box which imposes EQUAL variance on every axis regardless of the real
    PCA spectrum."""
    dim = model_cloud_post_pca.shape[1]
    stds = model_cloud_post_pca.std(axis=0)
    half_widths = np.sqrt(3.0) * stds
    out = np.zeros((n, dim))
    for i in range(dim):
        out[:, i] = rng.uniform(-half_widths[i], half_widths[i], size=n)
    return out


def dgm_array(by_dim, dim):
    arr = by_dim.get(str(dim)) if isinstance(by_dim, dict) and str(dim) in by_dim else by_dim.get(dim, [])
    return np.array(arr) if arr else np.empty((0, 2))


def bd(a, b, dim):
    A, B = dgm_array(a, dim), dgm_array(b, dim)
    try:
        return float(gudhi.bottleneck_distance(A, B))
    except Exception as exc:  # pragma: no cover
        return f"ERROR: {exc}"


def run():
    with open(os.path.join(HERE, "tda_prediction.json")) as f:
        prediction = json.load(f)

    X_model_raw, n_before_subsample = build_model_cloud()
    X_model_pca, pca_k, pca_frac = pca_reduce(X_model_raw, max_dim=8, var_target=0.95)
    res_model, bd_model = build_and_persist(X_model_pca, "model", sparse=0.2)

    X_sn, Om_used, offset_used = build_real_sn_cloud(n_target=res_model["n_points"])
    res_sn, bd_sn = build_and_persist(X_sn, "real_sn", sparse=0.2)

    X_circle = build_circle_cloud(res_model["n_points"])
    res_circle, bd_circle = build_and_persist(X_circle, "circle", sparse=0.2)
    h1_circle = dgm_array(bd_circle, 1)
    if len(h1_circle) >= 1:
        pers = np.sort(h1_circle[:, 1] - h1_circle[:, 0])[::-1]
        top = float(pers[0]); second = float(pers[1]) if len(pers) > 1 else 0.0
        ratio = top / second if second > 1e-12 else float("inf")
        circle_ok = bool(second == 0.0 or ratio > 5.0)
    else:
        top = second = ratio = 0.0; circle_ok = False

    X_null_box = build_null_uniform_box(res_model["n_points"], dim=pca_k)
    res_null_box, bd_null_box = build_and_persist(X_null_box, "null_uniform_box", sparse=0.2)

    X_null_vm = build_null_variance_matched(res_model["n_points"], X_model_pca)
    res_null_vm, bd_null_vm = build_and_persist(X_null_vm, "null_variance_matched", sparse=0.2)

    # load round-1's own model diagram for the direct round1-vs-round2 bottleneck check
    with open(os.path.join(ROUND1_DIR, "tda_persistence_model.json")) as f:
        r1_model = json.load(f)
    bd_r1_model = r1_model["diagram_by_dim"]

    effective_dim_reduced = None
    try:
        with open(os.path.join(HERE, "reduced_jacobian_report.json")) as f:
            jr = json.load(f)
        effective_dim_reduced = jr["probed_points"]["best_fit"]["effective_dimension_at_1e-3_of_max"]
    except Exception as exc:
        effective_dim_reduced = f"ABSENT: could not read reduced_jacobian_report.json ({exc})"

    betti_model = res_model["betti_numbers"]
    predicted_betti = [1, 0]
    betti_matches_corrected_prediction = bool(betti_model[:2] == predicted_betti)

    report = {
        "prediction_file": "tda_prediction.json (written before this script ran)",
        "model": {k: v for k, v in res_model.items() if k != "diagram_by_dim"},
        "model_n_before_subsample": n_before_subsample,
        "model_pca_raw_dim": int(X_model_raw.shape[1]),
        "model_pca_reduced_dim": pca_k,
        "model_pca_cumulative_variance_fraction": pca_frac[:pca_k],
        "real_sn": {k: v for k, v in res_sn.items() if k != "diagram_by_dim"},
        "real_sn_lcdm_Om_used": Om_used, "real_sn_lcdm_offset_used": offset_used,
        "circle_control": {k: v for k, v in res_circle.items() if k != "diagram_by_dim"},
        "circle_control_assertion": {
            "n_H1_bars": int(len(h1_circle)), "top_persistence": top, "second_persistence": second,
            "ratio": ratio, "ASSERTION_one_dominant_H1_bar": circle_ok,
        },
        "null_uniform_box_control": {k: v for k, v in res_null_box.items() if k != "diagram_by_dim"},
        "null_variance_matched_control": {k: v for k, v in res_null_vm.items() if k != "diagram_by_dim"},
        "bottleneck_distances": {
            "model_vs_real_sn_H0": bd(bd_model, bd_sn, 0),
            "model_vs_real_sn_H1": bd(bd_model, bd_sn, 1),
            "null_uniform_box_vs_real_sn_H1": bd(bd_null_box, bd_sn, 1),
            "null_variance_matched_vs_real_sn_H1": bd(bd_null_vm, bd_sn, 1),
            "model_reduced_vs_model_round1_H0": bd(bd_model, bd_r1_model, 0),
            "model_reduced_vs_model_round1_H1": bd(bd_model, bd_r1_model, 1),
        },
        "variance_matched_null_check": {
            "purpose": "round-1 next_round_priorities #3: was null_vs_real_sn_H1 > model_vs_real_sn_H1 a bounding-box artifact of the uniform-box null? This control uses a null matched to the model cloud's own per-axis variance instead.",
        },
        "predicted_topology_naive": prediction["naive_prediction_as_given"],
        "predicted_topology_corrected": prediction["corrected_prediction"],
        "betti_numbers_observed": betti_model,
        "betti_matches_corrected_prediction": betti_matches_corrected_prediction,
        "effective_dimension_svd_reduced_model": effective_dim_reduced,
        "effective_dimension_svd_round1_full_model": 3,
        "effective_dimension_matches_prediction": (effective_dim_reduced == 3),
    }

    vm_h1 = report["bottleneck_distances"]["null_variance_matched_vs_real_sn_H1"]
    box_h1 = report["bottleneck_distances"]["null_uniform_box_vs_real_sn_H1"]
    model_h1 = report["bottleneck_distances"]["model_vs_real_sn_H1"]
    if isinstance(vm_h1, float) and isinstance(model_h1, float):
        report["variance_matched_null_check"]["null_variance_matched_vs_real_sn_H1"] = vm_h1
        report["variance_matched_null_check"]["null_uniform_box_vs_real_sn_H1"] = box_h1
        report["variance_matched_null_check"]["model_vs_real_sn_H1"] = model_h1
        ordering_survives = bool(vm_h1 > model_h1)
        report["variance_matched_null_check"]["ordering_null_farther_than_model_survives_variance_matching"] = ordering_survives
        report["variance_matched_null_check"]["conclusion"] = (
            "Ordering (null farther from real-SN topology than the model is) SURVIVES variance-matching -- "
            "not purely a bounding-box artifact." if ordering_survives else
            "Ordering does NOT survive variance-matching -- round-1's null_vs_real_sn_H1 > model_vs_real_sn_H1 "
            "finding WAS a bounding-box artifact; no model-vs-real topology preference is supported."
        )

    with open(os.path.join(HERE, "reduced_tda_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: v for k, v in report.items()
                       if k not in ("model", "real_sn", "circle_control", "null_uniform_box_control", "null_variance_matched_control")},
                      indent=2))
    return report


if __name__ == "__main__":
    run()
