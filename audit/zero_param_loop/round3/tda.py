#!/usr/bin/env python3
"""
ROUND 3 forward loop, step 3: TDA with real INRIA GUDHI (`import gudhi`,
not scripts/tda_mapper.py). Same recipe as round2/tda_gudhi.py and
round2/reverse/tda_reduced.py (rescale-to-unit-diameter, median-based
max_edge_length, sparse Rips complex, homology dims 0-2), reused via
import rather than reimplemented.

Model side: standardized OBSERVABLE-space cloud built from round3/
sweep.csv's stable rows (18 columns: log10_ssf, log10_pcr, gamma_theta
x15, pta_max_dev), PCA-reduced to 95% variance.

Real-data side: the SAME genuine 3D comoving Pantheon+ SH0ES cloud
(RA/DEC/zHD, flat LCDM Om=0.31115) round 2 built, recomputed fresh here
via round2/tda_gudhi.build_real_3d_cloud (imported, not copied) -- this
satisfies the task's "3D comoving cloud with redshifts" requirement
without a new fetch (round 1's 2MRS catalogue had no z; this Pantheon+
cloud does).

Controls: known-answer noisy circle (round2's build_circle_cloud) and a
null Poisson cloud matched in N and bounding volume to the real cloud
(round2's build_null_uniform_box_matched). Both reused unmodified.

PREDICTION (stated before running, per advisor review): c4_pta_product
enters PTA linearly (all 15 gamma_theta bins move together along one
direction); mu_sym enters screening nonlinearly but through only 2
scalars. With 2 free parameters and no forcing to fill out any 3rd
independent direction, the model's observable cloud is expected to be a
2-D SHEET: Betti (1, 0), matching jacobian.py's effective_dimension=2 at
most points.

KNOWN BOTTLENECK-DEGENERACY, already diagnosed in round2/ledger.json
("round1_flag_unverified_resolution") and round2/reverse/tda_reduced.py's
descriptive-only comparison: gudhi.bottleneck_distance sends a diagram's
whole dominant unmatched bar to the diagonal (distance = half that bar's
persistence) whenever it exceeds twice the OTHER diagram's max
persistence in that dimension, regardless of what else is in the other
diagram. This round explicitly checks each reported bottleneck distance
against that formula and labels it when triggered, instead of
re-diagnosing it as if new.

Writes: tda_persistence_model.json/.png, tda_persistence_real_3d.json/.png,
        tda_persistence_circle.json/.png, tda_persistence_null_poisson.json/.png,
        tda_report.json
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import gudhi

HERE = os.path.dirname(os.path.abspath(__file__))
ROUND2_DIR = os.path.abspath(os.path.join(HERE, "..", "round2"))
sys.path.insert(0, ROUND2_DIR)

import tda_gudhi as tg  # noqa: E402 (round2 module: rescale_unit_diam, pca_reduce,
                          # build_real_3d_cloud, build_circle_cloud,
                          # build_null_uniform_box_matched, dgm_array, bd)

SEED = 42
rng = np.random.RandomState(SEED)
N_TARGET = 300


def build_model_cloud():
    df = pd.read_csv(os.path.join(HERE, "sweep.csv"))
    stable = df[df["numerically_stable"] == True].copy()  # noqa: E712
    cols = ["log10_ssf", "log10_pcr"] + [f"gamma_theta_{i}" for i in range(15)] + ["pta_max_deviation_from_hd"]
    X = stable[cols].to_numpy(dtype=float)
    finite_rows = np.all(np.isfinite(X), axis=1)
    X = X[finite_rows]
    n_before_subsample = X.shape[0]
    if X.shape[0] > N_TARGET:
        idx = rng.choice(X.shape[0], size=N_TARGET, replace=False)
        X = X[idx]
    mu, sd = X.mean(axis=0), X.std(axis=0)
    n_zero_var = int(np.sum(sd < 1e-12))
    sd_safe = np.where(sd < 1e-12, 1.0, sd)
    Xz = (X - mu) / sd_safe
    return Xz, n_before_subsample, n_zero_var, cols


def build_and_persist_local(X, label, sparse=0.2):
    """Local wrapper around tda_gudhi.build_and_persist that writes into
    THIS round's directory instead of round2's."""
    from scipy.spatial.distance import pdist
    X2, scale = tg.rescale_unit_diam(X)
    d = pdist(X2)
    med = float(np.median(d)) if len(d) else 1.0
    mel = 2.0 * med if med > 0 else 1.0
    rc = gudhi.RipsComplex(points=X2, max_edge_length=mel, sparse=sparse)
    st = rc.create_simplex_tree(max_dimension=2)
    diag = st.persistence(homology_coeff_field=2, min_persistence=0.0)
    betti = st.betti_numbers()
    by_dim = {0: [], 1: [], 2: []}
    for dim, (b, dth) in diag:
        if dim in by_dim and np.isfinite(dth):
            by_dim[dim].append([float(b), float(dth)])
    result = {
        "label": label, "n_points": int(X2.shape[0]), "dim": int(X2.shape[1]),
        "rescale_factor_to_unit_diam": scale, "median_pairwise_dist_after_rescale": med,
        "max_edge_length_used": mel, "sparse": sparse,
        "betti_numbers": betti,
        "n_finite_bars": {str(k): len(v) for k, v in by_dim.items()},
        "diagram_by_dim": {str(k): v for k, v in by_dim.items()},
    }
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(5, 5))
    colors = {0: "tab:blue", 1: "tab:orange", 2: "tab:green"}
    for dim in (0, 1, 2):
        pts = by_dim[dim]
        if pts:
            arr = np.array(pts)
            ax.scatter(arr[:, 0], arr[:, 1], s=10, color=colors[dim], label=f"H{dim} (n={len(pts)})")
    lim = max([p[1] for v in by_dim.values() for p in v] + [1e-6]) * 1.1
    ax.plot([0, lim], [0, lim], "k--", lw=0.5)
    ax.set_xlabel("birth"); ax.set_ylabel("death"); ax.set_title(f"round3 persistence: {label}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, f"tda_persistence_{label}.png"), dpi=120)
    plt.close(fig)
    with open(os.path.join(HERE, f"tda_persistence_{label}.json"), "w") as f:
        json.dump(result, f, indent=2)
    return result, by_dim


def check_half_bar_degeneracy(dist, dgm_a, dgm_b):
    """Returns True (with the triggering bar) if `dist` equals half the
    larger diagram's dominant unmatched-bar persistence -- the exact
    mechanism round2/ledger.json diagnosed, re-checked here per point
    rather than re-derived."""
    def max_pers(dgm):
        arr = np.asarray(dgm)
        if arr.size == 0:
            return 0.0
        return float(np.max(arr[:, 1] - arr[:, 0]))
    pa, pb = max_pers(dgm_a), max_pers(dgm_b)
    dominant = max(pa, pb)
    other = min(pa, pb)
    triggered = other > 0 and dominant > 2.0 * other and abs(dist - dominant / 2.0) < 1e-9
    return {
        "max_persistence_a": pa, "max_persistence_b": pb,
        "dominant_over_2": dominant / 2.0 if dominant else 0.0,
        "half_max_degeneracy_triggered": bool(triggered),
    }


def run():
    X_model, n_model_before, n_zero_var_model, model_cols = build_model_cloud()
    X_model_pca, pca_k, pca_frac = tg.pca_reduce(X_model, max_dim=8, var_target=0.95)
    res_model, bd_model = build_and_persist_local(X_model_pca, "model", sparse=0.2)

    X_real, n_real_before = tg.build_real_3d_cloud(n_target=res_model["n_points"])
    res_real, bd_real = build_and_persist_local(X_real, "real_3d_sn_comoving", sparse=0.2)

    X_circle = tg.build_circle_cloud(n=res_model["n_points"], noise_sigma=0.05)
    res_circle, bd_circle = build_and_persist_local(X_circle, "known_answer_circle", sparse=0.2)

    X_null = tg.build_null_uniform_box_matched(n=res_real["n_points"], X_ref=X_real)
    res_null, bd_null = build_and_persist_local(X_null, "null_poisson_matched_to_real", sparse=0.2)

    bd_h0_model_real = tg.bd(res_model["diagram_by_dim"], res_real["diagram_by_dim"], 0)
    bd_h1_model_real = tg.bd(res_model["diagram_by_dim"], res_real["diagram_by_dim"], 1)
    bd_h0_null_real = tg.bd(res_null["diagram_by_dim"], res_real["diagram_by_dim"], 0)
    bd_h1_null_real = tg.bd(res_null["diagram_by_dim"], res_real["diagram_by_dim"], 1)

    degeneracy_h1_model_real = check_half_bar_degeneracy(
        bd_h1_model_real, bd_model[1], bd_real[1])
    degeneracy_h1_null_real = check_half_bar_degeneracy(
        bd_h1_null_real, bd_null[1], bd_real[1])

    # circle self-check: one long H1 bar, ratio to next-longest
    circle_h1 = np.array(res_circle["diagram_by_dim"]["1"]) if res_circle["diagram_by_dim"]["1"] else np.empty((0, 2))
    if len(circle_h1) >= 1:
        pers = np.sort(circle_h1[:, 1] - circle_h1[:, 0])[::-1]
        circle_ratio = float(pers[0] / pers[1]) if len(pers) > 1 and pers[1] > 0 else float("inf")
    else:
        circle_ratio = None

    report = {
        "generated": "2026-09-18",
        "tool": f"gudhi {gudhi.__version__}",
        "model_cloud": {
            "n_stable_sweep_rows": n_model_before, "n_zero_variance_columns": n_zero_var_model,
            "columns": model_cols, "pca_k": pca_k, "pca_cumulative_variance_fraction_at_k": pca_frac[pca_k - 1] if pca_frac else None,
            "betti_numbers": res_model["betti_numbers"], "n_points_used": res_model["n_points"],
        },
        "real_cloud": {
            "source": "Pantheon+ SH0ES RA/DEC/zHD, flat LCDM Om=0.31115 comoving distance (sha256-verified, data/real/dark_energy/pantheon_plus_sh0es.dat)",
            "n_sn_before_subsample": n_real_before, "n_points_used": res_real["n_points"],
            "betti_numbers": res_real["betti_numbers"],
        },
        "known_answer_control_circle": {
            "n_points": res_circle["n_points"], "betti_numbers": res_circle["betti_numbers"],
            "n_H1_bars": res_circle["n_finite_bars"].get("1", 0),
            "dominant_to_second_H1_bar_ratio": circle_ratio,
            "passes_single_long_H1_bar_check": bool(circle_ratio is not None and circle_ratio > 5.0),
        },
        "null_control_poisson_matched_to_real": {
            "n_points": res_null["n_points"], "betti_numbers": res_null["betti_numbers"],
        },
        "bottleneck_distances": {
            "model_vs_real_H0": bd_h0_model_real, "model_vs_real_H1": bd_h1_model_real,
            "null_vs_real_H0": bd_h0_null_real, "null_vs_real_H1": bd_h1_null_real,
            "half_max_persistence_degeneracy_check_model_vs_real_H1": degeneracy_h1_model_real,
            "half_max_persistence_degeneracy_check_null_vs_real_H1": degeneracy_h1_null_real,
            "ordering_model_closer_than_null": (
                bd_h1_model_real < bd_h1_null_real if isinstance(bd_h1_model_real, float) and isinstance(bd_h1_null_real, float) else None
            ),
        },
        "prediction_check": {
            "predicted_model_betti": [1, 0],
            "observed_model_betti": res_model["betti_numbers"][:2],
            "matches_prediction": bool(res_model["betti_numbers"][:2] == [1, 0]),
        },
        "caveat": (
            "Per round2/ledger.json (tda.ordering_survives_variance_matching: "
            "false) and round2/reverse/tda_reduced.py's descriptive-only "
            "bottleneck comparison, TDA in this loop is a SANITY GATE, not "
            "discrimination evidence: model-cloud vs real-SN-cloud vs "
            "null-cloud bottleneck orderings have already been shown NOT to "
            "survive variance-matching once, so no reduction proposal in "
            "this round cites TDA as evidence for or against any parameter."
        ),
        "known_bottleneck_degeneracy_already_diagnosed_round2": (
            "round2/ledger.json results.round2.round1_flag_unverified_"
            "resolution: gudhi.bottleneck_distance returns half a diagram's "
            "own dominant unmatched bar whenever it exceeds 2x the other "
            "diagram's max persistence in that dimension, independent of "
            "the other diagram's content. Cited, not re-investigated from "
            "scratch; re-checked per-distance above via "
            "half_max_persistence_degeneracy_check_*."
        ),
    }
    with open(os.path.join(HERE, "tda_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    run()
