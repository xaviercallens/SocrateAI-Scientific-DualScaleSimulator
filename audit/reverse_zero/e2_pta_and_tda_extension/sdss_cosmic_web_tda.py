#!/usr/bin/env python3
"""
TDA extension (user request: "extend the data for the TDA to find
topology"). Tier X, exploratory, descriptive, UN-PREREGISTERED -- no
threshold for any Betti number or persistence statistic exists in
PRE_REGISTRATION.md, so nothing here is scored as pass/fail, and nothing
here is compared to the LeanMaster 8+16 Golay/Kummer octad structure
(that comparison is Tier C per the FRAMING RULE and per LeanMaster's own
E3/E4 verdicts, which explicitly say no cosmological prediction is
registered for it -- see improvement proposals in TDA_EXTENSION_RESULT.md).

Real data: sdss_dr17_cosmic_web_galaxies (193536 rows, ra/dec/z,
0.02<=z<=0.12, ra in [140,220], dec in [0,50]; sha256-verified per
manifest).
Matched null: sdss_dr17_random_shuffled_z (same 193536 rows, SAME
observed ra/dec, z permuted -- seed 20260919, per
data/real2/cosmic_web/make_random_catalogue.py). This is a STRICTLY
BETTER null than the Poisson/box nulls used in rounds 1-3: it preserves
the real angular selection function exactly and destroys ONLY the
radial (line-of-sight) structure, so any topological difference between
real and null isolates genuine 3D clustering, not a survey-geometry
artifact.

Method:
  - Comoving Cartesian conversion: flat LCDM, Om=0.31115 (frozen LeanMaster
    value), via comoving_distance_flat_lcdm imported UNCHANGED from
    round2/tda_gudhi.py (reuse, not reimplementation).
  - 193536 points is not Rips-able (H1/H2 blow up). Per advisor review:
    subsample (seeded) to N_SUB points and use gudhi.AlphaComplex (exact
    Delaunay-based persistence for 3D clouds, no max_edge_length fudge),
    not a sparse Rips complex.
  - Repeat over N_SEEDS independent subsample draws (real) and N_SEEDS
    independent subsample draws of the null, report the DISTRIBUTION of
    Betti numbers and total persistence per dimension, not a single draw.
  - Same pipeline, same N_SUB, same dimensions for real and null.
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import gudhi
from scipy import stats as scipy_stats

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse"
sys.path.insert(0, os.path.join(REPO_ROOT, "audit", "zero_param_loop", "round2"))
from tda_gudhi import comoving_distance_flat_lcdm  # noqa: E402  (reused unmodified)

REAL_PATH = os.path.join(REPO_ROOT, "data", "real2", "cosmic_web",
                          "sdss_dr17_galaxies_ra140_220_dec0_50.csv")
NULL_PATH = os.path.join(REPO_ROOT, "data", "real2", "cosmic_web",
                          "sdss_dr17_random_shuffled_z_ra140_220_dec0_50.csv")

N_SUB = 1200          # points per subsample draw (alpha-complex tractable in 3D)
N_SEEDS = 8           # independent subsample draws per catalogue
SUBSAMPLE_SEEDS = list(range(20260919, 20260919 + N_SEEDS))  # stated, deterministic
MAX_ALPHA_DIM = 2


def load_cloud(path):
    df = pd.read_csv(path, comment="#")
    df = df[np.isfinite(df["ra"]) & np.isfinite(df["dec"]) & np.isfinite(df["z"])].copy()
    ra = np.radians(df["ra"].to_numpy())
    dec = np.radians(df["dec"].to_numpy())
    z = df["z"].to_numpy()
    r = comoving_distance_flat_lcdm(z)  # dimensionless, c/H0 units
    x = r * np.cos(dec) * np.cos(ra)
    y = r * np.cos(dec) * np.sin(ra)
    zc = r * np.sin(dec)
    return np.column_stack([x, y, zc]), len(df)


def rescale_unit_diam(X):
    from scipy.spatial.distance import pdist
    d = pdist(X)
    maxd = float(np.max(d)) if len(d) else 1.0
    from scipy.spatial import cKDTree
    tree = cKDTree(X)
    nn_d, _ = tree.query(X, k=2)
    mean_nn = float(np.mean(nn_d[:, 1]))
    if maxd <= 0:
        maxd = 1.0
    return X / maxd, maxd, mean_nn


def alpha_persistence(X_sub):
    Xr, raw_diameter, raw_mean_nn = rescale_unit_diam(X_sub)
    ac = gudhi.AlphaComplex(points=Xr.tolist())
    st = ac.create_simplex_tree()
    st.persistence(homology_coeff_field=2, min_persistence=0.0)
    betti = st.betti_numbers()
    while len(betti) < (MAX_ALPHA_DIM + 1):
        betti.append(0)
    by_dim_total_persistence = {}
    for dim in range(MAX_ALPHA_DIM + 1):
        pairs = st.persistence_intervals_in_dimension(dim)
        finite = [(b, d) for (b, d) in pairs if np.isfinite(d)]
        total = float(sum(d - b for b, d in finite))
        by_dim_total_persistence[str(dim)] = {
            "n_finite_bars": len(finite),
            "total_persistence": total,
            "max_persistence": float(max((d - b for b, d in finite), default=0.0)),
        }
    return {"betti_numbers": betti[:MAX_ALPHA_DIM + 1],
            "by_dim": by_dim_total_persistence,
            "raw_diameter_pre_rescale": raw_diameter,
            "raw_mean_nn_distance_pre_rescale": raw_mean_nn,
            "mean_nn_over_diameter": raw_mean_nn / raw_diameter if raw_diameter else None}


def run_catalogue(label, path):
    X_full, n_before = load_cloud(path)
    results = []
    for seed in SUBSAMPLE_SEEDS:
        rng = np.random.default_rng(seed)
        idx = rng.choice(X_full.shape[0], size=N_SUB, replace=False)
        X_sub = X_full[idx]
        res = alpha_persistence(X_sub)
        res["subsample_seed"] = seed
        results.append(res)
    return {
        "label": label,
        "path": path,
        "n_rows_loaded": n_before,
        "n_sub": N_SUB,
        "n_seeds": N_SEEDS,
        "draws": results,
    }


def summarize(cat_result):
    b0 = [d["betti_numbers"][0] for d in cat_result["draws"]]
    b1 = [d["betti_numbers"][1] for d in cat_result["draws"]]
    b2 = [d["betti_numbers"][2] for d in cat_result["draws"]]
    tp1 = [d["by_dim"]["1"]["total_persistence"] for d in cat_result["draws"]]
    tp2 = [d["by_dim"]["2"]["total_persistence"] for d in cat_result["draws"]]
    diam = [d["raw_diameter_pre_rescale"] for d in cat_result["draws"]]
    mean_nn = [d["raw_mean_nn_distance_pre_rescale"] for d in cat_result["draws"]]
    nn_over_diam = [d["mean_nn_over_diameter"] for d in cat_result["draws"]]
    def stats(v):
        return {"mean": float(np.mean(v)), "std": float(np.std(v)), "values": v}
    return {
        "betti0": stats(b0), "betti1": stats(b1), "betti2": stats(b2),
        "total_persistence_H1": stats(tp1), "total_persistence_H2": stats(tp2),
        "raw_diameter_pre_rescale": stats(diam),
        "raw_mean_nn_distance_pre_rescale": stats(mean_nn),
        "mean_nn_over_diameter": stats(nn_over_diam),
    }


def main():
    real = run_catalogue("sdss_dr17_real", REAL_PATH)
    null = run_catalogue("sdss_dr17_shuffled_z_null", NULL_PATH)
    real_summary = summarize(real)
    null_summary = summarize(null)
    welch_tests = {}
    for key in ("total_persistence_H1", "total_persistence_H2",
                "raw_diameter_pre_rescale", "raw_mean_nn_distance_pre_rescale",
                "mean_nn_over_diameter"):
        a = real_summary[key]["values"]
        b = null_summary[key]["values"]
        t, p = scipy_stats.ttest_ind(a, b, equal_var=False)
        welch_tests[key] = {
            "t_statistic": float(t), "p_value": float(p),
            "real_mean": float(np.mean(a)), "null_mean": float(np.mean(b)),
            "note": "Welch's t-test, unequal variance, n=8 seeds per side. "
                "Exploratory (Tier X) -- NOT a pre-registered test; n=8 "
                "subsample draws from overlapping parent catalogues "
                "(draws are not fully independent of each other); direction "
                "and significance should be treated as a signal to "
                "pre-register and re-test, not a confirmed result.",
        }
    summary = {
        "real": real_summary,
        "null_shuffled_z": null_summary,
        "welch_t_tests_real_vs_null": welch_tests,
        "method": {
            "comoving_distance": "flat LCDM, Om=0.31115, imported unmodified "
                "from audit/zero_param_loop/round2/tda_gudhi.comoving_distance_flat_lcdm",
            "complex": "gudhi.AlphaComplex (exact Delaunay, 3D) -- NOT sparse Rips "
                "(193536 points is not Rips-tractable at H1/H2)",
            "n_sub": N_SUB, "n_seeds": N_SEEDS, "subsample_seeds": SUBSAMPLE_SEEDS,
            "null_construction": "sdss_dr17_random_shuffled_z: SAME observed "
                "ra/dec footprint as real, z permuted (seed 20260919, "
                "data/real2/cosmic_web/make_random_catalogue.py) -- preserves "
                "angular selection function exactly, destroys only radial "
                "(line-of-sight) structure. Strictly better than the Poisson/box "
                "nulls used in rounds 1-3 for isolating genuine 3D clustering "
                "from survey-geometry artifacts.",
        },
        "tier": "X (exploratory, descriptive, un-preregistered)",
        "framing_note": "No threshold for any Betti number or persistence "
            "statistic exists in PRE_REGISTRATION.md. This result is NOT "
            "compared to the LeanMaster Golay-octad/Kummer 8+16 split "
            "(Tier C, no registered cosmological prediction per Streams "
            "6-8 E3/E4 verdicts). A real-vs-null difference here would say "
            "only 'the survey shows more/less low-dimensional topological "
            "structure than a radially-scrambled version of itself' -- a "
            "descriptive fact about SDSS DR17 cosmic-web structure, not a "
            "test of K3xT2 or of mu_sym/c4_pta_product (parameter_effect: none).",
    }
    out = {"real_draws": real, "null_draws": null, "summary": summary}
    out_path = os.path.join(HERE, "sdss_cosmic_web_tda_result.json")
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2)
    print(json.dumps(summary, indent=2))
    print(f"\nwrote {out_path}")


if __name__ == "__main__":
    main()
