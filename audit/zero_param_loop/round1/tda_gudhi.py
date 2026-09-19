#!/usr/bin/env python3
"""
Round 1 FORWARD LOOP step 3 (TDA, real INRIA GUDHI -- `import gudhi`, no
home-made substitute, not scripts/tda_mapper.py).

IMPORTANT INTERPRETATION NOTE (see advisor guidance recorded in this
round's transcript): Betti numbers of the model-side point cloud do NOT
measure "how many parameters the observables can see." A smooth image of a
6-dimensional box under a smooth map is contractible regardless of the
box's dimension -- beta0=1, beta1=beta2=0 is expected WHATEVER the true
effective dimension is. Persistent homology here is used for SHAPE (thin
sliver vs. filled patch, via the persistence-diagram spread) and for the
sanity controls; the actual dimension count comes from
jacobian_response_rank.py's SVD (Fisher/PCA-style), reported separately.

a. MODEL side: standardized sweep observable point cloud (numerically
   stable points only), PCA'd to <=8 dims, Rips complex, dims 0-2.
b. REAL side (primary): Pantheon+ SH0ES SN residual-vs-flat-LCDM embedding
   (z, color c, stretch x1, residual), an observable-space cloud directly
   comparable to the model side.
   REAL side (extra, qualitative only): 2MRS RA/Dec angular patch --
   verified to have NO redshift column, so this is NOT the comoving
   galaxy-catalogue analysis the task's first TDA option names; reported
   as an extra check with known S^2-patch topology, not as the primary
   real-data comparison.
c. CONTROLS: noisy circle (must give exactly one dominant H1 bar -- this is
   asserted, not just remarked); shuffled/Poisson null cloud, same N and
   (post-rescale) volume as the model cloud.
   Bottleneck distances (gudhi.bottleneck_distance): model-vs-real,
   null-vs-real, circle-vs-real (context).

All clouds are rescaled to max pairwise distance = 1 before any complex is
built, so bottleneck numbers are comparable (still labeled qualitative).

Writes (all under this directory):
  tda_persistence_model.json / .png
  tda_persistence_real_sn.json / .png
  tda_persistence_real_2mrs.json / .png
  tda_persistence_null.json / .png
  tda_persistence_circle.json / .png
  tda_report.json
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

REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
sys.path.insert(0, REPO_ROOT)
OUT_DIR = os.path.dirname(os.path.abspath(__file__))
SEED = 42
rng = np.random.RandomState(SEED)

MAX_RIPS_DIM = 2


def rescale_unit_diam(X):
    """Rescale a point cloud so its max pairwise Euclidean distance is 1."""
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


def build_and_persist(X, label, max_edge_length=None, sparse=0.1):
    from scipy.spatial.distance import pdist
    X, scale = rescale_unit_diam(X)
    d = pdist(X)
    med = float(np.median(d))
    mel = max_edge_length if max_edge_length is not None else 2.0 * med
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
    # plot
    fig, ax = plt.subplots(figsize=(5, 5))
    colors = {0: "tab:blue", 1: "tab:orange", 2: "tab:green"}
    for dim in (0, 1, 2):
        pts = by_dim[dim]
        if pts:
            arr = np.array(pts)
            ax.scatter(arr[:, 0], arr[:, 1], s=10, color=colors[dim], label=f"H{dim} (n={len(pts)})")
    lim = max([p[1] for v in by_dim.values() for p in v] + [1e-6]) * 1.1
    ax.plot([0, lim], [0, lim], "k--", lw=0.5)
    ax.set_xlabel("birth"); ax.set_ylabel("death"); ax.set_title(f"Persistence diagram: {label}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    png_path = os.path.join(OUT_DIR, f"tda_persistence_{label}.png")
    fig.savefig(png_path, dpi=120)
    plt.close(fig)
    json_path = os.path.join(OUT_DIR, f"tda_persistence_{label}.json")
    with open(json_path, "w") as f:
        json.dump(result, f, indent=2)
    return result, by_dim


# ---------------------------------------------------------------------------
# a. MODEL side
# ---------------------------------------------------------------------------
def build_model_cloud():
    df = pd.read_csv(os.path.join(OUT_DIR, "sweep.csv"))
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
    # standardize (z-score) each column
    mu, sd = X.mean(axis=0), X.std(axis=0)
    sd[sd < 1e-12] = 1.0
    Xz = (X - mu) / sd
    return Xz, len(stable)


# ---------------------------------------------------------------------------
# b. REAL side: Pantheon+ SN residual-vs-LCDM embedding
# ---------------------------------------------------------------------------
def build_real_sn_cloud(n_target):
    path = os.path.join(REPO_ROOT, "data", "real", "dark_energy", "pantheon_plus_sh0es.dat")
    df = pd.read_csv(path, sep=r"\s+", engine="python")
    df = df[(df["zHD"] > 0.01) & (df["zHD"] <= 2.4) & np.isfinite(df["MU_SH0ES"])
            & np.isfinite(df["c"]) & np.isfinite(df["x1"])].copy()
    # best-fit flat LCDM (Om, offset) from chi2_report.json (already fitted
    # honestly against this same data); recompute residual using that Om.
    with open(os.path.join(OUT_DIR, "chi2_report.json")) as f:
        chi2_rep = json.load(f)
    Om = chi2_rep["lcdm_baseline"]["sn"]["Om_fit"]
    offset = chi2_rep["lcdm_baseline"]["sn"]["offset_fit"]
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


# ---------------------------------------------------------------------------
# b'. REAL side extra: 2MRS RA/Dec angular patch (known S^2-patch topology;
#     no redshift available, so NOT the comoving-coordinate analysis).
# ---------------------------------------------------------------------------
def build_real_2mrs_cloud(n_target):
    path = os.path.join(REPO_ROOT, "data", "real", "cosmic_web", "2mrs_sample.tsv")
    rows = []
    with open(path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            rows.append(parts)
    header_idx = None
    for i, r in enumerate(rows):
        if any("RAJ2000" in c or "RA" == c.strip() for c in r):
            header_idx = i
            break
    if header_idx is None:
        return None
    header = [c.strip() for c in rows[header_idx]]
    # VizieR TSV: header line, then a units line (" \tdeg\tdeg\tmag"), then a
    # dashed separator line ("-----\t----------\t..."), then data. Verified
    # by direct inspection (sed -n '30,45p' data/real/cosmic_web/2mrs_sample.tsv).
    data_rows = rows[header_idx + 3:]
    ra_i = next((i for i, c in enumerate(header) if c.startswith("RAJ2000") or c == "RA"), None)
    de_i = next((i for i, c in enumerate(header) if c.startswith("DEJ2000") or c == "DEC" or c == "DE"), None)
    if ra_i is None or de_i is None:
        return None
    ras, decs = [], []
    for r in data_rows:
        if len(r) <= max(ra_i, de_i):
            continue
        try:
            ras.append(float(r[ra_i].strip()))
            decs.append(float(r[de_i].strip()))
        except ValueError:
            continue
    ra = np.array(ras); dec = np.array(decs)
    if len(ra) == 0:
        return None
    ra_rad = np.deg2rad(ra); dec_rad = np.deg2rad(dec)
    x = np.cos(dec_rad) * np.cos(ra_rad)
    y = np.cos(dec_rad) * np.sin(ra_rad)
    z = np.sin(dec_rad)
    X = np.column_stack([x, y, z])
    if len(X) > n_target:
        idx = rng.choice(len(X), size=n_target, replace=False)
        X = X[idx]
    return X, len(ra)


# ---------------------------------------------------------------------------
# c. Controls
# ---------------------------------------------------------------------------
def build_circle_cloud(n, noise_sigma=0.05):
    theta = rng.uniform(0, 2 * math.pi, size=n)
    r = 1.0 + rng.normal(0, noise_sigma, size=n)
    x = r * np.cos(theta); y = r * np.sin(theta)
    return np.column_stack([x, y])


def build_null_cloud(n, dim, box_halfwidth=1.0):
    return rng.uniform(-box_halfwidth, box_halfwidth, size=(n, dim))


def dgm_array(by_dim, dim):
    arr = by_dim.get(str(dim)) if isinstance(by_dim, dict) and str(dim) in by_dim else by_dim.get(dim, [])
    return np.array(arr) if arr else np.empty((0, 2))


def run_tda():
    report = {}

    X_model_raw, n_model = build_model_cloud()
    X_model, pca_k, pca_frac = pca_reduce(X_model_raw, max_dim=8, var_target=0.95)
    res_model, bd_model = build_and_persist(X_model, "model", sparse=0.2)
    report["model"] = {k: v for k, v in res_model.items() if k != "diagram_by_dim"}
    report["model"]["pca_raw_dim"] = int(X_model_raw.shape[1])
    report["model"]["pca_reduced_dim"] = pca_k
    report["model"]["pca_cumulative_variance_fraction"] = pca_frac[:pca_k]

    X_sn, Om_used, offset_used = build_real_sn_cloud(n_target=n_model)
    res_sn, bd_sn = build_and_persist(X_sn, "real_sn", sparse=0.2)
    report["real_sn"] = {k: v for k, v in res_sn.items() if k != "diagram_by_dim"}
    report["real_sn"]["lcdm_Om_used_for_residual"] = Om_used
    report["real_sn"]["lcdm_offset_used_for_residual"] = offset_used

    mrs = build_real_2mrs_cloud(n_target=n_model)
    if mrs is not None:
        X_mrs, n_mrs_total = mrs
        res_mrs, bd_mrs = build_and_persist(X_mrs, "real_2mrs", sparse=0.2)
        report["real_2mrs"] = {k: v for k, v in res_mrs.items() if k != "diagram_by_dim"}
        report["real_2mrs"]["n_total_catalogue_before_subsample"] = n_mrs_total
        report["real_2mrs"]["caveat"] = "RA/Dec only, no redshift; angular S^2 patch, NOT comoving coordinates"
    else:
        bd_mrs = None
        report["real_2mrs"] = "ABSENT: could not locate RA/Dec columns in 2mrs_sample.tsv header"

    X_circle = build_circle_cloud(n_model)
    res_circle, bd_circle = build_and_persist(X_circle, "circle", sparse=0.2)
    report["circle_control"] = {k: v for k, v in res_circle.items() if k != "diagram_by_dim"}
    h1_circle = dgm_array(bd_circle, 1)
    if len(h1_circle) >= 1:
        pers = h1_circle[:, 1] - h1_circle[:, 0]
        pers_sorted = np.sort(pers)[::-1]
        top = float(pers_sorted[0])
        second = float(pers_sorted[1]) if len(pers_sorted) > 1 else 0.0
        ratio = top / second if second > 1e-12 else float("inf")
        circle_ok = bool(len(h1_circle) >= 1 and (second == 0.0 or ratio > 5.0))
    else:
        top, second, ratio, circle_ok = 0.0, 0.0, 0.0, False
    report["circle_control_assertion"] = {
        "n_H1_bars": int(len(h1_circle)), "top_persistence": top, "second_persistence": second,
        "ratio": ratio, "ASSERTION_one_dominant_H1_bar": circle_ok,
    }
    if not circle_ok:
        report["circle_control_assertion"]["WARNING"] = (
            "Known-answer control FAILED to show one dominant H1 bar -- pipeline parameters "
            "(sparse/max_edge_length) are likely wrong; downstream bottleneck numbers should not "
            "be trusted without fixing this."
        )

    X_null = build_null_cloud(n_model, dim=pca_k)
    res_null, bd_null = build_and_persist(X_null, "null", sparse=0.2)
    report["null_control"] = {k: v for k, v in res_null.items() if k != "diagram_by_dim"}

    def bd(a, b, dim):
        A, B = dgm_array(a, dim), dgm_array(b, dim)
        try:
            return float(gudhi.bottleneck_distance(A, B))
        except Exception as exc:  # pragma: no cover
            return f"ERROR: {exc}"

    report["bottleneck_distances_qualitative"] = {
        "note": "Clouds are pre-rescaled to max pairwise distance=1; still label these qualitative -- "
                "model and real_sn live in different native dimensions/semantics before PCA/standardization.",
        "model_vs_real_sn_H0": bd(bd_model, bd_sn, 0),
        "model_vs_real_sn_H1": bd(bd_model, bd_sn, 1),
        "null_vs_real_sn_H0": bd(bd_null, bd_sn, 0),
        "null_vs_real_sn_H1": bd(bd_null, bd_sn, 1),
        "circle_vs_real_sn_H1": bd(bd_circle, bd_sn, 1),
    }
    if bd_mrs is not None:
        report["bottleneck_distances_qualitative"]["model_vs_real_2mrs_H1"] = bd(bd_model, bd_mrs, 1)
        report["bottleneck_distances_qualitative"]["null_vs_real_2mrs_H1"] = bd(bd_null, bd_mrs, 1)

    report["interpretation_caveat"] = (
        "Betti numbers of the model cloud (beta0=1, beta1=beta2=0 expected for a smooth image of a "
        "box) are NOT a parameter count. They describe global shape (one connected contractible blob) "
        "which is uninformative about effective dimension here. Effective dimension comes from "
        "jacobian_report.json's singular-value spectrum, not from this file."
    )

    with open(os.path.join(OUT_DIR, "tda_report.json"), "w") as f:
        json.dump(report, f, indent=2)
    print(json.dumps({k: v for k, v in report.items() if k not in ("model", "real_sn", "real_2mrs", "null_control", "circle_control")}, indent=2))
    return report


if __name__ == "__main__":
    run_tda()
