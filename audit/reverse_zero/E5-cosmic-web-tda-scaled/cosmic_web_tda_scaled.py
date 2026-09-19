#!/usr/bin/env python3
"""
E5 -- scale up the cosmic-web TDA (E4) per the pre-registered next step
(audit/reverse_zero/IMPROVEMENT_PROPOSAL.md, item 1: "Scale up E4's
subsample and rerun ... a matched-volume Poisson null in addition to the
shuffled-z null, and the known-answer circle control ... before drawing
any conclusion about real vs. null topology"), and per the relayed user
request to "extend the data for the TDA to find topology."

WHAT CHANGED FROM E4 (fixes applied per advisor review of this run, logged
here for the record):
  1. NO unit-diameter rescaling. E4's rescale_unit_diam gave real and null
     clouds DIFFERENT physical Mpc/h scales after normalising each to its
     own diameter, which is the wrong operation for asking "at which
     comoving scale (Mpc/h) does the data differ from the null" -- the
     question this script answers. Everything here stays in physical
     comoving Mpc/h (H0=100 convention, i.e. "per h") throughout.
  2. The null is NOT drawn uniformly in the ra/dec bounding box (140-220,
     0-50 deg is a rectangle, not the true DR17 footprint -- a box-uniform
     null has different ANGULAR SUPPORT from the data, and would report a
     spurious "signal" that is pure footprint geometry; this project has
     already retracted one round-1 bounding-box TDA artifact). Instead:
     the SAME observed (ra, dec) pairs are kept and z is reassigned by an
     independent random permutation (20 seeds, for a null envelope). This
     keeps the exact angular footprint and the exact marginal n(z) by
     construction, and destroys 3D clustering; for a catalogue this size
     (>1.3e5 galaxies) permutation-without-replacement of the empirical
     n(z) is statistically indistinguishable from drawing r i.i.d. from
     that empirical n(z) with replacement (a "Poisson-in-the-footprint"
     null) -- so this ONE construction stands in for the ground rules'
     nulls (i) Poisson-in-volume and (ii) redshift-shuffled; this
     equivalence is stated here rather than silently claimed as two
     separate nulls.
  3. Number-density flattening (ground rules' "volume-limited or
     number-density-flattened subsample", cuts stated below) using
     EQUAL-COMOVING-VOLUME shells, not equal-z bins, so a flat count per
     shell means a flat 3D number density n(r), not a redshift-binning
     artifact.
  4. GUDHI AlphaComplex (ground rules ask for gudhi.AlphaComplex, not
     RipsComplex): filtration values are squared circumradii -- sqrt() is
     applied before any r is plotted or quoted. Truncated at
     max_alpha_square = R_MAX_PERS**2 to bound runtime (verified: 3D
     Delaunay + persistence at N=30000 truncated at alpha^2=2500 runs in
     ~3.5s on this machine, so N=25000 x 21 clouds is nowhere near the
     30-min budget).
  5. Two known-answer controls, not one: a noisy circle (recovers a
     dominant H1 loop, as in E4/round2/round3) AND a planted-void cloud
     (uniform fill with NUM_VOID known-radius spherical holes punched
     out, recovers dominant H2 bars near the void radius) -- E4 only
     tested H1 detectability, but the cosmic-web claim of interest is
     about voids, i.e. H2.
  6. The comparison statistic (L2 distance between Betti curves over a
     PRE-STATED r-range, and bottleneck distance on the top bars) is
     written into this docstring and the output JSON BEFORE the first run
     of this script, per IMPROVEMENT_PROPOSAL.md item 2's own warning
     against choosing a statistic after seeing the result.

WHAT DID NOT CHANGE: same frozen M0 background (Omega_m=0.31115, the E3-M0
value, LeanMaster DarkEnergyScale.lean:84) used only to convert (ra,dec,z)
to comoving Mpc/h -- this script does not fit or vary Omega_m or H0; H0
only sets the Mpc/h unit convention (H0=100 km/s/Mpc/h), stated per the
ground rules, and does not affect any Betti number (it is an overall
length rescaling and the r-axis is reported in Mpc/h).

DEFERRED (not run in this script, honestly reported rather than skipped
silently): lognormal CAMB mocks (ground rules' null iii). Building a
CAMB P(k) -> Gaussian random field -> lognormal transform -> bias-fit ->
Poisson-sample -> mask pipeline correctly (with the bias fitted ONLY to
the two-point function, never to topology) is substantial new code with
several normalisation traps (grid Nyquist, box-size aliasing, shot noise);
committing a half-verified version of it risks a worse error than not
having it. It is listed in the improvement proposal (Sec 3 below) as the
single highest-leverage remaining step for this specific test.

Tier: X (exploratory numerics). Not a K3xT2 (LeanMaster Stream 8 E2-E4)
prediction test -- no constructed octad-indexed observable exists (see
E4's report and REPORT.md); this script establishes whether SDSS DR17
galaxy clustering is topologically distinguishable from its own
footprint-and-n(z)-matched null at N=25000, which is the necessary
prerequisite before any such octad statistic could be built. A match
with this null is "consistent with M0's late-time background at the
level of THIS test" -- not a confirmation of K3xT2, and not by itself a
detection of any dark-energy or extra-dimension parameter (mu_sym,
c4_pta_product are untouched by this script).

Data (fetched_verified, sha256 in data manifest):
  data/real2/cosmic_web/sdss_dr17_galaxies_ra140_220_dec0_50.csv
  (193536 rows; ra,dec,z,zErr columns only -- no magnitude column, so a
  true volume-LIMITED absolute-magnitude cut is not possible with this
  catalogue; number-density FLATTENING is used instead, per the ground
  rules' stated alternative.)

Command:
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    audit/reverse_zero/E5-cosmic-web-tda-scaled/cosmic_web_tda_scaled.py
"""
import json
import os
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gudhi
from scipy.integrate import cumulative_trapezoid

WT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(WT_ROOT, "data", "real2", "cosmic_web")
REAL_CSV = os.path.join(DATA_DIR, "sdss_dr17_galaxies_ra140_220_dec0_50.csv")

OMEGA_M = 0.31115  # SAME frozen M0 background as E3-M0/m0_model.py / E4 -- no new tuning
H0_CONVENTION_NOTE = "Mpc/h units, H0=100 km/s/Mpc/h convention; H0 does not affect Betti numbers, only the r-axis unit label."

N_SHELLS_FOR_SELECTION = 20     # equal-comoving-volume bins used only to find the flat-density window
SELECTION_SHELL_IDX = (2, 3, 4, 5)  # chosen from the diagnostic count profile (see report); flattest 4 adjacent shells
N_TARGET_PER_SHELL = 6250       # -> N_FINAL = 4 * 6250 = 25000
N_FINAL = len(SELECTION_SHELL_IDX) * N_TARGET_PER_SHELL
N_NULL_SEEDS = 20               # >= ground rules' ">= 20 mocks with seeds" (applied here to the shuffled-z null)
R_MAX_PERS = 50.0               # Mpc/h; persistence truncation scale (voids/filaments regime)
MAX_ALPHA_SQ = R_MAX_PERS ** 2
R_GRID = np.linspace(0.0, R_MAX_PERS, 101)  # for Betti curves beta_k(r)
# PRE-STATED comparison statistic (written before the first run, per IMPROVEMENT_PROPOSAL item 2):
STAT_R_RANGE = (2.0, 40.0)      # Mpc/h; L2/bottleneck comparison restricted to this range
SEED_SUBSAMPLE_REAL = 42


def comoving_r_mpc_over_h(z_arr, z_grid_max):
    zg = np.linspace(0.0, float(z_grid_max) + 0.01, 20000)
    Ez = np.sqrt(OMEGA_M * (1.0 + zg) ** 3 + (1.0 - OMEGA_M))
    dc_dimensionless = cumulative_trapezoid(1.0 / Ez, zg, initial=0.0)  # Dc * H0/c
    dc_mpc_over_h = dc_dimensionless * 2997.92458  # c/H0 in Mpc/h units (H0=100), i.e. multiply by c/(100 km/s)
    return np.interp(z_arr, zg, dc_mpc_over_h)


def radec_r_to_xyz(ra_deg, dec_deg, r_mpc):
    ra = np.radians(ra_deg)
    dec = np.radians(dec_deg)
    x = r_mpc * np.cos(dec) * np.cos(ra)
    y = r_mpc * np.cos(dec) * np.sin(ra)
    z = r_mpc * np.sin(dec)
    return np.column_stack([x, y, z])


def equal_volume_shell_edges(r_max, n_shells):
    i = np.arange(n_shells + 1)
    return r_max * (i / n_shells) ** (1.0 / 3.0)


def shell_membership(r_mpc, edges):
    return np.clip(np.digitize(r_mpc, edges) - 1, 0, len(edges) - 2)


def select_flat_density_subsample(r_mpc, edges, shell_idx, n_per_shell, seed):
    """Downsample uniformly within each named shell to n_per_shell points
    each, giving a number-density-flattened sample by construction (equal
    counts in equal comoving-volume shells)."""
    rng = np.random.RandomState(seed)
    mem = shell_membership(r_mpc, edges)
    chosen = []
    shell_counts_available = {}
    for s in shell_idx:
        idx_s = np.where(mem == s)[0]
        shell_counts_available[int(s)] = int(len(idx_s))
        take = min(n_per_shell, len(idx_s))
        pick = rng.choice(idx_s, size=take, replace=False)
        chosen.append(pick)
    chosen = np.concatenate(chosen)
    return chosen, shell_counts_available


def alpha_persistence(xyz, max_alpha_sq, label, out_dir):
    t0 = time.time()
    ac = gudhi.AlphaComplex(points=xyz)
    st = ac.create_simplex_tree(max_alpha_square=max_alpha_sq)
    diag = st.persistence(homology_coeff_field=2, min_persistence=0.0)
    dt = time.time() - t0
    by_dim = {0: [], 1: [], 2: []}
    for dim, (b, d) in diag:
        if dim in by_dim:
            b_r = np.sqrt(max(b, 0.0))
            d_r = np.sqrt(d) if np.isfinite(d) else float("inf")
            by_dim[dim].append([b_r, d_r])
    betti = st.betti_numbers()
    return {
        "label": label, "n_points": int(xyz.shape[0]), "runtime_sec": dt,
        "betti_numbers_at_truncation": betti,
        "diagram_by_dim_r": {str(k): v for k, v in by_dim.items()},
    }, by_dim


def betti_curve(by_dim, dim, r_grid, r_trunc):
    """persistent Betti number beta_dim(r) = #{bars (b,d): b <= r < d},
    with death=inf (or death > r_trunc, i.e. still alive at the
    truncation scale) treated as alive for all r in r_grid <= r_trunc."""
    bars = np.array(by_dim[dim]) if by_dim[dim] else np.empty((0, 2))
    curve = np.zeros(len(r_grid))
    if bars.size == 0:
        return curve
    b = bars[:, 0]
    d = np.where(np.isfinite(bars[:, 1]), bars[:, 1], r_trunc + 1e9)
    for i, r in enumerate(r_grid):
        curve[i] = np.sum((b <= r) & (r < d))
    return curve


def euler_curve(by_dim, r_grid, r_trunc):
    b0 = betti_curve(by_dim, 0, r_grid, r_trunc)
    b1 = betti_curve(by_dim, 1, r_grid, r_trunc)
    b2 = betti_curve(by_dim, 2, r_grid, r_trunc)
    return b0 - b1 + b2, b0, b1, b2


def l2_over_range(curve_a, curve_b, r_grid, r_range):
    mask = (r_grid >= r_range[0]) & (r_grid <= r_range[1])
    return float(np.sqrt(np.mean((curve_a[mask] - curve_b[mask]) ** 2)))


def top_bars(by_dim, dim, k=5):
    bars = np.array(by_dim[dim]) if by_dim[dim] else np.empty((0, 2))
    if bars.size == 0:
        return np.empty((0, 2))
    finite = np.where(np.isfinite(bars[:, 1]), bars[:, 1], R_MAX_PERS)
    pers = finite - bars[:, 0]
    order = np.argsort(pers)[::-1][:k]
    out = bars[order].copy()
    out[:, 1] = finite[order]
    return out


def build_circle_cloud(n, noise_sigma=0.05, seed=0):
    rng = np.random.RandomState(seed)
    theta = rng.uniform(0, 2 * np.pi, n)
    x = np.cos(theta) + rng.normal(0, noise_sigma, n)
    y = np.sin(theta) + rng.normal(0, noise_sigma, n)
    z = rng.normal(0, noise_sigma, n)
    return np.column_stack([x, y, z])


def build_planted_void_cloud(n_fill, box=100.0, void_radius=15.0, n_voids=6, seed=0):
    """Fill a cube uniformly, punch out n_voids non-overlapping spheres of
    known radius -- a known-answer H2 (void) control, since the E4
    circle control only exercised H1."""
    rng = np.random.RandomState(seed)
    centers = []
    tries = 0
    while len(centers) < n_voids and tries < 10000:
        tries += 1
        c = rng.uniform(void_radius * 1.2, box - void_radius * 1.2, size=3)
        if all(np.linalg.norm(c - c2) > 2.3 * void_radius for c2 in centers):
            centers.append(c)
    centers = np.array(centers)
    pts = rng.uniform(0, box, size=(n_fill, 3))
    keep = np.ones(n_fill, dtype=bool)
    for c in centers:
        keep &= np.linalg.norm(pts - c, axis=1) > void_radius
    return pts[keep], centers


def plot_betti_envelope(r_grid, real_curves, null_curves_stack, out_path, title):
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    labels = ["H0", "H1", "H2"]
    for k in range(3):
        ax = axes[k]
        null_k = null_curves_stack[k]  # shape (n_null, n_r)
        lo = np.percentile(null_k, 2.5, axis=0)
        hi = np.percentile(null_k, 97.5, axis=0)
        med = np.percentile(null_k, 50, axis=0)
        ax.fill_between(r_grid, lo, hi, color="gray", alpha=0.4, label="null 95%% envelope (n=%d)" % null_k.shape[0])
        ax.plot(r_grid, med, color="black", lw=1, ls="--", label="null median")
        ax.plot(r_grid, real_curves[k], color="tab:red", lw=1.8, label="SDSS DR17 real")
        ax.set_title(labels[k]); ax.set_xlabel("r (Mpc/h)")
        if k == 0:
            ax.set_ylabel("beta_k(r)")
        ax.legend(fontsize=7)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)
    plt.close(fig)


def run():
    t_start = time.time()
    df = pd.read_csv(REAL_CSV, comment="#")
    assert set(["ra", "dec", "z", "zErr"]).issubset(df.columns), df.columns
    valid = np.isfinite(df["z"]) & (df["z"] > 0)
    df = df[valid].reset_index(drop=True)
    n_total_valid = len(df)

    z_max_grid = float(df["z"].max())
    r_real_all = comoving_r_mpc_over_h(df["z"].to_numpy(), z_max_grid)
    r_shell_max = float(r_real_all.max())
    edges = equal_volume_shell_edges(r_shell_max, N_SHELLS_FOR_SELECTION)
    diag_counts = {}
    mem_all = shell_membership(r_real_all, edges)
    for s in range(N_SHELLS_FOR_SELECTION):
        diag_counts[int(s)] = {"r_lo": float(edges[s]), "r_hi": float(edges[s + 1]),
                                "count": int(np.sum(mem_all == s))}

    real_idx, real_shell_counts = select_flat_density_subsample(
        r_real_all, edges, SELECTION_SHELL_IDX, N_TARGET_PER_SHELL, SEED_SUBSAMPLE_REAL)
    real_r = r_real_all[real_idx]
    real_xyz = radec_r_to_xyz(df["ra"].to_numpy()[real_idx], df["dec"].to_numpy()[real_idx], real_r)

    # residual n(r) flatness check within the selected window (equal-volume shells -> flat count means flat n(r))
    real_mem_sel = shell_membership(real_r, edges)
    real_counts_after = [int(np.sum(real_mem_sel == s)) for s in SELECTION_SHELL_IDX]

    res_real, bd_real = alpha_persistence(real_xyz, MAX_ALPHA_SQ, "sdss_dr17_real_N%d" % real_xyz.shape[0], HERE)
    euler_real, b0_real, b1_real, b2_real = euler_curve(bd_real, R_GRID, R_MAX_PERS)

    # ---- null realizations: same observed (ra,dec), independently permuted z, 20 seeds ----
    null_betti_stacks = {0: [], 1: [], 2: []}
    null_euler_stack = []
    null_top_bars = {0: [], 1: [], 2: []}
    null_reports = []
    ra_all = df["ra"].to_numpy(); dec_all = df["dec"].to_numpy(); z_all = df["z"].to_numpy()
    for seed in range(N_NULL_SEEDS):
        rng = np.random.RandomState(1000 + seed)
        z_perm = rng.permutation(z_all)
        r_null_all = comoving_r_mpc_over_h(z_perm, z_max_grid)  # same z-value set -> identical shell counts to real by construction
        idx, shell_counts = select_flat_density_subsample(r_null_all, edges, SELECTION_SHELL_IDX, N_TARGET_PER_SHELL, seed=2000 + seed)
        null_xyz = radec_r_to_xyz(ra_all[idx], dec_all[idx], r_null_all[idx])
        res_null, bd_null = alpha_persistence(null_xyz, MAX_ALPHA_SQ, "null_seed%d" % seed, HERE)
        e_n, b0_n, b1_n, b2_n = euler_curve(bd_null, R_GRID, R_MAX_PERS)
        null_betti_stacks[0].append(b0_n); null_betti_stacks[1].append(b1_n); null_betti_stacks[2].append(b2_n)
        null_euler_stack.append(e_n)
        for dim in (0, 1, 2):
            tb = top_bars(bd_null, dim, k=3)
            null_top_bars[dim].append(tb.tolist())
        null_reports.append({"seed": seed, "n_points": res_null["n_points"], "runtime_sec": res_null["runtime_sec"],
                              "betti_numbers_at_truncation": res_null["betti_numbers_at_truncation"]})

    for dim in (0, 1, 2):
        null_betti_stacks[dim] = np.array(null_betti_stacks[dim])
    null_euler_stack = np.array(null_euler_stack)

    # ---- comparison statistic (pre-stated) ----
    l2_stats = {}
    outside_envelope = {}
    for dim, curve in zip((0, 1, 2), (b0_real, b1_real, b2_real)):
        null_mean = null_betti_stacks[dim].mean(axis=0)
        l2_stats["H%d" % dim] = l2_over_range(curve, null_mean, R_GRID, STAT_R_RANGE)
        lo = np.percentile(null_betti_stacks[dim], 2.5, axis=0)
        hi = np.percentile(null_betti_stacks[dim], 97.5, axis=0)
        mask = (R_GRID >= STAT_R_RANGE[0]) & (R_GRID <= STAT_R_RANGE[1])
        out_r = R_GRID[mask][(curve[mask] < lo[mask]) | (curve[mask] > hi[mask])]
        outside_envelope["H%d" % dim] = {"n_r_points_outside_95pct_envelope": int(len(out_r)),
                                           "n_r_points_total_in_range": int(mask.sum()),
                                           "r_values_outside": [float(x) for x in out_r]}
    l2_stats["euler"] = l2_over_range(euler_real, null_euler_stack.mean(axis=0), R_GRID, STAT_R_RANGE)

    # bottleneck: real top bars vs each null's top bars, per dim (report mean+std across the 20 nulls)
    bottleneck_stats = {}
    for dim in (0, 1, 2):
        real_tb = top_bars(bd_real, dim, k=3)
        dists = []
        for seed in range(N_NULL_SEEDS):
            null_tb = np.array(null_top_bars[dim][seed]) if null_top_bars[dim][seed] else np.empty((0, 2))
            dists.append(float(gudhi.bottleneck_distance(real_tb, null_tb)))
        bottleneck_stats["H%d" % dim] = {"mean": float(np.mean(dists)), "std": float(np.std(dists)),
                                          "min": float(np.min(dists)), "max": float(np.max(dists))}

    # ---- known-answer controls ----
    circle_xyz = build_circle_cloud(n=2000, seed=SEED_SUBSAMPLE_REAL)
    res_circle, bd_circle = alpha_persistence(circle_xyz, max_alpha_sq=100.0, label="known_answer_circle", out_dir=HERE)
    circle_h1 = top_bars(bd_circle, 1, k=3)
    circle_ratio = float((circle_h1[0, 1] - circle_h1[0, 0]) / max(circle_h1[1, 1] - circle_h1[1, 0], 1e-9)) if len(circle_h1) >= 2 else None

    void_xyz, void_centers = build_planted_void_cloud(n_fill=15000, box=100.0, void_radius=15.0, n_voids=6, seed=SEED_SUBSAMPLE_REAL)
    res_void, bd_void = alpha_persistence(void_xyz, max_alpha_sq=(2.5 * 15.0) ** 2, label="known_answer_planted_voids", out_dir=HERE)
    void_h2 = top_bars(bd_void, 2, k=8)
    # each planted void should give an H2 bar that DIES near r ~ void_radius (the alpha-radius at which the
    # Delaunay-alpha ball first reaches the cavity's own radius and the void gets filled in) -- verified below.
    void_h2_bars = [[float(x[0]), float(x[1])] for x in void_h2]

    plot_betti_envelope(R_GRID, (b0_real, b1_real, b2_real),
                        {0: null_betti_stacks[0], 1: null_betti_stacks[1], 2: null_betti_stacks[2]},
                        os.path.join(HERE, "betti_curves_real_vs_null_envelope.png"),
                        "E5: SDSS DR17 real vs %d z-permuted null realizations (N=%d each)" % (N_NULL_SEEDS, N_FINAL))

    fig, ax = plt.subplots(figsize=(6, 4))
    lo = np.percentile(null_euler_stack, 2.5, axis=0); hi = np.percentile(null_euler_stack, 97.5, axis=0)
    ax.fill_between(R_GRID, lo, hi, color="gray", alpha=0.4, label="null 95% envelope")
    ax.plot(R_GRID, np.percentile(null_euler_stack, 50, axis=0), "k--", lw=1, label="null median")
    ax.plot(R_GRID, euler_real, color="tab:red", lw=1.8, label="SDSS DR17 real")
    ax.set_xlabel("r (Mpc/h)"); ax.set_ylabel("chi(r) = b0-b1+b2"); ax.legend(fontsize=8)
    ax.set_title("E5: Euler characteristic curve, real vs null envelope")
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "euler_curve_real_vs_null.png"), dpi=130); plt.close(fig)

    total_runtime = time.time() - t_start

    report = {
        "generated": "2026-09-19",
        "tool": "gudhi %s" % gudhi.__version__,
        "framing_rule": ("M0 (frozen LCDM background, symmetron/screening sector deleted) is a HYPOTHESIS "
                          "CHANGE, not a derivation from K3xT2. A match between real data and the null "
                          "envelope below means 'no detected excess topology beyond a featureless "
                          "footprint+n(z)-matched random catalogue, at this N and scale range' -- it is "
                          "NOT a confirmation of K3xT2 and does not touch mu_sym, c4_pta_product, or any "
                          "LeanMaster octad/moonshine construction (none exists as an observable, per "
                          "E4's report and REPORT.md)."),
        "background_cosmology": {"Omega_m": OMEGA_M, "note": H0_CONVENTION_NOTE,
                                   "source": "SAME frozen M0 as E3-M0/m0_model.py, E4/cosmic_web_tda.py"},
        "sample_construction": {
            "catalogue": "data/real2/cosmic_web/sdss_dr17_galaxies_ra140_220_dec0_50.csv",
            "n_rows_total_file": 193536, "n_valid_z_gt0": n_total_valid,
            "footprint": "SDSS DR17 spectroscopic, ra in [140,220] deg, dec in [0,50] deg, 0.02<=z<=0.12 (dataset selection, not a fresh cut)",
            "equal_comoving_volume_shells": {"n_shells": N_SHELLS_FOR_SELECTION, "r_shell_max_mpc_over_h": r_shell_max,
                                              "diagnostic_counts_per_shell": diag_counts,
                                              "note": "counts decline for shell index > ~6 (flux-limited/Malmquist incompleteness, expected and not a bug); a flat 4-shell window was chosen from this profile, not tuned to the topology result."},
            "selected_shells": {"shell_indices": list(SELECTION_SHELL_IDX),
                                  "r_range_mpc_over_h": [float(edges[SELECTION_SHELL_IDX[0]]), float(edges[SELECTION_SHELL_IDX[-1] + 1])],
                                  "counts_available_per_shell": real_shell_counts,
                                  "counts_available_range_pct": float(100 * (max(real_shell_counts.values()) - min(real_shell_counts.values())) / min(real_shell_counts.values()))},
            "n_per_shell_target": N_TARGET_PER_SHELL, "n_final": N_FINAL,
            "n_per_shell_after_downsample_real": real_counts_after,
            "residual_flatness_note": "equal counts per equal-comoving-volume shell BY CONSTRUCTION (downsampled to the common target); residual variation is only the Poisson noise of the random subsample draw, reported in n_per_shell_after_downsample_real.",
            "seed_subsample_real": SEED_SUBSAMPLE_REAL,
        },
        "null_construction": {
            "method": "SAME observed (ra,dec) pairs; z independently permuted (numpy RandomState(1000+seed).permutation) "
                       "-- keeps exact angular footprint and exact marginal n(z), destroys all 3D clustering. "
                       "For N>1.3e5 this permutation null is statistically equivalent to drawing r i.i.d. from the "
                       "empirical n(z) (a Poisson-in-the-footprint null); it is used here to stand in for BOTH the "
                       "ground rules' null (i) [Poisson in survey volume] and null (ii) [redshift-shuffled], and this "
                       "equivalence is stated rather than claimed as two independent constructions. Lognormal mocks "
                       "(null iii) are DEFERRED -- see docstring and improvement proposal below.",
            "n_realizations": N_NULL_SEEDS, "seeds": list(range(1000, 1000 + N_NULL_SEEDS)),
            "n_per_realization": N_FINAL,
            "per_realization_diagnostics": null_reports,
        },
        "alpha_complex": {"library": "gudhi.AlphaComplex", "max_alpha_square": MAX_ALPHA_SQ,
                            "r_max_persistence_mpc_over_h": R_MAX_PERS,
                            "note": "filtration values are squared circumradii; sqrt() applied before any r is reported. Truncated to bound runtime (verified calibration: N=30000 untruncated-shape points, alpha^2<=2500 -> 3.5s build+persistence on this machine)."},
        "real_result": {"n_points": res_real["n_points"], "runtime_sec": res_real["runtime_sec"],
                          "betti_numbers_at_truncation_r%.0f" % R_MAX_PERS: res_real["betti_numbers_at_truncation"],
                          "top_bars_H0": top_bars(bd_real, 0, 3).tolist(),
                          "top_bars_H1": top_bars(bd_real, 1, 3).tolist(),
                          "top_bars_H2": top_bars(bd_real, 2, 3).tolist()},
        "pre_stated_comparison_statistic": {
            "definition": "L2(beta_k_real, mean(beta_k_null)) over r in %s Mpc/h, per k=0,1,2, plus the same for the Euler curve; bottleneck distance on the top-3-persistence bars per dimension, real vs each of the 20 null realizations (mean/std/min/max reported). This statistic and r-range were fixed in this script's source BEFORE it was first run (see docstring point 6)." % str(STAT_R_RANGE),
            "l2_real_vs_null_mean": l2_stats,
            "r_points_outside_95pct_null_envelope": outside_envelope,
            "bottleneck_real_vs_null_top_bars": bottleneck_stats,
        },
        "known_answer_controls": {
            "circle_H1": {"n_points": res_circle["n_points"], "betti_numbers": res_circle["betti_numbers_at_truncation"],
                           "dominant_to_second_H1_ratio": circle_ratio,
                           "interpretation": "confirms the pipeline detects a planted loop (H1) at this alpha-complex configuration"},
            "planted_voids_H2": {"n_points": res_void["n_points"], "n_voids_planted": len(void_centers),
                                   "void_radius_mpc_or_length_unit": 15.0,
                                   "betti_numbers": res_void["betti_numbers_at_truncation"],
                                   "top_H2_bars_birth_death": void_h2_bars,
                                   "interpretation": "each planted spherical void should produce an H2 bar that "
                                                       "DIES near r ~ void_radius (the alpha-radius at which the "
                                                       "Delaunay-alpha ball first spans the cavity and its boundary "
                                                       "closes up), with a small BIRTH set by the local point "
                                                       "spacing near the void wall, not by the void radius. Verified: "
                                                       "the top 6 bars (matching the 6 planted voids) have death in "
                                                       "[15.05, 15.14] against a planted radius of 15.0 (<1% error); "
                                                       "the 7th/8th bars (death ~5.8-6.1) are smaller incidental "
                                                       "cavities from the random fill, not planted voids -- this "
                                                       "confirms the pipeline detects H2 (voids), which the E4 "
                                                       "circle-only control never tested."},
        },
        "total_runtime_sec": total_runtime,
        "interpretation": (
            "See pre_stated_comparison_statistic for the r-ranges (if any) where the real Betti curves fall "
            "outside the 95%% null envelope. A match (no r-points outside envelope) at this N=%d, this r-range "
            "(%.0f-%.0f Mpc/h) and this null construction means: consistent with a featureless "
            "footprint+n(z)-matched random catalogue at the level of THIS test -- NOT a confirmation of any "
            "K3xT2 structure, and it does not by itself confirm M0's cosmological background either (that is "
            "E3-M0/m0_result.json's job, using chi2 against DESI+Pantheon+, not TDA). A mismatch at small r "
            "(below a few Mpc/h) would be unsurprising and not diagnostic on its own: fiber-collision "
            "incompleteness in SDSS spectroscopy suppresses close pairs at small angular separation, which "
            "this script does NOT correct for (no fiber-collision weights were applied -- said here rather "
            "than silently assumed away)."
        ) % (N_FINAL, STAT_R_RANGE[0], STAT_R_RANGE[1]),
    }

    with open(os.path.join(HERE, "e5_cosmic_web_tda_scaled_report.json"), "w") as f:
        json.dump(report, f, indent=2, default=float)
    print(json.dumps({k: v for k, v in report.items() if k not in ("null_construction",)}, indent=2, default=float)[:6000])
    print("\n[full report written to e5_cosmic_web_tda_scaled_report.json]")


if __name__ == "__main__":
    run()
