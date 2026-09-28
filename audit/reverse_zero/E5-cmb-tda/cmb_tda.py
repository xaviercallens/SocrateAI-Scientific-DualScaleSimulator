#!/usr/bin/env python
"""
E5: CMB TDA (TDA T2 from the ground rules) -- persistent homology of the WMAP
9yr ILC temperature map, tested against a Gaussian-random-field null generated
from the SAME map's own masked power spectrum, plus a cosmic-string (Kaiser-
Stebbins step) injection sensitivity scan.

DATA:
  map : data/real2/cmb/wmap_ilc_9yr_v5.fits                (WMAP 9yr ILC, nside=512, mK)
  mask: data/real2/cmb/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits (binary, nside=512)

PIPELINE
  1. Read map+mask at native nside=512. ud_grade both to NSIDE_WORK=128
     (hp.ud_grade, power-preserving degrade). Mask thresholded at >=0.5 after
     ud_grade to recover a binary mask (stated). No additional smoothing
     beyond the ILC's native effective beam.
  2. Filtration: GUDHI SimplexTree on the unmasked-pixel adjacency graph
     (vertices, edges from hp.get_all_neighbours restricted to unmasked
     pairs, triangles = 3-cliques among mutual neighbours). This mask-only
     topology (edge/triangle index lists) is built ONCE per process and
     reused for every sim and both filtration directions -- profiling
     showed the neighbour/triangle *search* (set intersections) is
     temperature-independent and dominates cost (~8s of ~9.7s per build);
     caching it cuts wall-clock ~2.5x with IDENTICAL topology and
     filtration-value definitions (see build_topology /
     betti_curves_from_topology; this is a performance-only refactor of
     the smoke-test's build_complex_and_betti_curve, verified to give the
     same edge/triangle counts on the same mask).
     Lower-star filtration by T/sigma gives the SUBLEVEL curve; by -T/sigma
     the SUPERLEVEL (hot-spot) curve. Both are now null-tested (v2 fix:
     the smoke-test script computed a superlevel curve for the data only,
     with no ensemble behind it -- a "consistency check" with no p-value;
     the full run below null-tests both directions properly).
  3. Betti curves B0(nu), B1(nu) and Euler characteristic chi(nu)=B0-B1 are
     read off on a grid nu in [-4,4] sigma, both directions.
  4. NULL: N_SIMS Gaussian sims via healpy.synfast, seeded BASE_SEED+offset+i,
     using the pseudo-Cl/fsky spectrum estimated by hp.anafast on the SAME
     masked ILC map (stated: NOT CAMB, NOT MASTER-deconvolved). v2 fix: the
     monopole and dipole terms (ell=0,1) of that Cl estimate are explicitly
     zeroed before synfast -- a masked-map anafast estimate carries
     mask-induced residual power at ell=0,1 that is not physical CMB
     anisotropy and that the data pipeline itself removes (mean-subtracted
     before anafast); leaving it in the sim generator would inject a
     monopole/dipole into every sim that the data map does not have, biasing
     low-nu Betti values. Each sim uses the SAME mask and the SAME topology
     as the data (point 2); resolution is matched by construction since Cl
     is estimated at, and synfast draws directly at, NSIDE_WORK -- no
     second ud_grade is applied to sims (the smoke-test docstring's claim
     of a matching "ud_grade" step for sims was inaccurate and is corrected
     here: an ud_grade would be a redundant, lossy resample of an
     already-nside_work power spectrum, not a pipeline match). A
     variance-ratio calibration check (sim ensemble vs data, computed
     BEFORE any p-value) is still reported.
  5. STATISTIC: per curve (b0/b1/chi, sublevel AND superlevel), the coarse
     (N_BINS-bin) data vector's Hartlap-corrected Mahalanobis chi2 against
     the null ensemble is reported with its p-value, plus a covariance-free
     empirical rank p-value from a leave-one-out chi2 among the null sims
     themselves. v2 fix: the FULL (uncoarsened) ensemble mean and std curve
     and the standardized residual (data-mean)/std at every nu are now
     saved -- without them a chi2 alone cannot distinguish a broadband
     amplitude offset (spectrum mismatch) from a localized topological
     excursion.
  6. SENSITIVITY: a toy Kaiser-Stebbins cosmic-string injection (temperature
     step Delta T/T = 8 pi G mu v_gamma per crossing, v_gamma=0.6 fixed) is
     added to fresh Gaussian sims at several Gmu values; the smallest Gmu at
     which the injected ensemble separates from the pure-Gaussian null at
     95% is reported as a TDA-statistic sensitivity estimate -- NOT a
     cosmological bound, and NOT claimed equal to LeanMaster's
     planckGmuBound = 1.5e-7 (PRE_REGISTRATION.md:66: "nothing: it is an
     imported upper bound, not a prediction"). v2 fixes:
       (a) the 95th-percentile threshold is now taken from the LEAVE-ONE-OUT
           chi2 distribution of the null ensemble (each sim scored against
           the OTHER sims' mean/cov), not the in-sample chi2 (each sim
           scored partly against itself) -- the in-sample version is biased
           low and would report an artificially small (falsely sensitive)
           Gmu;
       (b) the injected "segment" is now truncated to a genuine arc (a
           random angular half-length in [1,10] deg along the great circle,
           in addition to the pre-existing ~0.5 deg perpendicular band) --
           the smoke-test code applied the temperature step across the
           WHOLE great circle (no along-circle truncation), which is far
           more coherent than a string segment and would make the
           sensitivity estimate optimistic by an unknown factor. Documented
           as a toy, non-network model regardless (see inject_ks_strings
           docstring).
       DISCRIMINATING CHECK stated in the ground rules' spirit: if this
       toy-model scan on an nside=128 WMAP ILC map returns a threshold
       BELOW 1.5e-7, that is evidence the toy is still too optimistic (a
       degraded WMAP map cannot out-sensitivity Planck-grade analyses), not
       evidence of new physics; the result section says so explicitly.

FRAMING (non-negotiable, per ground rules): M0 predicts Gaussian statistics
(no topological defects). A null result (data consistent with the Gaussian
ensemble) is consistent with M0 and with plain LCDM EQUALLY -- it is not
evidence for K3xT2, and changes the free-parameter count by exactly zero.

CHUNKING (v2, added for the N>=100 full run under the 10-minute-per-command
memory-safety cap): the mask-only topology cache (~5.8s) plus ~3.7s/sim/
direction means N=100 null sims x 2 directions ~ 100*2*3.7+6 ~ 746s, over the
600s cap in one call. Use:
  --dump-npz PATH             after generating this chunk's null sims, save
                               raw per-sim curves (not just the coarse chi2)
                               so multiple chunks can be concatenated exactly
                               (mean/cov/p-values computed once on the FULL
                               pooled ensemble, not chunk-by-chunk).
  --seed-offset N              shift this chunk's seeds so chunks are disjoint.
  --skip-null --null-cache F   skip null generation, load a merged null
                               ensemble (from merge_null.py) for the
                               sensitivity scan's reference mean/cov.
See merge_null.py (pools --dump-npz chunks into one null ensemble + report)
and merge_final.py (combines the null report with a sensitivity report).

Run (smoke test, already done, see smoke_test.json -- N=10, pre-v2-fix
script, kept for provenance):
  .venv-tda/bin/python audit/reverse_zero/E5-cmb-tda/cmb_tda.py --n-sims 10 --out smoke_test.json
"""
import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np
import healpy as hp
import gudhi

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]  # worktree root
MAP_PATH = REPO / "data/real2/cmb/wmap_ilc_9yr_v5.fits"
MASK_PATH = REPO / "data/real2/cmb/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits"

NSIDE_WORK = 128
NU_MIN, NU_MAX, NU_STEP_GRID = -4.0, 4.0, 41   # fine grid for the reported curve
N_BINS = 8                                      # coarse bins for chi2 (Hartlap-safe with N_SIMS=100)
BASE_SEED = 20260919


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def load_data_map_and_mask():
    m = hp.read_map(str(MAP_PATH), field=0)
    mask = hp.read_map(str(MASK_PATH), field=0)
    assert hp.get_nside(m) == hp.get_nside(mask) == 512
    m_dg = hp.ud_grade(m, NSIDE_WORK)
    mask_dg_frac = hp.ud_grade(mask, NSIDE_WORK)
    mask_dg = (mask_dg_frac >= 0.5).astype(np.uint8)
    return m_dg, mask_dg, mask_dg_frac


def estimate_cl(masked_map, mask, nside, lmax=None):
    """Pseudo-Cl / fsky approximation from the data's own masked map.
    v2: explicitly zero ell=0,1 (monopole/dipole) -- a masked-map anafast
    estimate carries mask-induced residual power there that is not physical
    CMB anisotropy; the data itself is mean-subtracted before anafast, so
    leaving ell=0,1 in the Cl used for synfast would inject a monopole/
    dipole into every sim that the data does not have."""
    if lmax is None:
        lmax = 2 * nside
    fsky = mask.mean()
    m_for_anafast = np.where(mask > 0, masked_map - masked_map[mask > 0].mean(), 0.0)
    cl = hp.anafast(m_for_anafast, lmax=lmax)
    cl_corrected = (cl / max(fsky, 1e-6)).copy()
    cl_corrected[0] = 0.0
    cl_corrected[1] = 0.0
    return cl_corrected, fsky


def build_topology(mask, nside):
    # =====================================================================
    # DEFECT D2 -- FOUND 2026-09-19 -- DO NOT USE FOR NEW WORK
    # ---------------------------------------------------------------------
    # This function builds a MALFORMED 2-complex.  It inserts every 3-clique
    # of the 8-neighbour pixel adjacency graph.  Each group of four pixels
    # meeting at a HEALPix grid corner is a 4-clique, and filling all four of
    # its triangles makes a hollow tetrahedron, i.e. a 2-sphere.
    #
    # Measured by the known-answer suite (branch loop/tda-simple, case F3topo
    # in audit/tda_validation/simple_suite/results.json), full sky, nside 64:
    #     V = 49152, E = 196596, F = 196592, V - E + F = 49148
    #     Betti = (1, 0, 49147)        expected for a triangulated S2: (1,0,1)
    # At nside 32: V - E + F = 12284, Betti = (1, 0, 12283).
    # With a 30 deg polar-cap mask at nside 32 it gives Betti (1, 0, 761) for
    # what is a disc, whose Betti numbers are (1, 0, 0).
    #
    # FIXED REPLACEMENT:
    #     audit/tda_validation/tda_fixed/cmb_topology.py
    #         :: build_topology_fixed(mask, nside)
    #   Same call shape.  Builds the HEALPix quad-corner 2-complex instead of
    #   filling cliques: V - E + F = 2 and Betti (1, 0, 1) at nside 8, 16, 32
    #   and 64; a polar cap gives (1, 0, 0).  Before/after numbers in
    #   audit/tda_validation/tda_fixed/validation_results.json; regression
    #   tests in audit/tda_validation/tda_fixed/tests/test_cmb_topology.py.
    #
    # THIS FUNCTION IS LEFT UNCHANGED ON PURPOSE.  It is the committed record
    # of the runs in audit/reverse_zero/ and other branches' results were
    # produced with it.  Only this comment block was added.
    # =====================================================================
    """Mask-only (temperature-independent) adjacency: unmasked pixel list,
    edge index-pairs, triangle index-triples among mutual neighbours. Built
    ONCE and reused for every sim and both filtration directions -- this is
    a performance-only refactor (verified against the smoke-test's per-call
    rebuild: identical edge_set/triangle count, identical Betti curves on
    the same input); it changes no result, only wall-clock (~2.5x faster,
    measured)."""
    npix = hp.nside2npix(nside)
    unmasked = np.where(mask > 0)[0]
    idx_map = -np.ones(npix, dtype=np.int64)
    idx_map[unmasked] = np.arange(unmasked.size)
    neighbours = hp.get_all_neighbours(nside, unmasked)
    edge_set = set()
    for k in range(8):
        nb = neighbours[k]
        valid = nb >= 0
        nb_idx = idx_map[nb[valid]]
        src_idx = np.arange(unmasked.size)[valid]
        keep = nb_idx >= 0
        for a, b in zip(src_idx[keep], nb_idx[keep]):
            if a == b:
                continue
            edge_set.add((min(int(a), int(b)), max(int(a), int(b))))
    edges_arr = np.array(sorted(edge_set), dtype=np.int64) if edge_set else np.zeros((0, 2), dtype=np.int64)
    adj = {}
    for a, b in edge_set:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    tris = []
    for a, b in edge_set:
        common = adj.get(a, set()) & adj.get(b, set())
        for c in common:
            if c > b:
                tris.append((a, b, c))
    tris_arr = np.array(tris, dtype=np.int64) if tris else np.zeros((0, 3), dtype=np.int64)
    return unmasked, edges_arr, tris_arr


def betti_curves_from_topology(temp, unmasked, edges_arr, tris_arr, nu_grid, sublevel=True):
    # =====================================================================
    # DEFECT D2, THIRD CLAUSE -- FOUND 2026-09-19 -- MISLABELLED RETURN VALUE
    # ---------------------------------------------------------------------
    # The third return value is chi = b0 - b1 and is called the "Euler
    # characteristic" here, in the module docstring (point 3) and in the
    # 'euler_chi' key of every report this file wrote.  For a 2-complex the
    # Euler characteristic is b0 - b1 + b2, so the name is wrong whenever
    # b2 != 0 -- which is ALWAYS on a closed surface, and catastrophically so
    # on the complex build_topology returns (b2 = 49147 at nside 64).
    # Measured (audit/tda_validation/tda_fixed/parts/smooth.json, nside 64,
    # full sky, seed 30000): at the top of the filtration this function
    # returns 1, while the true Euler characteristic of the sphere is 2.
    #
    # The FILTRATION ITSELF IS CORRECT and was verified: on the F1 grid case
    # (256x256 Freudenthal triangulation, 7 Gaussian wells, seed 21) the
    # fixed implementation reproduces this one's b0 and b1 curves EXACTLY
    # (np.array_equal, 4001 grid points).  Only the Euler key is wrong, and
    # the complex it is usually fed (see build_topology above).
    #
    # FIXED REPLACEMENT:
    #     audit/tda_validation/tda_fixed/cmb_topology.py
    #         :: betti_curves_from_topology(...)   -- same signature
    #   Returns a dict: 'b0', 'b1', 'euler_char_true' (the true #V-#E+#F of
    #   the filtered subcomplex), 'b0_minus_b1' (what this function returns)
    #   and 'b2_implied'.
    #
    # THIS FUNCTION IS LEFT UNCHANGED ON PURPOSE (committed record of past
    # runs).  Only this comment block was added.
    # =====================================================================
    """Given the cached topology, compute Betti-0/1 and Euler-char curves for
    ONE map on nu_grid, for the given filtration direction. Same lower-star
    filtration definition (max of endpoint T/sigma) as the smoke-test's
    per-call build_complex_and_betti_curve."""
    sigma = temp[unmasked].std()
    t = temp / sigma
    if not sublevel:
        t = -t
    verts_f = t[unmasked]

    st = gudhi.SimplexTree()
    for i in range(unmasked.size):
        st.insert([i], filtration=float(verts_f[i]))
    if edges_arr.size:
        ef = np.maximum(verts_f[edges_arr[:, 0]], verts_f[edges_arr[:, 1]])
        for i in range(edges_arr.shape[0]):
            st.insert([int(edges_arr[i, 0]), int(edges_arr[i, 1])], filtration=float(ef[i]))
    if tris_arr.size:
        trf = np.maximum(np.maximum(verts_f[tris_arr[:, 0]], verts_f[tris_arr[:, 1]]), verts_f[tris_arr[:, 2]])
        for i in range(tris_arr.shape[0]):
            st.insert([int(tris_arr[i, 0]), int(tris_arr[i, 1]), int(tris_arr[i, 2])], filtration=float(trf[i]))

    st.make_filtration_non_decreasing()
    st.compute_persistence(persistence_dim_max=True)
    diag0 = st.persistence_intervals_in_dimension(0)
    diag1 = st.persistence_intervals_in_dimension(1)

    b0 = np.array([np.sum((diag0[:, 0] <= nu) & ((diag0[:, 1] > nu) | np.isinf(diag0[:, 1]))) for nu in nu_grid])
    b1 = np.array([np.sum((diag1[:, 0] <= nu) & ((diag1[:, 1] > nu) | np.isinf(diag1[:, 1]))) for nu in nu_grid]) if diag1.size else np.zeros_like(nu_grid, dtype=int)
    chi = b0 - b1
    return b0, b1, chi


def coarsen(curve, n_bins):
    n = len(curve)
    edges = np.linspace(0, n, n_bins + 1).astype(int)
    return np.array([curve[edges[i]:edges[i + 1]].mean() for i in range(n_bins)])


def chi2_stat(vec, mean, cov_inv):
    d = vec - mean
    return float(d @ cov_inv @ d)


def coarse_stats(sim_curves, data_curve, n_bins=N_BINS):
    # =====================================================================
    # DEFECT D1 -- FOUND 2026-09-19 -- MISCALIBRATED p-VALUES
    # ---------------------------------------------------------------------
    # The chi2 p-values this function returns are not uniform under the null.
    # Measured by the known-answer suite (branch loop/tda-simple, case F3 in
    # audit/tda_validation/simple_suite/results.json): 100 independent test
    # maps against a 100-map same-C_ell ensemble, HEALPix nside 64, KS test
    # of the 100 p-values against U(0,1):
    #     b0   KS p = 0.0205        #{p < 0.05} = 3
    #     b1   KS p = 0.000109      #{p < 0.05} = 9
    #     chi  KS p = 0.397         #{p < 0.05} = 9
    # Three causes, all measured:
    #   (a) df is ALWAYS n_bins (8) below, even when a coarse bin is
    #       identically constant across the whole ensemble (b0 bin 7 = 1.0 in
    #       100/100 sims; b1 bin 0 = 0.0 in 100/100).
    #   (b) `cov = np.cov(...) + 1e-8 * np.eye(n_bins)` followed by
    #       `np.linalg.pinv` INVERTS the ridge instead of dropping the dead
    #       direction, giving that direction a weight of about 1e8.
    #   (c) bins that are atomic but not dead are treated as Gaussian: b1
    #       bin 1 is 0.0 in 74/100 sims (std 0.1202), so a single extra count
    #       there gives |z| ~ 6 and a chi2 survival p ~ 1e-6.
    # The same three causes apply to the leave-one-out branch below, which
    # the sensitivity scan uses for its 95th-percentile threshold.
    #
    # FIXED REPLACEMENT:
    #     audit/tda_validation/tda_fixed/stats.py
    #         :: coarse_stats_fixed(sim_curves, data_curve, n_bins=...)
    #   Drops dead and atomic bins using the ensemble only, conditions the
    #   kept block by eigenvalue truncation (no ridge), sets df to the
    #   retained rank, and returns BOTH the chi2 p-value and a
    #   pooled-exchangeable rank p-value, with the kept-bin count and the
    #   dropped-bin list (and whether the test vector differs in a dropped
    #   bin) in the output.  On the same F3 case: b0 KS p = 0.386
    #   (#{p<0.05} = 4), b1 KS p = 0.222 (9), chi KS p = 0.397 (9).
    #   Before/after in audit/tda_validation/tda_fixed/validation_results.json;
    #   regression tests in audit/tda_validation/tda_fixed/tests/test_stats.py.
    #
    # EFFECT ON THIS DIRECTORY'S COMMITTED RESULTS: recomputing null_report
    # .json's six p-values from chunks/null_merged.npz with the fixed
    # statistic moves them, e.g. sublevel b0 0.0650 -> 0.0428 (crosses 0.05)
    # and superlevel b1 0.0109 -> 0.0054.  See validation_results.json,
    # section e5 -- and note that those curves were themselves produced by
    # the defective build_topology above, which is NOT corrected there.
    #
    # THIS FUNCTION IS LEFT UNCHANGED ON PURPOSE (committed record of past
    # runs).  Only this comment block was added.
    # =====================================================================
    """Hartlap-corrected Mahalanobis chi2 of the data's coarse curve against
    the null ensemble, plus a leave-one-out chi2 distribution among the null
    sims (used both as a covariance-free empirical rank p-value AND, by the
    sensitivity scan, as the reference distribution for the 95th-percentile
    threshold -- v2 fix: NOT the in-sample chi2, which is biased low)."""
    sim_curves = np.asarray(sim_curves)
    sim_coarse = np.array([coarsen(c, n_bins) for c in sim_curves])
    data_coarse = coarsen(np.asarray(data_curve), n_bins)
    mean = sim_coarse.mean(axis=0)
    cov = np.cov(sim_coarse, rowvar=False) + 1e-8 * np.eye(n_bins)
    n = sim_coarse.shape[0]
    p = n_bins
    hartlap = (n - p - 2) / (n - 1) if n > p + 2 else float("nan")
    cov_inv = np.linalg.pinv(cov) * (hartlap if not np.isnan(hartlap) else 1.0)
    chi2 = chi2_stat(data_coarse, mean, cov_inv)
    from scipy.stats import chi2 as chi2dist
    pval = float(chi2dist.sf(chi2, df=p)) if not np.isnan(hartlap) else None

    loo_chi2 = []
    for j in range(n):
        others = np.delete(sim_coarse, j, axis=0)
        m2 = others.mean(axis=0)
        c2 = np.cov(others, rowvar=False) + 1e-8 * np.eye(p)
        n2 = others.shape[0]
        hl2 = (n2 - p - 2) / (n2 - 1) if n2 > p + 2 else 1.0
        ci2 = np.linalg.pinv(c2) * hl2
        loo_chi2.append(chi2_stat(sim_coarse[j], m2, ci2))
    loo_chi2 = np.array(loo_chi2)
    rank_p = float((np.sum(loo_chi2 >= chi2) + 1) / (n + 1))

    return {
        "data_chi2_hartlap": chi2, "p_value_chi2_survival": pval,
        "hartlap_factor": hartlap, "n_bins": p, "n_sims": n,
        "empirical_rank_p": rank_p,
        "sim_loo_chi2_mean": float(loo_chi2.mean()), "sim_loo_chi2_std": float(loo_chi2.std()),
        "loo_chi2_values": loo_chi2.tolist(),
        "sim_ensemble_full_curve_mean": sim_curves.mean(axis=0).tolist(),
        "sim_ensemble_full_curve_std": sim_curves.std(axis=0).tolist(),
        "standardized_residual_full_curve":
            ((np.asarray(data_curve) - sim_curves.mean(axis=0)) /
             np.where(sim_curves.std(axis=0) == 0, np.nan, sim_curves.std(axis=0))).tolist(),
    }


def inject_ks_strings(temp_map, nside, mask, n_segments, gmu, v_gamma=0.6, rng=None):
    """
    Toy Kaiser-Stebbins injection. A cosmic string segment is modeled as a
    great-circle ARC on the sphere: a random plane (normal n) sets the step
    discontinuity's orientation; within that plane, a random angular
    half-length in [1,10] deg (uniform(2,20)/2) and a random phase center
    truncate the step to a genuine segment (v2 fix -- the smoke-test code
    applied the step across the WHOLE great circle with no along-circle
    truncation, i.e. it injected a full great-circle step, not a segment;
    that is far more coherent than a string segment and would bias the
    sensitivity estimate optimistic). Each pixel within the along-circle arc
    AND within a +/-0.5 deg perpendicular band of the great circle (localizes
    the step to a thin strip straddling the string, comparable to the
    nside=128 pixel scale of ~0.46 deg) receives Delta T/T = 8 pi G mu *
    v_gamma. This is a documented toy model (no network evolution, no
    correlated segment population, no scaling density of strings): it is
    used only to calibrate the TDA statistic's SENSITIVITY to a KS-like step
    pattern, not to forecast a physical string signal or a bound.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    injected = temp_map.copy()
    T_CMB = 2.7255  # K, standard CMB monopole temperature used to convert dT/T to mK
    step_mK = 8 * np.pi * gmu * v_gamma * T_CMB * 1000.0  # mK
    unmasked = np.where(mask > 0)[0]
    vecs = np.array(hp.pix2vec(nside, unmasked)).T  # (N,3)
    perp_halfwidth = np.radians(0.5)
    sin_perp = np.sin(perp_halfwidth)
    for _ in range(n_segments):
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        u = rng.normal(size=3)
        u -= (u @ n) * n
        u /= np.linalg.norm(u)
        w = np.cross(n, u)
        half_len = np.radians(rng.uniform(2, 20)) / 2.0
        center_phase = rng.uniform(0, 2 * np.pi)

        side = vecs @ n
        in_band = np.abs(side) < sin_perp
        proj_u = vecs @ u
        proj_w = vecs @ w
        phase = np.arctan2(proj_w, proj_u)
        dphase = np.angle(np.exp(1j * (phase - center_phase)))
        in_arc = np.abs(dphase) < half_len
        near_segment = in_band & in_arc

        step_sign = np.sign(side)
        injected[unmasked[near_segment]] += step_mK * step_sign[near_segment]
    return injected, step_mK


def run_null(n_sims, out_path, dump_npz=None, seed_base=BASE_SEED, seed_offset=0, lmax_override=None):
    """Generate n_sims Gaussian null sims (this chunk), compute data + sim
    Betti/Euler curves for BOTH sublevel and superlevel, and (if dump_npz)
    save everything needed to pool this chunk with others."""
    log("loading data map + mask")
    data_map, mask, mask_frac = load_data_map_and_mask()
    fsky = float(mask.mean())
    log(f"ud_graded to nside={NSIDE_WORK}, unmasked pixels={int(mask.sum())}/{mask.size} (fsky={fsky:.4f})")

    cl, fsky_cl = estimate_cl(data_map, mask, NSIDE_WORK, lmax=lmax_override)
    log(f"estimated Cl via anafast/fsky on masked map (ell=0,1 zeroed), lmax={len(cl)-1}, fsky_used={fsky_cl:.4f}")

    nu_grid = np.linspace(NU_MIN, NU_MAX, NU_STEP_GRID)

    log("building mask-only topology cache (once)")
    unmasked, edges_arr, tris_arr = build_topology(mask, NSIDE_WORK)
    log(f"topology: n_unmasked={unmasked.size} n_edges={edges_arr.shape[0]} n_triangles={tris_arr.shape[0]}")

    log("computing data Betti/Euler curves (sublevel + superlevel)")
    b0_d, b1_d, chi_d = betti_curves_from_topology(data_map, unmasked, edges_arr, tris_arr, nu_grid, sublevel=True)
    b0_d_s, b1_d_s, chi_d_s = betti_curves_from_topology(data_map, unmasked, edges_arr, tris_arr, nu_grid, sublevel=False)

    seed_start = seed_base + seed_offset
    log(f"generating {n_sims} Gaussian null sims (synfast, seeds {seed_start}..{seed_start+n_sims-1})")
    sim_b0, sim_b1, sim_chi = [], [], []
    sim_b0_s, sim_b1_s, sim_chi_s = [], [], []
    sim_var = []
    for i in range(n_sims):
        rng_seed = seed_start + i
        np.random.seed(rng_seed)
        sim = hp.synfast(cl, nside=NSIDE_WORK, new=True)
        sim_var.append(float(sim[mask > 0].var()))
        b0, b1, chi = betti_curves_from_topology(sim, unmasked, edges_arr, tris_arr, nu_grid, sublevel=True)
        b0s, b1s, chis = betti_curves_from_topology(sim, unmasked, edges_arr, tris_arr, nu_grid, sublevel=False)
        sim_b0.append(b0); sim_b1.append(b1); sim_chi.append(chi)
        sim_b0_s.append(b0s); sim_b1_s.append(b1s); sim_chi_s.append(chis)
        if (i + 1) % max(1, n_sims // 10) == 0:
            log(f"  sim {i+1}/{n_sims} done")
    sim_b0 = np.array(sim_b0); sim_b1 = np.array(sim_b1); sim_chi = np.array(sim_chi)
    sim_b0_s = np.array(sim_b0_s); sim_b1_s = np.array(sim_b1_s); sim_chi_s = np.array(sim_chi_s)

    data_var = float(data_map[mask > 0].var())
    sim_var_mean = float(np.mean(sim_var))
    var_ratio = sim_var_mean / data_var if data_var else float("nan")
    log(f"CALIBRATION CHECK: data unmasked-pixel variance={data_var:.6g}, "
        f"sim ensemble mean variance={sim_var_mean:.6g}, ratio={var_ratio:.4f}")

    chunk = {
        "n_sims": n_sims, "seeds": list(range(seed_start, seed_start + n_sims)),
        "nu_grid_sigma": nu_grid.tolist(),
        "fsky": fsky, "n_unmasked_pixels_work_res": int(unmasked.size),
        "n_edges": int(edges_arr.shape[0]), "n_triangles": int(tris_arr.shape[0]),
        "lmax_used": len(cl) - 1,
        "data_var": data_var, "sim_var_mean": sim_var_mean, "var_ratio": var_ratio,
        "data_betti_curve": {
            "sublevel": {"b0": b0_d.tolist(), "b1": b1_d.tolist(), "euler_chi": chi_d.tolist()},
            "superlevel": {"b0": b0_d_s.tolist(), "b1": b1_d_s.tolist(), "euler_chi": chi_d_s.tolist()},
        },
    }
    if dump_npz:
        # save the raw arrays FIRST -- these are the expensive-to-recompute
        # result; if the (cheap) summary JSON write below fails for a path
        # reason, the sim ensemble itself is not lost.
        np.savez(dump_npz,
                  sim_b0=sim_b0, sim_b1=sim_b1, sim_chi=sim_chi,
                  sim_b0_super=sim_b0_s, sim_b1_super=sim_b1_s, sim_chi_super=sim_chi_s,
                  b0_d=b0_d, b1_d=b1_d, chi_d=chi_d,
                  b0_d_super=b0_d_s, b1_d_super=b1_d_s, chi_d_super=chi_d_s,
                  nu_grid=nu_grid, seeds=np.array(list(range(seed_start, seed_start + n_sims))),
                  data_var=data_var, sim_var=np.array(sim_var), fsky=fsky,
                  n_unmasked=unmasked.size, n_edges=edges_arr.shape[0], n_tri=tris_arr.shape[0],
                  lmax=len(cl) - 1)
        log(f"dumped raw arrays to {dump_npz}")

    with open(out_path, "w") as f:
        json.dump(chunk, f, indent=2)
    log(f"wrote chunk summary {out_path}")
    return chunk


def run_sensitivity(null_cache_npz, n_segments_inj, gmu_list, n_inj_sims, out_path, seed_base=BASE_SEED, lmax_override=None):
    """Load a merged null ensemble (sim_chi, sublevel) from null_cache_npz,
    build its reference mean/cov (LOO-based 95th-percentile threshold --
    v2 fix, see coarse_stats), then inject KS-string steps into FRESH
    Gaussian sims (same Cl re-derived from the data map, same topology) at
    each Gmu and report the fraction exceeding the LOO-95th-percentile
    threshold."""
    cache = np.load(null_cache_npz)
    sim_chi = cache["sim_chi"]
    chi_d = cache["chi_d"]
    nu_grid = cache["nu_grid"]

    log("loading data map + mask (for Cl + topology, sensitivity run)")
    data_map, mask, mask_frac = load_data_map_and_mask()
    cl, fsky_cl = estimate_cl(data_map, mask, NSIDE_WORK, lmax=lmax_override)
    unmasked, edges_arr, tris_arr = build_topology(mask, NSIDE_WORK)

    null_stats = coarse_stats(sim_chi, chi_d, n_bins=N_BINS)
    loo_chi2 = np.array(null_stats["loo_chi2_values"])
    thresh_95 = float(np.percentile(loo_chi2, 95))
    log(f"LOO-based null 95th percentile chi2 threshold = {thresh_95:.3f} "
        f"(v2 fix: NOT the in-sample percentile, which is biased low)")

    sim_coarse = np.array([coarsen(c, N_BINS) for c in sim_chi])
    mean_chi = sim_coarse.mean(axis=0)
    cov_chi = np.cov(sim_coarse, rowvar=False) + 1e-8 * np.eye(N_BINS)
    n_null = sim_coarse.shape[0]
    hartlap_null = (n_null - N_BINS - 2) / (n_null - 1) if n_null > N_BINS + 2 else 1.0
    cov_inv_chi = np.linalg.pinv(cov_chi) * hartlap_null

    log(f"running cosmic-string sensitivity scan: Gmu = {gmu_list}, n_segments={n_segments_inj}, n_inj_sims={n_inj_sims}")
    sens = {}
    for gmu in gmu_list:
        # v2 fix: the injection-GEOMETRY rng is seeded IDENTICALLY for every
        # Gmu (only the Gaussian background draw and the step amplitude
        # differ across Gmu, both already controlled -- background by
        # np.random.seed(seed_base+5000+i) below, amplitude by gmu itself).
        # Previously this rng was seeded differently per Gmu
        # (1000 + round(log10(gmu)*-1000)), so each Gmu scanned a DIFFERENT
        # set of string segment positions on the same backgrounds -- a
        # non-monotonic fraction_exceeding across Gmu would then be
        # confounded between amplitude and geometry, not just amplitude.
        rng = np.random.default_rng(7777)
        inj_chi2 = []
        step_mK = None
        for i in range(n_inj_sims):
            np.random.seed(seed_base + 5000 + i)
            sim = hp.synfast(cl, nside=NSIDE_WORK, new=True)
            sim_inj, step_mK = inject_ks_strings(sim, NSIDE_WORK, mask, n_segments_inj, gmu, rng=rng)
            _, _, chi_i = betti_curves_from_topology(sim_inj, unmasked, edges_arr, tris_arr, nu_grid, sublevel=True)
            chi_coarse_i = coarsen(chi_i, N_BINS)
            inj_chi2.append(chi2_stat(chi_coarse_i, mean_chi, cov_inv_chi))
        inj_chi2 = np.array(inj_chi2)
        frac_separated = float(np.mean(inj_chi2 > thresh_95))
        sens[str(gmu)] = {
            "n_injected_sims": n_inj_sims, "n_segments": n_segments_inj,
            "step_mK_example": step_mK,
            "fraction_exceeding_null_95th_pct_LOO": frac_separated,
            "null_95th_pct_chi2_LOO": thresh_95,
            "mean_injected_chi2": float(inj_chi2.mean()),
            "injected_chi2_values": inj_chi2.tolist(),
        }
        log(f"  Gmu={gmu}: frac exceeding LOO null 95th pct = {frac_separated:.2f}")

    result = {
        "null_reference_stats_euler_chi_sublevel": null_stats,
        "cosmic_string_sensitivity": {
            "model": "toy Kaiser-Stebbins: great-circle ARC segments (v2: truncated "
                     "along-circle to a random 1-10 deg half-length, not a full "
                     "great circle), Delta T/T = 8 pi G mu v_gamma per crossing, "
                     "v_gamma=0.6 fixed -- NOT a scaling-solution network -- "
                     "sensitivity of the TDA statistic to this injection only, "
                     "not a cosmological bound",
            "statistic_used": "euler_chi coarse chi2 vs LOO null 95th percentile (v2 fix)",
            "planck_gmu_bound_comparison": "LeanMaster planckGmuBound = 1.5e-7 "
                "(PRE_REGISTRATION.md:66: 'nothing: it is an imported upper bound, "
                "not a prediction'). Any Gmu found here is a TDA-pipeline "
                "sensitivity threshold for THIS toy injection on a degraded "
                "(nside=128) WMAP ILC map, compared in order of magnitude only, "
                "never claimed equal or as a new bound. DISCRIMINATING CHECK: a "
                "threshold below 1.5e-7 here would indicate the toy model is "
                "still optimistic (a degraded WMAP map should not out-sensitivity "
                "Planck-grade analyses), not new physics.",
            "gmu_scan": sens,
        },
    }
    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    log(f"wrote {out_path}")
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-sims", type=int, default=100)
    ap.add_argument("--seed-offset", type=int, default=0)
    ap.add_argument("--dump-npz", type=str, default=None)
    ap.add_argument("--skip-null", action="store_true")
    ap.add_argument("--null-cache", type=str, default=None)
    ap.add_argument("--n-segments-inj", type=int, default=30)
    ap.add_argument("--n-inj-sims", type=int, default=20)
    ap.add_argument("--gmu", type=str, default="1e-6,1e-7,1e-8,1e-9")
    ap.add_argument("--no-sensitivity", action="store_true")
    ap.add_argument("--out", type=str, default="e5_cmb_tda_report.json")
    ap.add_argument("--lmax", type=int, default=None)
    args = ap.parse_args()
    out_path = HERE / args.out

    if args.skip_null:
        assert args.null_cache, "--skip-null requires --null-cache"
        gmu_list = [float(x) for x in args.gmu.split(",")]
        run_sensitivity(args.null_cache, args.n_segments_inj, gmu_list, args.n_inj_sims, out_path, lmax_override=args.lmax)
    else:
        run_null(args.n_sims, out_path, dump_npz=args.dump_npz, seed_offset=args.seed_offset, lmax_override=args.lmax)
