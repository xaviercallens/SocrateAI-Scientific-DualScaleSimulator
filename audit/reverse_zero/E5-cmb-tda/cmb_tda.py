#!/usr/bin/env python
"""
E5: CMB TDA (TDA T2 from the ground rules) -- persistent homology of the WMAP
9yr ILC temperature map, tested against a Gaussian-random-field null generated
from the SAME map's own masked power spectrum, plus a cosmic-string (Kaiser-
Stebbins step) injection sensitivity scan.

DATA (from the fetched, hash-verified manifest -- see MANIFEST.json in the
worktree root, ids wmap_ilc_9yr_v5 / wmap_kq85_mask):
  map : data/real2/cmb/wmap_ilc_9yr_v5.fits                (WMAP 9yr ILC, nside=512, mK)
  mask: data/real2/cmb/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits (binary, nside=512)

PIPELINE
  1. Read map+mask at native nside=512. ud_grade (power-preserving degrade,
     hp.ud_grade with pess=False default) both to NSIDE_WORK=128 -- stated
     explicitly. The mask, a binary map, becomes fractional after ud_grade;
     threshold at >=0.5 to recover a binary mask at the working resolution
     (stated). No additional smoothing is applied beyond the ILC's native
     ~1 deg effective beam (WMAP ILC is a foreground-cleaned combination, not
     a single-beam map; we do not deconvolve or reconvolve it -- "no
     additional smoothing, ILC native effective beam" per the ground rules).
  2. Filtration: GUDHI SimplexTree built from the unmasked-pixel adjacency
     graph (vertices = unmasked HEALPix pixels, edges from
     hp.get_all_neighbours restricted to unmasked pairs, triangles = closed
     3-cliques among mutual neighbours so H1 is a genuine 2-complex
     invariant, not a graph-only artifact). Lower-star filtration by
     T/sigma (sigma = masked-map std) gives the SUBLEVEL filtration;
     lower-star filtration by -T/sigma gives the SUPERLEVEL filtration
     (hot-spot topology). H0/H1 persistence computed with GUDHI for both.
  3. Betti curves B0(nu), B1(nu) and the Euler characteristic curve
     chi(nu) = B0(nu) - B1(nu) are read off the persistence diagrams on a
     grid of nu in [-4,4] sigma (NU_GRID points), for the sublevel
     filtration only (the superlevel filtration of the same map gives the
     mirror-image curve by construction and is reported as a consistency
     check, not a second independent statistic).
  4. NULL: N_SIMS Gaussian sims via healpy.synfast, seeded 20260919+i, using
     the power spectrum C_ell estimated by hp.anafast on the SAME masked ILC
     map, corrected for the mask's f_sky (a standard, stated
     pseudo-C_ell/fsky approximation -- NOT MASTER-deconvolved; this is
     named explicitly as an approximation, per the ground rules' NULL
     clause "or the CAMB M0 spectrum; say which" -- we use the data-derived
     anafast/fsky spectrum, not CAMB, because it best matches the ILC map's
     own effective transfer function/residual-foreground envelope without a
     model instrument response; this choice is stated in the output JSON).
     Each sim gets the SAME mask, SAME ud_grade, SAME pipeline as the data,
     and (calibration check, run BEFORE any p-value is computed) the
     unmasked-pixel variance of the sim ensemble is compared to the data's.
  5. STATISTIC: per nu-bin, the vector of (B0,B1,chi) curves is coarsened to
     N_BINS bins; the sample covariance across the N_SIMS null sims is
     computed per curve, Hartlap-corrected ((N-p-2)/(N-1)), and a
     Mahalanobis chi2 of the data curve against the null mean is reported
     with the corresponding p-value (chi2 survival function, dof=N_BINS).
     An empirical rank p-value (fraction of null sims with a summary
     statistic >= the data's, using a leave-one-out chi2 against the
     other N_SIMS-1 sims for each null member) is reported alongside as a
     covariance-free cross-check, per the ground rules' emphasis on not
     reporting a single fragile number.
  6. SENSITIVITY: a toy Kaiser-Stebbins cosmic-string injection (straight
     segments on the sphere, temperature step Delta T/T = 8 pi G mu v gamma
     per crossing, v gamma ~ 0.6 fixed, random great-circle segments) is
     added to fresh Gaussian sims at several Gmu values; the same statistic
     is recomputed and the smallest Gmu at which >=95% of injected sims separate
     from the pure-Gaussian null ensemble (by the chi2 statistic's 95th
     percentile) is reported as a TDA-statistic sensitivity estimate -- NOT
     a cosmological bound, and NOT claimed equal to LeanMaster's
     planckGmuBound = 1.5e-7 (PRE_REGISTRATION.md:66: "nothing: it is an
     imported upper bound, not a prediction").

FRAMING (non-negotiable, per ground rules): M0 predicts Gaussian statistics
(no topological defects). A null result (data consistent with the Gaussian
ensemble) is consistent with M0 and with plain LCDM EQUALLY -- it is not
evidence for K3xT2, and changes the free-parameter count by exactly zero.

Run:
  /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    audit/reverse_zero/E5-cmb-tda/cmb_tda.py --n-sims 10 --out smoke_test.json
  (smoke test first, per advisor guidance, before the N_SIMS=100 run)
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
    """Pseudo-Cl / fsky approximation from the data's own masked map."""
    if lmax is None:
        lmax = 2 * nside
    fsky = mask.mean()
    m0 = np.where(mask > 0, masked_map, hp.UNSEEN)
    m_for_anafast = np.where(mask > 0, masked_map - masked_map[mask > 0].mean(), 0.0)
    cl = hp.anafast(m_for_anafast, lmax=lmax)
    cl_corrected = cl / max(fsky, 1e-6)
    return cl_corrected, fsky


def build_complex_and_betti_curve(temp, mask, nside, nu_grid, sublevel=True):
    """
    Build the unmasked-pixel adjacency simplicial complex (vertices, edges,
    triangles among mutual neighbours), run a lower-star filtration by
    T/sigma (or -T/sigma for superlevel), and return Betti-0/1 curves on
    nu_grid plus the Euler characteristic curve, using GUDHI's persistence.
    """
    npix = hp.nside2npix(nside)
    unmasked = np.where(mask > 0)[0]
    sigma = temp[unmasked].std()
    t = temp / sigma
    if not sublevel:
        t = -t

    idx_map = -np.ones(npix, dtype=np.int64)
    idx_map[unmasked] = np.arange(unmasked.size)

    neighbours = hp.get_all_neighbours(nside, unmasked)  # shape (8, N)
    st = gudhi.SimplexTree()
    for i, p in enumerate(unmasked):
        st.insert([i], filtration=t[p])

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
            e = (min(a, b), max(a, b))
            if e in edge_set:
                continue
            edge_set.add(e)
            f = max(t[unmasked[a]], t[unmasked[b]])
            st.insert([int(a), int(b)], filtration=float(f))

    # triangles: 3-cliques among mutual neighbours already in edge_set
    adj = {}
    for a, b in edge_set:
        adj.setdefault(a, set()).add(b)
        adj.setdefault(b, set()).add(a)
    n_tri = 0
    for a, b in edge_set:
        common = adj.get(a, set()) & adj.get(b, set())
        for c in common:
            if c > b:
                tri_pts = [unmasked[a], unmasked[b], unmasked[c]]
                f = max(t[p] for p in tri_pts)
                st.insert([int(a), int(b), int(c)], filtration=float(f))
                n_tri += 1

    st.make_filtration_non_decreasing()
    st.compute_persistence(persistence_dim_max=True)
    diag0 = st.persistence_intervals_in_dimension(0)
    diag1 = st.persistence_intervals_in_dimension(1)

    b0 = np.array([np.sum((diag0[:, 0] <= nu) & ((diag0[:, 1] > nu) | np.isinf(diag0[:, 1]))) for nu in nu_grid])
    b1 = np.array([np.sum((diag1[:, 0] <= nu) & ((diag1[:, 1] > nu) | np.isinf(diag1[:, 1]))) for nu in nu_grid]) if diag1.size else np.zeros_like(nu_grid, dtype=int)
    chi = b0 - b1
    return b0, b1, chi, int(len(unmasked)), int(len(edge_set)), n_tri


def coarsen(curve, n_bins):
    n = len(curve)
    edges = np.linspace(0, n, n_bins + 1).astype(int)
    return np.array([curve[edges[i]:edges[i + 1]].mean() for i in range(n_bins)])


def chi2_stat(vec, mean, cov_inv):
    d = vec - mean
    return float(d @ cov_inv @ d)


def inject_ks_strings(temp_map, nside, mask, n_segments, gmu, v_gamma=0.6, rng=None):
    """
    Toy Kaiser-Stebbins injection. A cosmic string segment is modeled as a
    great-circle arc on the sphere between two random points within an
    angular separation drawn uniform in [2, 20] deg (a coherence-length
    proxy, NOT a scaling-solution network). Each pixel whose position is on
    one side of the arc's dividing great circle (within an angular band of
    +/-0.5 deg of the arc, to keep the step localized to the string) receives
    a temperature step Delta T/T = 8 pi G mu * v_gamma. This is a documented
    toy model (no network evolution, no correlated segment population, no
    scaling density of strings): it is used only to calibrate the TDA
    statistic's SENSITIVITY to a KS-like step pattern, not to forecast a
    physical string signal.
    """
    if rng is None:
        rng = np.random.default_rng(0)
    npix = hp.nside2npix(nside)
    injected = temp_map.copy()
    T_CMB = 2.7255  # K, standard CMB monopole temperature used to convert mK map + dimensionless step
    step_mK = 8 * np.pi * gmu * v_gamma * T_CMB * 1000.0  # mK, using dT/T * T_CMB in mK
    unmasked = np.where(mask > 0)[0]
    theta, phi = hp.pix2ang(nside, unmasked)
    vecs = hp.pix2vec(nside, unmasked)
    vecs = np.array(vecs).T  # (N,3)
    for _ in range(n_segments):
        # random great circle: normal vector n
        n = rng.normal(size=3)
        n /= np.linalg.norm(n)
        # random arc extent (angular half-length) -- restrict step to a band
        # near the great circle to localize the "segment"
        band = np.radians(rng.uniform(2, 20))
        side = vecs @ n
        near_band = np.abs(side) < np.sin(band) * 0.15  # localizes to a strip
        step_sign = np.sign(side)
        injected[unmasked[near_band]] += step_mK * step_sign[near_band]
    return injected, step_mK


def run_pipeline(n_sims, n_segments_inj, gmu_list, out_path, seed_base=BASE_SEED, lmax_override=None):
    log("loading data map + mask")
    data_map, mask, mask_frac = load_data_map_and_mask()
    fsky = float(mask.mean())
    log(f"ud_graded to nside={NSIDE_WORK}, unmasked pixels={int(mask.sum())}/{mask.size} (fsky={fsky:.4f})")

    cl, fsky_cl = estimate_cl(data_map, mask, NSIDE_WORK, lmax=lmax_override)
    log(f"estimated Cl via anafast/fsky on masked map, lmax={len(cl)-1}, fsky_used={fsky_cl:.4f}")

    nu_grid = np.linspace(NU_MIN, NU_MAX, NU_STEP_GRID)

    log("computing data Betti/Euler curves (sublevel)")
    b0_d, b1_d, chi_d, n_unmasked, n_edges, n_tri = build_complex_and_betti_curve(
        data_map, mask, NSIDE_WORK, nu_grid, sublevel=True)
    log("computing data Betti/Euler curves (superlevel, consistency check)")
    b0_d_super, b1_d_super, chi_d_super, *_ = build_complex_and_betti_curve(
        data_map, mask, NSIDE_WORK, nu_grid, sublevel=False)

    data_sigma_unmasked = float(data_map[mask > 0].std())

    log(f"generating {n_sims} Gaussian null sims (synfast, seeds {seed_base}..{seed_base+n_sims-1})")
    sim_b0, sim_b1, sim_chi = [], [], []
    sim_var = []
    for i in range(n_sims):
        rng_seed = seed_base + i
        np.random.seed(rng_seed)
        sim = hp.synfast(cl, nside=NSIDE_WORK, new=True, verbose=False)
        sim_var.append(float(sim[mask > 0].var()))
        b0, b1, chi, *_ = build_complex_and_betti_curve(sim, mask, NSIDE_WORK, nu_grid, sublevel=True)
        sim_b0.append(b0); sim_b1.append(b1); sim_chi.append(chi)
        if (i + 1) % max(1, n_sims // 10) == 0:
            log(f"  sim {i+1}/{n_sims} done")
    sim_b0 = np.array(sim_b0); sim_b1 = np.array(sim_b1); sim_chi = np.array(sim_chi)

    data_var = float(data_map[mask > 0].var())
    sim_var_mean = float(np.mean(sim_var))
    var_ratio = sim_var_mean / data_var if data_var else float("nan")
    log(f"CALIBRATION CHECK: data unmasked-pixel variance={data_var:.6g}, "
        f"sim ensemble mean variance={sim_var_mean:.6g}, ratio={var_ratio:.4f}")

    def coarse_stats(sim_curves, data_curve):
        sim_coarse = np.array([coarsen(c, N_BINS) for c in sim_curves])
        data_coarse = coarsen(data_curve, N_BINS)
        mean = sim_coarse.mean(axis=0)
        cov = np.cov(sim_coarse, rowvar=False) + 1e-8 * np.eye(N_BINS)
        n = sim_coarse.shape[0]
        p = N_BINS
        hartlap = (n - p - 2) / (n - 1) if n > p + 2 else float("nan")
        cov_inv = np.linalg.pinv(cov) * (hartlap if not np.isnan(hartlap) else 1.0)
        chi2 = chi2_stat(data_coarse, mean, cov_inv)
        from scipy.stats import chi2 as chi2dist
        pval = float(chi2dist.sf(chi2, df=p)) if not np.isnan(hartlap) else None
        # empirical rank p: leave-one-out chi2 for each sim against the rest
        loo_chi2 = []
        for j in range(n):
            others = np.delete(sim_coarse, j, axis=0)
            m2 = others.mean(axis=0)
            c2 = np.cov(others, rowvar=False) + 1e-8 * np.eye(N_BINS)
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
        }

    stats_b0 = coarse_stats(sim_b0, b0_d)
    stats_b1 = coarse_stats(sim_b1, b1_d)
    stats_chi = coarse_stats(sim_chi, chi_d)

    result = {
        "tier": "X (exploratory numerics)",
        "framing": "M0 predicts Gaussian statistics (no topological defects). A null "
                    "result here is consistent with M0 and plain LCDM equally; it is "
                    "NOT evidence for K3xT2 and changes the free-parameter count by 0.",
        "data_map": str(MAP_PATH.relative_to(REPO)),
        "mask": str(MASK_PATH.relative_to(REPO)),
        "nside_native": 512,
        "nside_work": NSIDE_WORK,
        "smoothing": "no additional smoothing; WMAP 9yr ILC native effective beam (foreground-cleaned combination map)",
        "mask_threshold_after_udgrade": 0.5,
        "fsky": fsky,
        "n_unmasked_pixels_work_res": n_unmasked,
        "n_edges": n_edges,
        "n_triangles": n_tri,
        "cl_source": "hp.anafast on the SAME masked ILC map at nside_work, divided by fsky (pseudo-Cl/fsky approximation, NOT CAMB, NOT MASTER-deconvolved -- stated explicitly)",
        "lmax_used": len(cl) - 1,
        "nu_grid_sigma": nu_grid.tolist(),
        "n_sims": n_sims,
        "seeds": list(range(seed_base, seed_base + n_sims)),
        "calibration_check": {
            "data_unmasked_pixel_variance": data_var,
            "sim_ensemble_mean_variance": sim_var_mean,
            "ratio_sim_over_data": var_ratio,
            "note": "ratio should be close to 1; if not, the anafast/fsky Cl "
                    "estimate does not reproduce the ILC map's own variance and "
                    "the null is not directly interpretable as p-value-grade"
        },
        "data_betti_curve": {
            "sublevel": {"b0": b0_d.tolist(), "b1": b1_d.tolist(), "euler_chi": chi_d.tolist()},
            "superlevel_consistency_check": {"b0": b0_d_super.tolist(), "b1": b1_d_super.tolist(), "euler_chi": chi_d_super.tolist()},
        },
        "statistic_b0": stats_b0,
        "statistic_b1": stats_b1,
        "statistic_euler_chi": stats_chi,
        "n_coarse_bins": N_BINS,
    }

    if gmu_list:
        log(f"running cosmic-string sensitivity scan: Gmu = {gmu_list}, n_segments={n_segments_inj}")
        sens = {}
        # use the fitted null (sim_chi coarse stats) as the reference ensemble's 95th percentile
        sim_chi_coarse = np.array([coarsen(c, N_BINS) for c in sim_chi])
        mean_chi = sim_chi_coarse.mean(axis=0)
        cov_chi = np.cov(sim_chi_coarse, rowvar=False) + 1e-8 * np.eye(N_BINS)
        n_null = sim_chi_coarse.shape[0]
        hartlap_null = (n_null - N_BINS - 2) / (n_null - 1)
        cov_inv_chi = np.linalg.pinv(cov_chi) * hartlap_null
        null_chi2_values = np.array([chi2_stat(v, mean_chi, cov_inv_chi) for v in sim_chi_coarse])
        thresh_95 = float(np.percentile(null_chi2_values, 95))

        n_inj_sims = min(20, n_sims)  # smaller ensemble per Gmu point (budget)
        for gmu in gmu_list:
            rng = np.random.default_rng(1000 + int(round(np.log10(gmu) * -1000)))
            inj_chi2 = []
            for i in range(n_inj_sims):
                np.random.seed(seed_base + 5000 + i)
                sim = hp.synfast(cl, nside=NSIDE_WORK, new=True, verbose=False)
                sim_inj, step_mK = inject_ks_strings(sim, NSIDE_WORK, mask, n_segments_inj, gmu, rng=rng)
                _, _, chi_i, *_ = build_complex_and_betti_curve(sim_inj, mask, NSIDE_WORK, nu_grid, sublevel=True)
                chi_coarse_i = coarsen(chi_i, N_BINS)
                inj_chi2.append(chi2_stat(chi_coarse_i, mean_chi, cov_inv_chi))
            inj_chi2 = np.array(inj_chi2)
            frac_separated = float(np.mean(inj_chi2 > thresh_95))
            sens[str(gmu)] = {
                "n_injected_sims": n_inj_sims, "n_segments": n_segments_inj,
                "step_mK_example": step_mK,
                "fraction_exceeding_null_95th_pct": frac_separated,
                "null_95th_pct_chi2": thresh_95,
                "mean_injected_chi2": float(inj_chi2.mean()),
            }
            log(f"  Gmu={gmu}: frac exceeding null 95th pct = {frac_separated:.2f}")
        result["cosmic_string_sensitivity"] = {
            "model": "toy Kaiser-Stebbins: straight great-circle segments, "
                     "Delta T/T = 8 pi G mu v_gamma per crossing, v_gamma=0.6 fixed, "
                     "NOT a scaling-solution network -- sensitivity of the TDA "
                     "statistic to this injection only, not a cosmological bound",
            "statistic_used": "euler_chi coarse chi2 vs null 95th percentile",
            "planck_gmu_bound_comparison": "LeanMaster planckGmuBound = 1.5e-7 "
                "(PRE_REGISTRATION.md:66: 'nothing: it is an imported upper bound, "
                "not a prediction'). Any Gmu found here is a TDA-pipeline "
                "sensitivity threshold for THIS toy injection, compared in order "
                "of magnitude only, never claimed equal or as a new bound.",
            "gmu_scan": sens,
        }

    with open(out_path, "w") as f:
        json.dump(result, f, indent=2)
    log(f"wrote {out_path}")
    return result


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-sims", type=int, default=100)
    ap.add_argument("--n-segments-inj", type=int, default=30)
    ap.add_argument("--gmu", type=str, default="1e-6,1e-7,1e-8,1e-9")
    ap.add_argument("--no-sensitivity", action="store_true")
    ap.add_argument("--out", type=str, default="e5_cmb_tda_report.json")
    ap.add_argument("--lmax", type=int, default=None)
    args = ap.parse_args()
    gmu_list = [] if args.no_sensitivity else [float(x) for x in args.gmu.split(",")]
    out_path = HERE / args.out
    run_pipeline(args.n_sims, args.n_segments_inj, gmu_list, out_path, lmax_override=args.lmax)
