#!/usr/bin/env python3
"""
E5b -- fills the two gaps the E5 resume note named as missing (ground
rules' nulls (i) Poisson-in-volume and (iii) >=20 CAMB lognormal mocks
with bias fitted ONLY to xi(r)), against the SAME N=25000 SDSS DR17
subsample E5 (cosmic_web_tda_scaled.py) already built and tested. It does
NOT redo E5's z-shuffle null or known-answer controls (both already run
and, after this session's top_bars() r_trunc fix, both pass -- see E5's
report and its "gate" field on planted_voids_H2).

WHY A GENUINE POISSON NULL, SEPARATE FROM E5'S Z-SHUFFLE NULL: E5's
docstring claims the z-permutation null "stands in for" a Poisson-in-
volume null. That claim is only half right. Z-permutation keeps the
OBSERVED (ra,dec) pairs exactly -- so it preserves the data's full
ANGULAR two-point clustering and only randomises the radial direction.
A true Poisson-in-volume null must randomise angle too (within the same
footprint), which is a different, strictly weaker null: any TDA signal
found against the z-shuffle null could be pure angular clustering
(gravity acting on the sky, not a 3D void/filament signal), while a
mismatch against BOTH nulls is closer to a genuine 3D-structure claim.
This script builds that missing angular-random null explicitly, from a
HEALPix (nside=64) occupied-pixel mask of the DATA's own (ra,dec)
footprint -- not the ra/dec bounding rectangle, which the docstring of
E5 and this project's round-1 retraction already flagged as the wrong
angular support to sample from.

THE ACTUAL M0 TEST (per the ground rules, this is the decisive one, not
E5's z-shuffle-envelope numbers): >=20 CAMB lognormal mock catalogues,
built as
    P_target(k) = b^2 * P_lin(k)      [CAMB linear P(k) at z_eff, M0 cosmology]
    xi_target(r) = pk2xi(P_target)
    xi_G(r)      = ln(1 + xi_target(r))          <-- REQUIRED step, see below
    P_G(k)       = xi2pk(xi_G)
    delta_G(x)   ~ GaussianRandomField(P_G)      [on a grid covering the
                                                    survey wedge's bounding
                                                    box, NOT a symmetric cube]
    delta_LN(x)  = exp(delta_G(x) - sigma_G^2/2) - 1
    galaxies     ~ Poisson( n_bar * (1+delta_LN) ) per cell, masked to the
                   SAME occupied-HEALPix-footprint x SAME comoving-shell
                   window as the real E5 subsample, then pushed through
                   the SAME select_flat_density_subsample() call E5 used
                   on the data, to N_FINAL=25000.
The ln(1+xi) step is not optional bookkeeping: it is what makes the
EXPONENTIATED field's own two-point function equal xi_target (a Gaussian
field's exponential has a two-point function that is NOT its own P(k)
transform; matching P_G to P_target directly and then exponentiating
would give the wrong amplitude and shape -- a known trap, stated here
because skipping it silently is the single most common lognormal-mock
bug).

BIAS: b^2 is fit by ordinary least squares to the DATA's own Landy-
Szalay xi(r) over r in [20,60] Mpc/h (r < 20 is the nonlinear/fibre-
collision regime, r > 60 Mpc/h is noise-dominated at this survey volume
-- BOTH bounds fixed here, in the source, BEFORE the fit is run). b is
NEVER fit to any topological (Betti/persistence) statistic -- stated
explicitly per the ground rules' requirement.

VALIDATION GATE (must pass before any mock Betti curve is trusted):
the mock galaxies' OWN Landy-Szalay xi(r), averaged over the mocks, is
compared to the DATA's xi(r) over the same [20,60] Mpc/h fit range. A
gate failure (ratio far from 1) means the lognormal pipeline's
normalisation is broken and the Betti-curve comparison below is not
diagnostic of anything -- this is reported as a hard pass/fail, not
folded into the interpretation prose.

FRAMING (identical to E5, repeated per ground rules): M0 (frozen LCDM
background, Omega_m=0.31115, screening/symmetron sector deleted) is a
HYPOTHESIS CHANGE, not a derivation from K3xT2. A lognormal-mock match
means "SDSS DR17's galaxy topology at this N and scale range is
consistent with M0's linear P(k) plus a linearly biased lognormal field"
-- it is NOT a confirmation of K3xT2, and it says nothing about
mu_sym, c4_pta_product, or any octad/moonshine construction (none of
those touch a galaxy-clustering observable; see E4's report and
REPORT.md). A mismatch at r below a few Mpc/h is expected and NOT
diagnostic: lognormal transforms reproduce the 1-point PDF and the
input 2-point function but not higher-order correlations, filament
shapes, or deep-void statistics, and SDSS fibre-collision incompleteness
suppresses close pairs at small separation (not corrected for here,
stated rather than silently assumed away, same as E5).

Tier: X (exploratory numerics) throughout. Not a K3xT2 (LeanMaster
Stream 8) prediction test -- no constructed octad-indexed observable
exists for a galaxy catalogue (E4's report, REPORT.md).

MEMORY SAFETY: this script never builds more than ONE alpha complex, ONE
grid field, and ONE random catalogue at a time; each large array is
`del`-ed and gc.collect()-ed before the next mock. The grid covers only
the survey wedge's own (padded) bounding box (~94x169x109 cells at a
2.0 Mpc/h cell -- NOT a symmetric 256^3 cube, which would waste >75% of
its cells outside the wedge for no benefit) -- reported below.

Command:
  timeout 1700 prlimit --as=10737418200 -- \
    /home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python \
    audit/reverse_zero/E5-cosmic-web-tda-scaled/e5b_lognormal_and_poisson_null.py
"""
import gc
import json
import os
import sys
import time

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import gudhi
import healpy as hp
import camb
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import cosmic_web_tda_scaled as base  # noqa: E402  (reuse E5's exact selection machinery)

# ---------------------------------------------------------------------------
# Shared constants (imported from E5 so the two scripts cannot silently drift)
# ---------------------------------------------------------------------------
OMEGA_M = base.OMEGA_M
N_SHELLS = base.N_SHELLS_FOR_SELECTION
SELECTION_SHELL_IDX = base.SELECTION_SHELL_IDX
N_TARGET_PER_SHELL = base.N_TARGET_PER_SHELL
N_FINAL = base.N_FINAL
R_MAX_PERS = base.R_MAX_PERS
MAX_ALPHA_SQ = base.MAX_ALPHA_SQ
R_GRID = base.R_GRID
STAT_R_RANGE = base.STAT_R_RANGE

# ---- this script's own pre-registered constants (fixed before first run) ----
NSIDE_MASK = 64                    # healpix resolution for the angular footprint mask
N_POISSON_NULL_SEEDS = 20          # ground rules ">= 20 mocks" applied to the Poisson null too
N_LOGNORMAL_MOCKS = 22             # >= 20 (ground rules' floor), small margin for any residual failures
XI_FIT_RANGE_MPC_H = (20.0, 60.0)  # PRE-STATED bias-fit range: below is nonlinear/fibre-collision,
                                    # above is noise-dominated at this survey volume
XI_RANDOM_OVERSAMPLE = 4           # random-catalogue size = this x N_FINAL, for the LS estimator
                                    # (reduced from an initial 8x after this run's own timing test showed
                                    # RR pair-counting cost scales as N_rand^2 and dominates runtime; 4x is
                                    # still within the standard 3-10x LS practice, see report's runtime note)
GRID_CELL_MPC_H = 2.0              # lognormal field grid cell size
GRID_PAD_MPC_H = 25.0              # padding beyond the data's own bounding box, each side
Z_EFF = 0.08                       # effective redshift of the selected shell window (see report)
PLANCK18_PARAMS = dict(H0=67.66, ombh2=0.02237, omch2=0.1200, mnu=0.06, omk=0.0, tau=0.0544,
                        As_raw=2.100e-9, ns=0.9649, sigma8_target=0.8102)
GATE_TOLERANCE = 0.5               # informational only now (raw ratio still reported) -- see GATE_Z_THRESHOLD for the actual pass/fail
N_GATE_MOCKS = N_LOGNORMAL_MOCKS   # gate EVERY mock's own xi(r) (RR is shared/precomputed, so this is cheap
                                    # per mock; using all of them, not a subset of 5-8, gives a much tighter
                                    # mock-to-mock standard-error estimate for the z-score gate below)
GATE_Z_THRESHOLD = 3.0             # PASS iff RMS over the fit range of (xi_mock_mean-xi_data)/sigma_total < this,
                                    # sigma_total = data's OWN Poisson shot-noise error (+) mock-to-mock standard
                                    # error of the mean, combined in quadrature. A raw ratio gate was tried first
                                    # and rejected: xi(r) crosses zero within the pre-registered [20,60] Mpc/h fit
                                    # range at this survey volume (mean galaxy spacing ~5 Mpc/h, N=25000), so a
                                    # ratio blows up near every zero-crossing regardless of whether the mocks
                                    # actually match the data -- not a meaningful pass/fail criterion. The z-score
                                    # version stays well-behaved there because sigma_total also grows near the
                                    # noise floor.
BOOT_SEED_BASE = 5000


def log(msg):
    print("[e5b %6.1fs] %s" % (time.time() - T0, msg), flush=True)


# ---------------------------------------------------------------------------
# 1. Load the SAME real subsample E5 used (re-derive with E5's own functions,
#    not re-read from E5's json, so a confound cannot enter through re-typing
#    numbers by hand).
# ---------------------------------------------------------------------------
def load_real_selection():
    df = pd.read_csv(base.REAL_CSV, comment="#")
    valid = np.isfinite(df["z"]) & (df["z"] > 0)
    df = df[valid].reset_index(drop=True)
    z_max_grid = float(df["z"].max())
    r_all = base.comoving_r_mpc_over_h(df["z"].to_numpy(), z_max_grid)
    r_shell_max = float(r_all.max())
    edges = base.equal_volume_shell_edges(r_shell_max, N_SHELLS)
    real_idx, real_shell_counts = base.select_flat_density_subsample(
        r_all, edges, SELECTION_SHELL_IDX, N_TARGET_PER_SHELL, base.SEED_SUBSAMPLE_REAL)
    ra = df["ra"].to_numpy(); dec = df["dec"].to_numpy()
    real_ra, real_dec, real_r = ra[real_idx], dec[real_idx], r_all[real_idx]
    real_xyz = base.radec_r_to_xyz(real_ra, real_dec, real_r)
    r_lo, r_hi = float(edges[SELECTION_SHELL_IDX[0]]), float(edges[SELECTION_SHELL_IDX[-1] + 1])
    return {
        "df_ra": ra, "df_dec": dec, "df_r": r_all, "edges": edges,
        "real_xyz": real_xyz, "real_ra": real_ra, "real_dec": real_dec, "real_r": real_r,
        "r_lo": r_lo, "r_hi": r_hi, "real_shell_counts": real_shell_counts,
    }


# ---------------------------------------------------------------------------
# 2. HEALPix occupied-footprint mask + angular/radial samplers
# ---------------------------------------------------------------------------
def build_mask(ra_deg, dec_deg, nside):
    theta = np.radians(90.0 - dec_deg)
    phi = np.radians(ra_deg)
    pix = hp.ang2pix(nside, theta, phi)
    occupied = np.unique(pix)
    occ_set = np.zeros(hp.nside2npix(nside), dtype=bool)
    occ_set[occupied] = True
    frac = occupied.size / hp.nside2npix(nside)
    return occ_set, frac


def sample_mask_angles(occ_set, nside, ra_min, ra_max, dec_min, dec_max, n, rng):
    """Rejection-sample (ra,dec) uniformly on the sphere, restricted to the
    occupied-HEALPix-pixel set. Isotropic sampling: ra ~ U(ra_min,ra_max),
    sin(dec) ~ U(sin(dec_min), sin(dec_max))."""
    out_ra = np.empty(n); out_dec = np.empty(n)
    filled = 0
    sdec_min, sdec_max = np.sin(np.radians(dec_min)), np.sin(np.radians(dec_max))
    while filled < n:
        batch = max(int((n - filled) * 1.4) + 1000, 2000)
        ra_try = rng.uniform(ra_min, ra_max, batch)
        dec_try = np.degrees(np.arcsin(rng.uniform(sdec_min, sdec_max, batch)))
        theta = np.radians(90.0 - dec_try); phi = np.radians(ra_try)
        pix = hp.ang2pix(nside, theta, phi)
        keep = occ_set[pix]
        k = min(keep.sum(), n - filled)
        out_ra[filled:filled + k] = ra_try[keep][:k]
        out_dec[filled:filled + k] = dec_try[keep][:k]
        filled += k
    return out_ra, out_dec


def sample_uniform_r_shell(r_lo, r_hi, n, rng):
    u = rng.uniform(0.0, 1.0, n)
    return (u * (r_hi ** 3 - r_lo ** 3) + r_lo ** 3) ** (1.0 / 3.0)


def poisson_footprint_catalogue(occ_set, nside, ra_min, ra_max, dec_min, dec_max, r_lo, r_hi, n, seed):
    rng = np.random.RandomState(seed)
    ra, dec = sample_mask_angles(occ_set, nside, ra_min, ra_max, dec_min, dec_max, n, rng)
    r = sample_uniform_r_shell(r_lo, r_hi, n, rng)
    return ra, dec, r


# ---------------------------------------------------------------------------
# 3. Landy-Szalay xi(r) via cKDTree.count_neighbors (exact pair counts)
# ---------------------------------------------------------------------------
def precompute_random_tree(rand_xyz, r_edges, verbose=False):
    """RR pair counts depend ONLY on the random catalogue, never on the data
    or any mock -- computed ONCE here and reused by every ls_xi_with_tR()
    call below. Skipping this and rebuilding RR per-call (the initial, naive
    version of this script) made RR's O(nR^2)-pair-count cost repeat once
    per mock in the xi(r) gate, which is what made the first timing test
    look 'stuck': RR at nR=200000 took ~7 minutes, times 6 calls (data + 5
    mocks) would have been ~42 minutes, over the ground rules' per-script
    budget. Reusing one RR (and dropping the oversample from 8x to 4x,
    which cuts RR's cost ~4x again since it scales as nR^2) is the fix."""
    nR = rand_xyz.shape[0]
    t0 = time.time()
    tR = cKDTree(rand_xyz)
    RR_cum = np.array(tR.count_neighbors(tR, r_edges), dtype=float)
    RR = np.diff(RR_cum) / 2.0
    nRR = nR * (nR - 1) / 2.0
    if verbose:
        log("  precompute_random_tree: nR=%d, RR done (%.2fs)" % (nR, time.time() - t0))
    return {"tR": tR, "RR": RR, "nRR": nRR, "nR": nR, "r_edges": r_edges}


def ls_xi_with_tR(data_xyz, rtree, verbose=False):
    """LS xi(r) using a PRECOMPUTED random tree + RR counts (see
    precompute_random_tree). Only DD (data self-pairs) and DR (data-random
    cross-pairs) are computed per call -- both far cheaper than RR since
    nD << nR is NOT assumed; DD scales as nD^2 (fixed cost, nD=N_FINAL for
    every catalogue compared here) and DR as nD*nR (linear in nR)."""
    nD = data_xyz.shape[0]
    r_edges = rtree["r_edges"]
    t0 = time.time()
    tD = cKDTree(data_xyz)
    DD_cum = np.array(tD.count_neighbors(tD, r_edges), dtype=float)
    DD = np.diff(DD_cum) / 2.0
    if verbose:
        log("  ls_xi_with_tR: DD done (%.2fs), counting DR..." % (time.time() - t0))
    t0 = time.time()
    DR_cum = np.array(tD.count_neighbors(rtree["tR"], r_edges), dtype=float)
    DR = np.diff(DR_cum)
    if verbose:
        log("  ls_xi_with_tR: DR done (%.2fs)" % (time.time() - t0))
    nDD = nD * (nD - 1) / 2.0
    nDR = float(nD) * float(rtree["nR"])
    dd = DD / nDD; rr = rtree["RR"] / rtree["nRR"]; dr = DR / nDR
    with np.errstate(divide="ignore", invalid="ignore"):
        xi = np.where(rr > 0, (dd - 2 * dr + rr) / rr, np.nan)
    return xi, DD, rtree["RR"], DR


# ---------------------------------------------------------------------------
# 4. CAMB linear P(k) at z_eff, Hankel-type transforms P(k)<->xi(r)
# ---------------------------------------------------------------------------
def camb_linear_pk(z_eff):
    p = PLANCK18_PARAMS
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=p["H0"], ombh2=p["ombh2"], omch2=p["omch2"], mnu=p["mnu"], omk=p["omk"], tau=p["tau"])
    pars.InitPower.set_params(As=p["As_raw"], ns=p["ns"])
    pars.set_matter_power(redshifts=[z_eff, 0.0], kmax=15.0)
    pars.NonLinear = camb.model.NonLinear_none
    results = camb.get_results(pars)
    kh, zs, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=12.0, npoints=600)
    sigma8_raw = results.get_sigma8()
    omega_m_arith = (p["ombh2"] + p["omch2"]) / (p["H0"] / 100.0) ** 2
    omega_m_pars = pars.omegam
    # BUG FOUND AND FIXED (2026-09-19, caught by an advisor-review sigma8 self-check landing at ratio
    # 1.043 instead of ~1.000): get_matter_power_spectrum()'s returned `zs` is SORTED (increasing), so
    # `pk` is indexed by sorted-z position -- but CAMBdata.get_sigma8() returns values in
    # pars.Transfer.PK_redshifts's ORIGINAL INPUT order, which is UNSORTED here ([z_eff, 0.0] =
    # [0.08, 0.0], already decreasing). Indexing sigma8_raw with idx_z0 found from the SORTED `zs`
    # array silently pulled sigma8(z_eff=0.08)=0.778 instead of sigma8(z=0)=0.812 -- an 8.8% error in
    # sigma8_z0, and hence in `rescale` and every downstream P(k) amplitude (xi_lin, b_fit, sigma_G^2).
    # Confirmed by comparing a single-redshift CAMB call (unambiguous, ratio 1.0003) against this
    # multi-redshift call: verified sigma8_raw is ordered by pars.Transfer.PK_redshifts, NOT by `zs`.
    input_z_order = list(pars.Transfer.PK_redshifts)[:len(sigma8_raw)]
    idx_zeff = int(np.argmin(np.abs(np.array(zs) - z_eff)))       # indexes `pk` (sorted-z order)
    idx_z0 = int(np.argmin(np.abs(np.array(zs) - 0.0)))           # indexes `pk` (sorted-z order)
    idx_z0_sigma8 = int(np.argmin(np.abs(np.array(input_z_order) - 0.0)))  # indexes sigma8_raw (input order)
    pk_zeff = pk[idx_zeff]
    pk_z0 = pk[idx_z0]
    sigma8_z0 = float(sigma8_raw[idx_z0_sigma8])
    rescale = (p["sigma8_target"] / sigma8_z0) ** 2
    pk_zeff_rescaled = pk_zeff * rescale
    # SELF-CHECK (added after an advisor review flagged that pk2xi/xi2pk's round-trip test, done elsewhere
    # in this script, validates SHAPE but not ABSOLUTE normalisation -- a constant multiplicative error would
    # cancel in a round trip). Compute sigma8 directly from kh/pk_z0 with a top-hat window at R=8 Mpc/h, using
    # the SAME quadrature style as pk2xi/xi2pk, and compare to CAMB's own reported sigma8(z=0). Run once
    # (2026-09-19): ratio 1.0003 (0.03% agreement) -- rules out a normalisation bug in the shared transform.
    def _tophat_w(x):
        return np.where(np.abs(x) < 1e-4, 1.0, 3 * (np.sin(x) - x * np.cos(x)) / x ** 3)
    sigma8_selfcheck = float(np.sqrt(np.trapezoid(kh ** 2 * pk_z0 * _tophat_w(kh * 8.0) ** 2, kh) / (2 * np.pi ** 2)))
    return {
        "kh": kh, "pk_zeff": pk_zeff_rescaled, "z_eff_actual": float(zs[idx_zeff]),
        "sigma8_z0_raw": sigma8_z0, "sigma8_z0_rescale_factor": float(rescale),
        "sigma8_z0_selfcheck_from_own_quadrature": sigma8_selfcheck,
        "sigma8_selfcheck_ratio_to_camb": sigma8_selfcheck / sigma8_z0,
        "omega_m_arith_check": float(omega_m_arith), "omega_m_pars": float(omega_m_pars),
        "planck18_params_imported": {k: v for k, v in p.items()},
    }


def pk2xi(k, pk, r):
    r = np.atleast_1d(r).astype(float)
    xi = np.empty_like(r)
    for i, rr in enumerate(r):
        integrand = k ** 2 * pk * np.sinc(k * rr / np.pi)
        xi[i] = np.trapezoid(integrand, k) / (2 * np.pi ** 2)
    return xi


def xi2pk(r, xi, k):
    k = np.atleast_1d(k).astype(float)
    pk = np.empty_like(k)
    for i, kk in enumerate(k):
        integrand = r ** 2 * xi * np.sinc(r * kk / np.pi)
        pk[i] = 4 * np.pi * np.trapezoid(integrand, r)
    return pk


# ---------------------------------------------------------------------------
# 5. Gaussian random field on the wedge's bounding box (validated normalisation,
#    see this run's scratch calibration: measured P(k) matches input to <6%)
# ---------------------------------------------------------------------------
def gaussian_random_field(k_of_pg, pg_of_k, shape, cell, seed):
    rng = np.random.RandomState(seed)
    noise = rng.standard_normal(shape)
    kx = 2 * np.pi * np.fft.fftfreq(shape[0], d=cell)
    ky = 2 * np.pi * np.fft.fftfreq(shape[1], d=cell)
    kz = 2 * np.pi * np.fft.rfftfreq(shape[2], d=cell)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    kmag = np.sqrt(KX ** 2 + KY ** 2 + KZ ** 2)
    del KX, KY, KZ
    pk_grid = np.interp(kmag, k_of_pg, pg_of_k, left=pg_of_k[0], right=0.0)
    pk_grid[kmag == 0] = 0.0
    field_k = np.fft.rfftn(noise) * np.sqrt(pk_grid / cell ** 3)
    del noise, pk_grid, kmag
    field_x = np.fft.irfftn(field_k, s=shape, axes=(0, 1, 2))
    del field_k
    return field_x


def lognormal_sample_points(field_g, sigma_g2, n_bar, box_origin, cell, rng):
    """Poisson-sample points from n_bar*(1+delta_LN) per cell; return
    absolute Cartesian coordinates (box_origin + cell-local uniform offset)."""
    delta_ln = np.exp(field_g - 0.5 * sigma_g2) - 1.0
    lam = np.clip(n_bar * cell ** 3 * (1.0 + delta_ln), 0.0, None)
    counts = rng.poisson(lam)
    idx = np.nonzero(counts)
    if idx[0].size == 0:
        return np.empty((0, 3))
    reps = counts[idx]
    ix = np.repeat(idx[0], reps); iy = np.repeat(idx[1], reps); iz = np.repeat(idx[2], reps)
    n_pts = ix.size
    offsets = rng.uniform(0.0, cell, size=(n_pts, 3))
    xyz = np.column_stack([ix, iy, iz]).astype(float) * cell + offsets + np.asarray(box_origin)
    return xyz


def xyz_to_radecr(xyz):
    x, y, z = xyz[:, 0], xyz[:, 1], xyz[:, 2]
    r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
    dec = np.degrees(np.arcsin(np.clip(z / np.maximum(r, 1e-12), -1, 1)))
    ra = np.degrees(np.arctan2(y, x))
    ra = np.mod(ra, 360.0)
    return ra, dec, r


def mask_and_shell_keep(xyz, occ_set, nside, edges):
    """Boolean mask: which Cartesian points fall in the occupied-footprint
    HEALPix mask AND the selected comoving-shell window. Shared by the
    calibration step (below) and mock_to_final_selection, so 'what fraction
    of the box survives' is measured with the EXACT same code that later
    does the real filtering -- no separate analytic formula to drift out
    of sync with it (see calibrate_n_bar's docstring for why an earlier
    analytic solid-angle estimate was wrong by a factor of ~3)."""
    ra, dec, r = xyz_to_radecr(xyz)
    theta = np.radians(90.0 - dec); phi = np.radians(ra)
    pix = hp.ang2pix(nside, theta, phi)
    in_mask = occ_set[pix]
    r_lo, r_hi = float(edges[SELECTION_SHELL_IDX[0]]), float(edges[SELECTION_SHELL_IDX[-1] + 1])
    in_shell = (r >= r_lo) & (r < r_hi)
    return in_mask & in_shell


def calibrate_n_bar(bbox_min, box_size, occ_set, nside, edges, n_target_prewedge, n_probe=300000, seed=8000):
    """Measure (not assume) what fraction of the grid's bounding box
    actually falls inside the mask+shell window, by drawing n_probe points
    UNIFORMLY over the box and running them through the SAME
    mask_and_shell_keep() the mocks are filtered with, then solve for the
    n_bar that gives n_target_prewedge expected survivors.

    WHY THIS REPLACES AN ANALYTIC FORMULA: an earlier version of this
    script computed the wedge volume analytically (solid angle x radial
    shell volume) and derived n_bar from that. Every mock then had only
    ~25000-27000 points survive the mask+shell filter, not the intended
    ~75000 (n_target_prewedge=3xN_FINAL) -- roughly a factor of 3 short.
    The margin was too thin: several mocks then failed the
    >=N_FINAL-available check outright (seeds 4003, 4005 in this run's
    log). Rather than debug the analytic solid-angle integral (a wide,
    non-rectangular ra/dec footprint over a finite radial shell does not
    reduce to a clean closed form once combined with the AXIS-ALIGNED
    CARTESIAN bounding box the grid actually uses), this replaces it with
    a direct Monte Carlo measurement against the exact filter function."""
    rng = np.random.RandomState(seed)
    pts = rng.uniform(0.0, 1.0, size=(n_probe, 3)) * np.asarray(box_size) + np.asarray(bbox_min)
    keep = mask_and_shell_keep(pts, occ_set, nside, edges)
    frac = float(keep.mean())
    box_volume = float(np.prod(box_size))
    n_bar = n_target_prewedge / max(frac * box_volume, 1e-9)
    return n_bar, frac


def mock_to_final_selection(xyz, occ_set, nside, edges, r_lo_full=None):
    """Mask a lognormal-mock Cartesian point cloud to the survey footprint
    and shell window, then run it through E5's OWN select_flat_density_subsample
    to N_FINAL, exactly reproducing the data's selection function."""
    keep = mask_and_shell_keep(xyz, occ_set, nside, edges)
    ra, dec, r = xyz_to_radecr(xyz)
    n_avail = int(keep.sum())
    if n_avail < len(SELECTION_SHELL_IDX) * N_TARGET_PER_SHELL:
        return None, n_avail
    r_kept = r[keep]
    idx_kept = np.nonzero(keep)[0]
    idx_sel, shell_counts = base.select_flat_density_subsample(
        r_kept, edges, SELECTION_SHELL_IDX, N_TARGET_PER_SHELL, seed=int(rng_seed_for_select))
    final_xyz = xyz[idx_kept[idx_sel]]
    return final_xyz, n_avail


rng_seed_for_select = 777  # fixed, stated: subsample-within-selection seed for every mock/poisson realization


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------
T0 = time.time()


def run():
    report = {"generated": "2026-09-19",
              "tool": "gudhi %s, camb %s, healpy %s" % (gudhi.__version__, camb.__version__, hp.__version__),
              "framing_rule": ("M0 (frozen LCDM background, Omega_m=0.31115, screening/symmetron sector "
                                "deleted) is a HYPOTHESIS CHANGE, not a derivation from K3xT2. Agreement "
                                "below means 'consistent with M0's linear P(k) + linear bias + lognormal "
                                "field at the level of this test' -- NOT a K3xT2 confirmation, and it does "
                                "not touch mu_sym, c4_pta_product or any LeanMaster octad/moonshine "
                                "construction (none exists as a galaxy-clustering observable; see E4's "
                                "report and REPORT.md)."),
              "seeds": {"mask_nside": NSIDE_MASK, "poisson_null_seeds": list(range(3000, 3000 + N_POISSON_NULL_SEEDS)),
                        "lognormal_mock_seeds": list(range(4000, 4000 + N_LOGNORMAL_MOCKS)),
                        "ls_random_seed": 6000, "select_within_mock_seed": rng_seed_for_select}}

    log("loading real E5 subsample (re-derived, not re-typed from E5's json)")
    sel = load_real_selection()
    real_xyz = sel["real_xyz"]; edges = sel["edges"]
    r_lo, r_hi = sel["r_lo"], sel["r_hi"]
    ra_min, ra_max = 140.0, 220.0
    dec_min, dec_max = 0.0, 50.0

    log("building HEALPix nside=%d occupied-footprint mask from the data's own (ra,dec)" % NSIDE_MASK)
    occ_set, occ_frac = build_mask(sel["df_ra"], sel["df_dec"], NSIDE_MASK)
    report["mask"] = {"nside": NSIDE_MASK, "npix_total": int(hp.nside2npix(NSIDE_MASK)),
                       "npix_occupied": int(occ_set.sum()), "occupied_fraction_of_full_sky": float(occ_frac),
                       "note": "mask built from the OBSERVED galaxy (ra,dec), not the ra/dec bounding "
                               "rectangle -- avoids this project's earlier round-1 bounding-box artifact."}

    # ================= (i) Poisson-in-footprint null =================
    log("Poisson-in-footprint null: %d realizations, N=%d each" % (N_POISSON_NULL_SEEDS, N_FINAL))
    poisson_betti = {0: [], 1: [], 2: []}
    poisson_euler = []
    poisson_diag = []
    for j, seed in enumerate(range(3000, 3000 + N_POISSON_NULL_SEEDS)):
        ra_p, dec_p, r_p = poisson_footprint_catalogue(occ_set, NSIDE_MASK, ra_min, ra_max, dec_min, dec_max,
                                                         r_lo, r_hi, N_FINAL, seed)
        xyz_p = base.radec_r_to_xyz(ra_p, dec_p, r_p)
        res_p, bd_p = base.alpha_persistence(xyz_p, MAX_ALPHA_SQ, "poisson_seed%d" % seed, HERE)
        e_p, b0_p, b1_p, b2_p = base.euler_curve(bd_p, R_GRID, R_MAX_PERS)
        poisson_betti[0].append(b0_p); poisson_betti[1].append(b1_p); poisson_betti[2].append(b2_p)
        poisson_euler.append(e_p)
        poisson_diag.append({"seed": seed, "betti_numbers_at_truncation": res_p["betti_numbers_at_truncation"],
                              "runtime_sec": res_p["runtime_sec"]})
        del xyz_p, bd_p
        gc.collect()
        if (j + 1) % 5 == 0:
            log("  poisson null %d/%d done" % (j + 1, N_POISSON_NULL_SEEDS))
    for d in (0, 1, 2):
        poisson_betti[d] = np.array(poisson_betti[d])
    poisson_euler = np.array(poisson_euler)

    res_real, bd_real = base.alpha_persistence(real_xyz, MAX_ALPHA_SQ, "real_for_e5b", HERE)
    e_real, b0_real, b1_real, b2_real = base.euler_curve(bd_real, R_GRID, R_MAX_PERS)

    l2_vs_poisson = {}
    outside_poisson = {}
    for d, curve in zip((0, 1, 2), (b0_real, b1_real, b2_real)):
        mean_p = poisson_betti[d].mean(axis=0)
        l2_vs_poisson["H%d" % d] = base.l2_over_range(curve, mean_p, R_GRID, STAT_R_RANGE)
        lo = np.percentile(poisson_betti[d], 2.5, axis=0); hi = np.percentile(poisson_betti[d], 97.5, axis=0)
        mask_r = (R_GRID >= STAT_R_RANGE[0]) & (R_GRID <= STAT_R_RANGE[1])
        out_r = R_GRID[mask_r][(curve[mask_r] < lo[mask_r]) | (curve[mask_r] > hi[mask_r])]
        outside_poisson["H%d" % d] = {"n_outside": int(len(out_r)), "n_total": int(mask_r.sum()),
                                        "r_values_outside": [float(x) for x in out_r]}
    l2_vs_poisson["euler"] = base.l2_over_range(e_real, poisson_euler.mean(axis=0), R_GRID, STAT_R_RANGE)

    report["poisson_null_i"] = {
        "method": "angles drawn uniformly (in solid angle) from the HEALPix-occupied footprint via "
                  "rejection sampling; r drawn uniform-in-comoving-volume within the SAME shell window "
                  "[r_lo,r_hi] used to select the data (E5's SELECTION_SHELL_IDX). Independent of the "
                  "data's actual (ra,dec) pairs and independent of its actual n(z) beyond the flat-in-"
                  "volume shell window -- unlike E5's z-shuffle null, this ALSO randomises angle.",
        "n_realizations": N_POISSON_NULL_SEEDS, "n_per_realization": N_FINAL,
        "per_realization_diagnostics": poisson_diag,
        "l2_real_vs_poisson_mean": l2_vs_poisson,
        "r_points_outside_95pct_envelope": outside_poisson,
    }
    log("Poisson null done. L2 vs poisson: %s" % l2_vs_poisson)

    # ================= xi(r): data, LS random, CAMB fit =================
    log("building LS random catalogue (%dx N_FINAL) for xi(r)" % XI_RANDOM_OVERSAMPLE)
    n_rand = XI_RANDOM_OVERSAMPLE * N_FINAL
    ra_r, dec_r, r_r = poisson_footprint_catalogue(occ_set, NSIDE_MASK, ra_min, ra_max, dec_min, dec_max,
                                                     r_lo, r_hi, n_rand, seed=6000)
    rand_xyz = base.radec_r_to_xyz(ra_r, dec_r, r_r)

    r_edges_xi = np.linspace(2.0, 70.0, 35)
    r_mid_xi = 0.5 * (r_edges_xi[:-1] + r_edges_xi[1:])
    log("precomputing RR (random-random pairs, done ONCE, reused for data + every mock's xi(r))")
    rtree = precompute_random_tree(rand_xyz, r_edges_xi, verbose=True)
    log("computing data LS xi(r) (cKDTree pair counts, N_data=%d, N_rand=%d)" % (real_xyz.shape[0], rand_xyz.shape[0]))
    xi_data, DD, RR, DR = ls_xi_with_tR(real_xyz, rtree, verbose=True)

    log("CAMB linear P(k) at z_eff=%.2f (Planck 2018 params, imported)" % Z_EFF)
    camb_out = camb_linear_pk(Z_EFF)
    kh, pk_lin = camb_out["kh"], camb_out["pk_zeff"]
    xi_lin = pk2xi(kh, pk_lin, r_mid_xi)

    fit_mask = (r_mid_xi >= XI_FIT_RANGE_MPC_H[0]) & (r_mid_xi <= XI_FIT_RANGE_MPC_H[1]) & np.isfinite(xi_data)
    b2 = float(np.sum(xi_data[fit_mask] * xi_lin[fit_mask]) / np.sum(xi_lin[fit_mask] ** 2))
    b2 = max(b2, 1e-6)
    b_fit = float(np.sqrt(b2))
    log("bias fit on xi(r), r in %s Mpc/h: b = %.4f (b^2=%.4f)" % (XI_FIT_RANGE_MPC_H, b_fit, b2))

    report["xi_r_and_bias"] = {
        "estimator": "Landy-Szalay, cKDTree.count_neighbors exact pair counts",
        "r_edges_mpc_over_h": r_edges_xi.tolist(), "xi_data": [None if not np.isfinite(x) else float(x) for x in xi_data],
        "xi_lin_matter_camb": xi_lin.tolist(),
        "fit_range_mpc_over_h": list(XI_FIT_RANGE_MPC_H),
        "fit_method": "b^2 = sum(xi_data*xi_lin)/sum(xi_lin^2) over the fit range ONLY, ordinary least squares; "
                      "NEVER fit to any topological statistic (per ground rules).",
        "b_fit": b_fit, "b2_fit": b2,
        "sanity_note": ("CORRECTED after an advisor review of this run's first draft, which quoted the WRONG "
                         "band (b~1.1-1.3 is for the SDSS main sample as a WHOLE; this shell, z~0.055-0.08, "
                         "is flux-limited to sub-L* (M_r~-19.5), faint-dominated galaxies, for which "
                         "b~0.9-1.0 is the more relevant literature range, L, quoted as a sanity band not a "
                         "fetched table). b_fit=%.3f is %s that (corrected) band. SEPARATELY: xi_data here is "
                         "measured in REDSHIFT SPACE (positions built from observed z, not a real-space "
                         "reconstruction), while xi_lin is a REAL-SPACE model -- the Kaiser boost factor "
                         "(1+2beta/3+beta^2/5) at beta=f/b~0.59 is ~1.47, so the REAL-SPACE-equivalent bias "
                         "is b_fit/sqrt(1.47)~%.2f, lower still than the number fit here. This mismatch (real-"
                         "space model vs redshift-space measurement) is a known simplification of this "
                         "pipeline, not corrected for, and is the more likely explanation for a low b_fit "
                         "than a normalisation error: this run's own sigma8 self-check (computing sigma8 from "
                         "the SAME kh/pk arrays and quadrature used for xi_lin, via a top-hat window, and "
                         "comparing to CAMB's own reported sigma8(z=0)) matched to 0.03%%, ruling out a gross "
                         "normalisation bug in pk2xi/xi2pk.") %
                       (b_fit, "inside" if 0.85 <= b_fit <= 1.3 else "OUTSIDE -- treat mock comparison with caution",
                        b_fit / np.sqrt(1.47)),
        "camb": {"z_eff_requested": Z_EFF, "z_eff_actual_from_camb": camb_out["z_eff_actual"],
                 "sigma8_z0_raw_from_As": camb_out["sigma8_z0_raw"],
                 "sigma8_z0_rescale_factor_to_hit_planck_target": camb_out["sigma8_z0_rescale_factor"],
                 "sigma8_selfcheck_own_quadrature_vs_camb_ratio": camb_out["sigma8_selfcheck_ratio_to_camb"],
                 "sigma8_selfcheck_note": "sigma8 recomputed from this script's own kh/pk(z=0) arrays with a "
                                          "top-hat window, R=8 Mpc/h, using the SAME trapezoidal quadrature "
                                          "style as pk2xi/xi2pk, and compared to CAMB's own reported sigma8. "
                                          "This checks ABSOLUTE normalisation, which the pk2xi<->xi2pk round-"
                                          "trip test elsewhere in this script cannot (a constant multiplicative "
                                          "error cancels in a round trip). A ratio far from 1 would mean b_fit "
                                          "and sigma_G^2 are both built on a mis-normalised P(k).",
                 "omega_m_arithmetic_check_ombh2_plus_omch2_over_h2": camb_out["omega_m_arith_check"],
                 "omega_m_pars_object_incl_neutrino": camb_out["omega_m_pars"],
                 "omega_m_m0_target": OMEGA_M,
                 "omega_m_consistency_pct": 100 * abs(camb_out["omega_m_arith_check"] - OMEGA_M) / OMEGA_M,
                 "planck18_params_imported": camb_out["planck18_params_imported"]},
    }

    # ================= lognormal field: xi_G, P_G, grid =================
    log("computing xi_G = ln(1+xi_target) and P_G(k) via Hankel transform (grid-band-limited)")
    k_nyq = 1.6 * np.pi / GRID_CELL_MPC_H  # grid resolution limit -- see smoothing note below
    R_SMOOTH_MPC_H = GRID_CELL_MPC_H       # Gaussian smoothing scale = 1 grid cell, see below
    # SMOOTH P_lin at the grid's own cell scale BEFORE transforming to xi_target(r), with a Gaussian
    # window exp(-(k*R_smooth)^2/2), R_smooth = 1 grid cell. Two earlier attempts without this both
    # failed: (1) feeding the FULL (unsmoothed) linear P(k) up to CAMB's own kmax=12 h/Mpc into
    # pk2xi gave xi_target(r->0.3 Mpc/h) of order thousands and sigma_G^2~8.4 (realistic galaxy-
    # scale sigma_G^2 is order 0.1-3); (2) a HARD cutoff at k_nyq alone only reduced sigma_G^2 to
    # ~7.4, because a sharp cutoff leaves a spike of un-damped power right at k_nyq that still
    # dominates the k^2 P(k) variance integral. In BOTH cases, a single finite-grid realization of
    # that skewed a field under-samples the lognormal transform's rare huge peaks, so the REALISED
    # mean density came out ~3x below the target n_bar every time -- exactly the shortfall this
    # run's own calibration test still showed even after Monte-Carlo-calibrating n_bar for the pure
    # geometric mask+shell fraction. The Gaussian window (a standard finite-cell smoothing choice)
    # gives sigma_G^2~2.7 here (this run's scratch calibration) while leaving xi_target UNCHANGED
    # at the r=20-60 Mpc/h scales the actual comparison uses (verified: <1% change there).
    smoothing_window = np.exp(-0.5 * (kh * R_SMOOTH_MPC_H) ** 2)
    pk_lin_band = pk_lin * smoothing_window * (kh <= k_nyq)
    r_fine = np.logspace(np.log10(0.3), np.log10(600.0), 700)
    xi_target_fine = b2 * pk2xi(kh, pk_lin_band, r_fine)
    xi_g_fine = np.log1p(np.clip(xi_target_fine, -0.999, None))
    k_for_pg = np.logspace(-3, np.log10(k_nyq), 300)
    pg_of_k = xi2pk(r_fine, xi_g_fine, k_for_pg)
    pg_of_k = np.clip(pg_of_k, 0.0, None)
    sigma_g2 = float(np.trapezoid(k_for_pg ** 2 * pg_of_k, k_for_pg) / (2 * np.pi ** 2))
    report["lognormal_field"] = {
        "recipe": "P_target=b^2*P_lin*exp(-(k*R_smooth)^2/2)*[k<=k_nyq], R_smooth=1 grid cell -> "
                  "xi_target=pk2xi(P_target) -> xi_G=ln(1+xi_target) -> "
                  "P_G=xi2pk(xi_G) -> delta_G~GRF(P_G) -> delta_LN=exp(delta_G-sigma_G^2/2)-1",
        "sigma_g2": sigma_g2,
        "k_nyquist_band_limit_h_over_mpc": k_nyq,
        "gaussian_smoothing_scale_mpc_over_h": R_SMOOTH_MPC_H,
        "transform_method": "1D trapezoidal quadrature of the k^2 P(k) sinc(kr) Hankel integral (validated "
                             "in this run's scratch calibration against a toy power-law P(k): round-trip "
                             "recovery within a few percent; GRF normalisation validated the same way, "
                             "measured P(k) of a realized field matched the input P(k) to <6% on a single "
                             "64^3 realization).",
    }

    # ---- grid covering the survey wedge's own bounding box (not a symmetric cube) ----
    xyz_full = base.radec_r_to_xyz(sel["df_ra"], sel["df_dec"], sel["df_r"])
    in_window = (sel["df_r"] >= r_lo) & (sel["df_r"] < r_hi)
    bbox_min = xyz_full[in_window].min(axis=0) - GRID_PAD_MPC_H
    bbox_max = xyz_full[in_window].max(axis=0) + GRID_PAD_MPC_H
    box_size = bbox_max - bbox_min
    shape = tuple(int(np.ceil(s / GRID_CELL_MPC_H)) for s in box_size)
    box_volume = float(np.prod(box_size))
    n_target_prewedge = 4 * N_FINAL  # oversample margin (raised from an initial 3x, see calibrate_n_bar docstring)
    log("calibrating n_bar by Monte Carlo (measuring the mask+shell survival fraction directly, not analytically)")
    n_bar, measured_wedge_fraction = calibrate_n_bar(bbox_min, box_size, occ_set, NSIDE_MASK, edges, n_target_prewedge)
    report["lognormal_grid"] = {
        "cell_mpc_over_h": GRID_CELL_MPC_H, "shape": list(shape), "n_cells": int(np.prod(shape)),
        "box_size_mpc_over_h": box_size.tolist(), "box_volume_mpc_over_h3": box_volume,
        "measured_wedge_fraction_of_box": measured_wedge_fraction,
        "measured_wedge_fraction_method": "Monte Carlo, 300000 uniform probe points over the box run through "
                                           "the SAME mask_and_shell_keep() filter the mocks use (replaces an "
                                           "earlier analytic solid-angle formula that was wrong by a factor of "
                                           "~3 -- see calibrate_n_bar's docstring).",
        "n_bar_per_mpc_over_h3": n_bar, "n_expected_points_in_box": n_bar * box_volume,
        "n_target_prewedge": n_target_prewedge,
        "note": "grid covers only the survey wedge's padded bounding box (not a symmetric NxNxN cube spanning "
                "the whole survey depth), which keeps n_cells small and avoids wasting >75% of a cube's "
                "volume outside the wedge.",
    }
    log("grid shape %s (%d cells), n_bar=%.5f /Mpc^3h3, expected box points ~%.0f" %
        (shape, np.prod(shape), n_bar, n_bar * box_volume))

    # ================= mocks =================
    log("running %d lognormal mocks" % N_LOGNORMAL_MOCKS)
    mock_betti = {0: [], 1: [], 2: []}
    mock_euler = []
    mock_diag = []
    mock_xi_stack = []
    r_edges_xi_mock = r_edges_xi  # same bins as data, for the gate
    for j, seed in enumerate(range(4000, 4000 + N_LOGNORMAL_MOCKS)):
        field_g = gaussian_random_field(k_for_pg, pg_of_k, shape, GRID_CELL_MPC_H, seed)
        rng = np.random.RandomState(seed + 1)
        pts = lognormal_sample_points(field_g, sigma_g2, n_bar, bbox_min, GRID_CELL_MPC_H, rng)
        del field_g
        gc.collect()
        final_xyz, n_avail = mock_to_final_selection(pts, occ_set, NSIDE_MASK, edges)
        del pts
        gc.collect()
        if final_xyz is None:
            mock_diag.append({"seed": seed, "status": "FAILED_insufficient_points", "n_available_in_window": n_avail})
            log("  mock seed=%d FAILED: only %d points available in mask+shell window (< %d needed)" %
                (seed, n_avail, len(SELECTION_SHELL_IDX) * N_TARGET_PER_SHELL))
            continue
        res_m, bd_m = base.alpha_persistence(final_xyz, MAX_ALPHA_SQ, "mock_seed%d" % seed, HERE)
        e_m, b0_m, b1_m, b2_m = base.euler_curve(bd_m, R_GRID, R_MAX_PERS)
        mock_betti[0].append(b0_m); mock_betti[1].append(b1_m); mock_betti[2].append(b2_m)
        mock_euler.append(e_m)
        mock_diag.append({"seed": seed, "n_points": final_xyz.shape[0], "n_available_in_window": n_avail,
                           "betti_numbers_at_truncation": res_m["betti_numbers_at_truncation"],
                           "runtime_sec": res_m["runtime_sec"]})
        if j < N_GATE_MOCKS:  # xi(r) gate computed on the first N_GATE_MOCKS mocks only (RR is shared, cheap)
            xi_m, _, _, _ = ls_xi_with_tR(final_xyz, rtree)
            mock_xi_stack.append(xi_m)
        del bd_m, final_xyz
        gc.collect()
        if (j + 1) % 4 == 0:
            log("  mock %d/%d done" % (j + 1, N_LOGNORMAL_MOCKS))

    n_ok = len(mock_betti[0])
    for d in (0, 1, 2):
        mock_betti[d] = np.array(mock_betti[d])
    mock_euler = np.array(mock_euler)
    mock_xi_stack = np.array(mock_xi_stack) if mock_xi_stack else np.empty((0, len(r_mid_xi)))
    # save the full mock Betti-curve stacks (n_ok x 3 x len(R_GRID), negligible size) for full provenance --
    # an advisor review of an earlier run noted only aggregate L2/envelope statistics were kept, so the
    # rank p-value's margin could not be inspected or recomputed after the fact.
    np.savez(os.path.join(HERE, "e5b_mock_betti_stacks.npz"), r_grid=R_GRID,
             b0=mock_betti[0], b1=mock_betti[1], b2=mock_betti[2], euler=mock_euler,
             b0_real=b0_real, b1_real=b1_real, b2_real=b2_real, euler_real=e_real)

    # ================= gate: recovered xi(r) vs data xi(r) =================
    n_gated = mock_xi_stack.shape[0]
    xi_mock_mean = np.nanmean(mock_xi_stack, axis=0) if mock_xi_stack.size else np.full_like(r_mid_xi, np.nan)
    # sigma_mock_field: the MOCK-TO-MOCK (per-realization) scatter -- how much ONE realization of xi(r)
    # varies at fixed cosmology/volume (cosmic variance). This, NOT the standard error of the mock mean,
    # is the right yardstick for comparing the DATA (itself a single realization) to the mock ensemble.
    # An EARLIER version of this gate used xi_mock_sem = sigma_mock_field/sqrt(n_gated) instead, which
    # tests only "is the ensemble MEAN biased relative to the data", not "is the data's one realization
    # consistent with the ensemble" -- an advisor review of that first full run caught this: RMS(z) rose
    # from 3.23 (6 mocks) to 3.81 (22 mocks), the diagnostic signature of a 1/sqrt(n) statistic (a
    # correctly calibrated consistency test should be roughly n-independent), and the corrected version
    # below gives RMS(z)~0.83 on the SAME 22 mocks -- a factor of sqrt(22)~4.7, exactly as the sem-based
    # bug predicts. Kept conservative: lognormal mocks are known to UNDERESTIMATE real non-Gaussian
    # covariance, so sigma_mock_field here is if anything too small, making the corrected z too large
    # (i.e. this gate is not tuned to pass easily).
    sigma_mock_field = (np.nanstd(mock_xi_stack, axis=0, ddof=1)
                         if n_gated > 1 else np.full_like(r_mid_xi, np.nan))
    sigma_data_shotnoise = (1.0 + xi_data) / np.sqrt(np.maximum(DD, 1.0))  # standard LS shot-noise error (Landy & Szalay 1993 form)
    # Var(mock_mean - data) = Var(one mock realization)*(1 + 1/n_gated) [finite-sample correction for
    # comparing to a MEAN of n_gated draws from the same distribution the data itself is assumed to be
    # one draw from] + data's own shot-noise variance.
    sigma_total = np.sqrt(sigma_mock_field ** 2 * (1.0 + 1.0 / max(n_gated, 1)) + sigma_data_shotnoise ** 2)
    ratio = xi_mock_mean[fit_mask] / xi_data[fit_mask]  # kept for the plot/record; NOT the pass/fail criterion (see GATE_Z_THRESHOLD)
    z = (xi_mock_mean[fit_mask] - xi_data[fit_mask]) / sigma_total[fit_mask]
    z_finite = z[np.isfinite(z)]
    rms_z = float(np.sqrt(np.mean(z_finite ** 2))) if z_finite.size else float("nan")
    gate_pass = bool(np.isfinite(rms_z) and rms_z < GATE_Z_THRESHOLD)
    report["lognormal_gate"] = {
        "definition": "RMS over r in %s Mpc/h of z_i=(xi_mock_mean_i-xi_data_i)/sigma_total_i, where "
                      "sigma_total_i = sqrt(sigma_mock_field_i^2*(1+1/n_gated) + sigma_data_shotnoise_i^2), "
                      "sigma_mock_field is the mock-TO-MOCK (per-realization, NOT standard-error-of-the-mean) "
                      "scatter over the %d gated mocks, and sigma_data_shotnoise is the data's own Poisson "
                      "shot-noise error ((1+xi)/sqrt(DD), Landy & Szalay 1993). PASS iff RMS(z) < %.1f. This "
                      "gate MUST pass for the Betti-curve comparison below to be treated as diagnostic of "
                      "anything (ground rules / advisor review requirement). A raw amplitude ratio is also "
                      "reported for the record but is NOT the pass/fail criterion (it diverges at this fit "
                      "range's zero-crossings)." % (XI_FIT_RANGE_MPC_H, n_gated, GATE_Z_THRESHOLD),
        "n_mocks_used_for_gate": n_gated,
        "r_mid_fit_range": r_mid_xi[fit_mask].tolist(),
        "xi_data_fit_range": xi_data[fit_mask].tolist(),
        "xi_mock_mean_fit_range": xi_mock_mean[fit_mask].tolist(),
        "sigma_data_shotnoise_fit_range": sigma_data_shotnoise[fit_mask].tolist(),
        "sigma_mock_field_scatter_fit_range": sigma_mock_field[fit_mask].tolist(),
        "z_score_fit_range": z.tolist(),
        "rms_z": rms_z,
        "ratio_mock_over_data_for_record_only": ratio.tolist(),
        "pass": gate_pass,
    }
    log("lognormal xi(r) gate: PASS=%s, RMS(z)=%.3f (n_gate_mocks=%d)" % (gate_pass, rms_z, n_gated))

    # ================= Betti comparison + rank p-value (pre-stated stat) =================
    l2_vs_mock = {}
    outside_mock = {}
    for d, curve in zip((0, 1, 2), (b0_real, b1_real, b2_real)):
        if mock_betti[d].shape[0] == 0:
            l2_vs_mock["H%d" % d] = None
            continue
        mean_m = mock_betti[d].mean(axis=0)
        l2_vs_mock["H%d" % d] = base.l2_over_range(curve, mean_m, R_GRID, STAT_R_RANGE)
        lo = np.percentile(mock_betti[d], 2.5, axis=0); hi = np.percentile(mock_betti[d], 97.5, axis=0)
        mask_r = (R_GRID >= STAT_R_RANGE[0]) & (R_GRID <= STAT_R_RANGE[1])
        out_r = R_GRID[mask_r][(curve[mask_r] < lo[mask_r]) | (curve[mask_r] > hi[mask_r])]
        outside_mock["H%d" % d] = {"n_outside": int(len(out_r)), "n_total": int(mask_r.sum()),
                                     "r_values_outside": [float(x) for x in out_r]}
    if mock_euler.shape[0]:
        l2_vs_mock["euler"] = base.l2_over_range(e_real, mock_euler.mean(axis=0), R_GRID, STAT_R_RANGE)

    # rank p-value on H1 (the dimension the "filament/void topology" claim is about), pre-stated stat = L2 over STAT_R_RANGE
    # ALSO split into r<R_SPACING (shot-noise/discreteness dominated at this N, per the mean inter-galaxy
    # spacing computed below) vs r>=R_SPACING, so a reader can see how much of the headline p-value comes
    # from a scale below this sample's own resolution -- added after an advisor review noted this was an
    # unquantified caveat in an earlier version of this report.
    R_SPACING_DIAGNOSTIC = float((measured_wedge_fraction * box_volume / N_FINAL) ** (1 / 3))
    rank_pvalues = {}
    for d, curve in zip((0, 1, 2), (b0_real, b1_real, b2_real)):
        stack = mock_betti[d]
        if stack.shape[0] < 2:
            rank_pvalues["H%d" % d] = None
            continue
        loo_l2 = []
        loo_l2_lo = []
        loo_l2_hi = []
        r_range_lo = (STAT_R_RANGE[0], R_SPACING_DIAGNOSTIC)
        r_range_hi = (R_SPACING_DIAGNOSTIC, STAT_R_RANGE[1])
        for i in range(stack.shape[0]):
            others_mean = np.delete(stack, i, axis=0).mean(axis=0)
            loo_l2.append(base.l2_over_range(stack[i], others_mean, R_GRID, STAT_R_RANGE))
            loo_l2_lo.append(base.l2_over_range(stack[i], others_mean, R_GRID, r_range_lo))
            loo_l2_hi.append(base.l2_over_range(stack[i], others_mean, R_GRID, r_range_hi))
        loo_l2 = np.array(loo_l2)
        real_l2 = base.l2_over_range(curve, stack.mean(axis=0), R_GRID, STAT_R_RANGE)
        real_l2_lo = base.l2_over_range(curve, stack.mean(axis=0), R_GRID, r_range_lo)
        real_l2_hi = base.l2_over_range(curve, stack.mean(axis=0), R_GRID, r_range_hi)
        rank = int(np.sum(loo_l2 >= real_l2))  # how many mocks are AT LEAST as extreme as the real data
        rank_lo = int(np.sum(np.array(loo_l2_lo) >= real_l2_lo))
        rank_hi = int(np.sum(np.array(loo_l2_hi) >= real_l2_hi))
        p = (rank + 1) / (stack.shape[0] + 1)
        # margin diagnostic (added after advisor review: rank=0 alone does not say whether real_l2 is
        # barely above the mock spread or many multiples of it -- report the ratio so a reader can tell)
        margin_ratio = float(real_l2 / np.median(loo_l2)) if np.median(loo_l2) > 0 else float("inf")
        rank_pvalues["H%d" % d] = {"real_l2": real_l2, "n_mocks_at_least_as_extreme": rank,
                                     "n_mocks": int(stack.shape[0]), "p_value": p,
                                     "resolution_floor": 1.0 / (stack.shape[0] + 1),
                                     "loo_l2_min_median_max": [float(loo_l2.min()), float(np.median(loo_l2)), float(loo_l2.max())],
                                     "real_l2_over_median_loo_l2": margin_ratio,
                                     "r_split_diagnostic": {
                                         "r_spacing_mpc_over_h": R_SPACING_DIAGNOSTIC,
                                         "r_range_below_spacing": list(r_range_lo), "p_value_below_spacing": (rank_lo + 1) / (stack.shape[0] + 1),
                                         "r_range_at_or_above_spacing": list(r_range_hi), "p_value_at_or_above_spacing": (rank_hi + 1) / (stack.shape[0] + 1),
                                     }}

    report["lognormal_mocks_iii"] = {
        "n_requested": N_LOGNORMAL_MOCKS, "n_succeeded": n_ok, "n_failed": N_LOGNORMAL_MOCKS - n_ok,
        "per_mock_diagnostics": mock_diag,
        "l2_real_vs_mock_mean": l2_vs_mock,
        "r_points_outside_95pct_envelope": outside_mock,
        "rank_pvalue_leave_one_out": rank_pvalues,
        "interpretation": ("THIS (not the z-shuffle or Poisson-null envelopes above) is the actual M0 test "
                            "per the ground rules. p close to the resolution floor (rank=0, i.e. the real "
                            "L2 exceeds every leave-one-out mock L2) is the strongest evidence THIS TEST can "
                            "give against 'consistent with M0's linear P(k)+bias+lognormal field'; a mid-"
                            "range p is a match at this N, this scale range, and this (approximate) mock "
                            "recipe -- consistent with M0 at the level of this test, NOT a K3xT2 "
                            "confirmation. A mismatch concentrated below ~4-6 Mpc/h (this sample's mean "
                            "inter-galaxy spacing, ~%.1f Mpc/h) is expected: lognormal mocks reproduce the "
                            "1-point PDF and 2-point function but not higher-order correlations or true "
                            "void/filament shapes, and that scale is shot-noise dominated at N=%d in this "
                            "volume regardless." % ((measured_wedge_fraction * box_volume / N_FINAL) ** (1 / 3), N_FINAL)),
    }

    # ---- plot ----
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    labels = ["H0", "H1", "H2"]
    curves_real = (b0_real, b1_real, b2_real)
    for k in range(3):
        ax = axes[k]
        if mock_betti[k].shape[0]:
            lo = np.percentile(mock_betti[k], 2.5, axis=0); hi = np.percentile(mock_betti[k], 97.5, axis=0)
            ax.fill_between(R_GRID, lo, hi, color="tab:blue", alpha=0.3, label="lognormal mocks 95%% (n=%d)" % mock_betti[k].shape[0])
        if poisson_betti[k].shape[0]:
            lo2 = np.percentile(poisson_betti[k], 2.5, axis=0); hi2 = np.percentile(poisson_betti[k], 97.5, axis=0)
            ax.fill_between(R_GRID, lo2, hi2, color="gray", alpha=0.25, label="Poisson-footprint 95%%")
        ax.plot(R_GRID, curves_real[k], color="tab:red", lw=1.8, label="SDSS DR17 real")
        ax.set_title(labels[k]); ax.set_xlabel("r (Mpc/h)")
        if k == 0:
            ax.set_ylabel("beta_k(r)")
        ax.legend(fontsize=6)
    fig.suptitle("E5b: real vs lognormal-mock and Poisson-footprint null envelopes")
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "e5b_betti_real_vs_lognormal_and_poisson.png"), dpi=130)
    plt.close(fig)

    fig2, ax2 = plt.subplots(figsize=(6.5, 4.5))
    ax2.plot(r_mid_xi, xi_data, "o-", color="tab:red", label="data LS xi(r)")
    if mock_xi_stack.size:
        ax2.plot(r_mid_xi, xi_mock_mean, "s--", color="tab:blue", label="mean mock LS xi(r) (n=%d)" % mock_xi_stack.shape[0])
    ax2.plot(r_mid_xi, b2 * xi_lin, "k:", label="b^2 * xi_lin (CAMB, fit)")
    ax2.axvspan(*XI_FIT_RANGE_MPC_H, color="gray", alpha=0.15, label="fit range")
    ax2.set_xscale("log"); ax2.set_yscale("log")
    ax2.set_xlabel("r (Mpc/h)"); ax2.set_ylabel("xi(r)"); ax2.legend(fontsize=8)
    ax2.set_title("E5b: two-point function, data vs CAMB-lognormal mocks")
    fig2.tight_layout()
    fig2.savefig(os.path.join(HERE, "e5b_xi_r_data_vs_mocks.png"), dpi=130)
    plt.close(fig2)

    report["total_runtime_sec"] = time.time() - T0
    out_path = os.path.join(HERE, "e5b_lognormal_and_poisson_null_report.json")
    with open(out_path, "w") as f:
        json.dump(report, f, indent=2)
    log("wrote %s" % out_path)
    log("TOTAL RUNTIME %.1f sec" % report["total_runtime_sec"])
    return report


if __name__ == "__main__":
    run()
