#!/usr/bin/env python3
"""Shared library for audit/cosmic_vorticity (registration.json id CV).

Implements the registered detectors CV-CMB-D1 / CV-DESI-D1 and the three
registered statistics S1 (H0 death-radius spread), S2 (count), S3
(orientational order).  Every hyperparameter here is the one frozen in
registration.json, which was committed (7ecf289) before any map or catalogue
was read.

The FIXED TDA library is at ../lib/ (vendored from loop/tda-simple commit
d8175f1).  audit/reverse_zero/E5-cmb-tda/cmb_tda.py is NOT imported.

Tier of every number: X (numerics).
"""
from __future__ import annotations

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
CV = os.path.dirname(HERE)
WT_ROOT = os.path.dirname(os.path.dirname(CV))
sys.path.insert(0, CV)          # so `import lib.stats` works
RESULTS = os.path.join(CV, "results")
CACHE = os.path.join(CV, "cache")

# ---------------------------------------------------------------- registered constants
NSIDE = 128
LMAX = 3 * NSIDE                       # 384, X1's registered lmax
SIGMA_S_DEG = 1.0                      # registered; chosen on sims 9000001/2
R_DISK_DEG = 2.0 * SIGMA_S_DEG         # = a/3 when sigma = a/6  (Re6Zr ratio exactly 2)
NU_PRIMARY = 1.0
NU_SECONDARY = 1.5
A_REF_DEG = 6.0 * SIGMA_S_DEG          # inverse of Re6Zr's sigma = a/6
GATE_THRESH = 2.8070337683438042       # scipy.stats.norm.isf(0.025/10), X1's threshold
PIXEL_SCALE_DEG = float(np.degrees(np.sqrt(4 * np.pi / (12 * NSIDE ** 2))))  # 0.45815

# DESI (registered)
L_CELL = 8.0                           # Mpc/h
SIGMA_G = 16.0                         # Mpc/h
R_EX = 2.0 * SIGMA_G                   # 32 Mpc/h
A_REF_MPC = 6.0 * SIGMA_G              # 96 Mpc/h
NU_DESI = 1.0

X1 = os.path.join(WT_ROOT, "audit/reverse_zero_r2/X1-cmb")
MAPS = {
    "wmap": dict(map=os.path.join(WT_ROOT, "data/real2/cmb/wmap_ilc_9yr_v5.fits"),
                 mask=os.path.join(WT_ROOT, "data/real2/cmb/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits")),
    "smica": dict(map=os.path.join(WT_ROOT, "data/real2/cmb/planck_smica_2048_R3.fits"),
                  mask=os.path.join(WT_ROOT, "data/real2/cmb/planck_common_mask_int_2048_R3.fits")),
}


def deg2chord(d_deg):
    """chord length on the unit sphere for a great-circle separation in degrees."""
    return 2.0 * np.sin(np.radians(d_deg) / 2.0)


def chord2deg(c):
    return np.degrees(2.0 * np.arcsin(np.clip(np.asarray(c, float) / 2.0, -1.0, 1.0)))


MAX_ALPHA_SQ_CMB = (3.0 * deg2chord(A_REF_DEG)) ** 2
R_TRUNC_CMB_CHORD = 3.0 * deg2chord(A_REF_DEG)
MAX_ALPHA_SQ_DESI = (3.0 * A_REF_MPC) ** 2


# ---------------------------------------------------------------- CMB I/O (X1 verbatim)
def load_data(which):
    """X1's load_data, verbatim: ud_grade the map to nside 128 and threshold the
    ud_graded mask at 0.5."""
    import healpy as hp
    cfg = MAPS[which]
    m = hp.read_map(cfg["map"], field=0)           # read_map converts NESTED -> RING
    mk = hp.read_map(cfg["mask"], field=0)
    m = np.asarray(m, dtype=np.float64)
    m[m < -1e20] = 0.0
    m128 = hp.ud_grade(m, NSIDE)
    mask128 = (hp.ud_grade(np.asarray(mk, dtype=np.float64), NSIDE) >= 0.5).astype(np.uint8)
    return m128, mask128


def prep_map(t, mask):
    """X1's prep_map, verbatim: remove monopole+dipole fitted on unmasked pixels.
    The analogue of Re6Zr's plane subtraction."""
    import healpy as hp
    unm = mask > 0
    x, y, z = hp.pix2vec(NSIDE, np.arange(hp.nside2npix(NSIDE)))
    A = np.stack([np.ones_like(x), x, y, z], axis=1)
    coef, *_ = np.linalg.lstsq(A[unm], t[unm], rcond=None)
    return t - A @ coef


def make_sim(cl_in, seed):
    """X1's make_sim, verbatim."""
    import healpy as hp
    np.random.seed(seed)
    alm = hp.synalm(cl_in, lmax=LMAX, new=True)
    return hp.alm2map(alm, NSIDE, lmax=LMAX, pixwin=False)


def load_cl_in(which):
    d = np.load(os.path.join(X1, which, "cl_in.npz"))
    return d["cl_in"], float(d["fsky"])


def smooth_map(t, mask, sigma_deg=SIGMA_S_DEG):
    """Mask applied ONCE, then band-limited Gaussian smoothing at lmax=384."""
    import healpy as hp
    tm = np.where(mask > 0, t, 0.0)
    alm = hp.map2alm(tm, lmax=LMAX, iter=1)
    alm = hp.smoothalm(alm, sigma=np.radians(sigma_deg), inplace=False)
    return hp.alm2map(alm, NSIDE, lmax=LMAX, pixwin=False)


# ---------------------------------------------------------------- disk neighbourhood cache
def disk_table(radius_deg=R_DISK_DEG):
    """Ragged->padded table of the HEALPix pixels inside a geodesic disk of the
    given radius around every pixel.  Purely geometric, so it is computed once
    and cached; the local-max test is then fully vectorised and identical for
    the data, every null, every injection and every control."""
    import healpy as hp
    os.makedirs(CACHE, exist_ok=True)
    fn = os.path.join(CACHE, "disk_nside%d_r%.3fdeg.npz" % (NSIDE, radius_deg))
    if os.path.exists(fn):
        d = np.load(fn)
        return d["tab"], d["nnb"]
    npix = hp.nside2npix(NSIDE)
    vec = np.array(hp.pix2vec(NSIDE, np.arange(npix))).T
    R = np.radians(radius_deg)
    lists = [hp.query_disc(NSIDE, vec[p], R) for p in range(npix)]
    w = max(len(x) for x in lists)
    tab = np.full((npix, w), -1, dtype=np.int32)
    nnb = np.zeros(npix, dtype=np.int32)
    for p, x in enumerate(lists):
        tab[p, :len(x)] = x
        nnb[p] = len(x)
    np.savez_compressed(fn, tab=tab, nnb=nnb)
    return tab, nnb


def eroded_mask(mask, tab):
    """A pixel survives erosion iff EVERY pixel of its R_disk disk is unmasked,
    so no candidate is an artefact of a truncated neighbourhood."""
    m = (mask > 0)
    # tab is padded with -1 (disks hold 56..64 pixels at nside 128, R = 2 deg).
    # The sentinel must be True so that PADDING is not mistaken for a masked
    # neighbour; real masked neighbours are False and still veto the pixel.
    padded = np.concatenate([m, [True]])
    return padded[tab].all(axis=1) & m


# ---------------------------------------------------------------- detector CV-CMB-D1
def detect_cmb(t_prepped, mask, tab, eroded, nu=NU_PRIMARY, sigma_deg=SIGMA_S_DEG,
               sign=0, ts_precomputed=None):
    """Registered detector CV-CMB-D1.

    A pixel p inside the eroded valid region is a candidate iff
      (i)  |T_s(p)| > nu * sigma(T_s over unmasked pixels of THIS realisation)
      (ii) |T_s(p)| >= max |T_s| over the geodesic disk of radius R_disk around p.

    sign = 0 -> |T_s| (registered primary); +1 -> maxima of T_s only;
    -1 -> minima only (registered secondary diagnostics).
    Returns (pixel indices, sigma_Ts, T_s).
    """
    ts = smooth_map(t_prepped, mask, sigma_deg) if ts_precomputed is None else ts_precomputed
    sd = float(ts[mask > 0].std())
    if sign == 0:
        f = np.abs(ts)
    elif sign > 0:
        f = ts.copy()
    else:
        f = -ts
    pad = np.concatenate([f, [-np.inf]])
    local_max = f >= pad[tab].max(axis=1)
    sel = eroded & local_max & (f > nu * sd)
    return np.where(sel)[0], sd, ts


# ---------------------------------------------------------------- statistics
def alpha_h0_deaths(points, max_alpha_sq):
    """gudhi.AlphaComplex, Z/2, min_persistence=0.0, bars in RADIUS units --
    the call quoted verbatim from cosmic_web_tda_scaled.alpha_persistence."""
    import gudhi
    pts = np.asarray(points, float)
    if len(pts) < 4:
        return np.empty(0), [0, 0, 0]
    ac = gudhi.AlphaComplex(points=pts.tolist())
    st = ac.create_simplex_tree(max_alpha_square=float(max_alpha_sq))
    diag = st.persistence(homology_coeff_field=2, min_persistence=0.0)
    d0 = np.array([np.sqrt(max(d, 0.0)) for dim, (b, d) in diag if dim == 0 and np.isfinite(d)])
    return d0, list(st.betti_numbers())


def spread_stats(x):
    """qf_common.spread_stats, verbatim."""
    x = np.asarray(x, float)
    if x.size == 0:
        return {"n": 0}
    q1, med, q3 = np.percentile(x, [25, 50, 75])
    return {"n": int(x.size), "median": float(med), "q1": float(q1), "q3": float(q3),
            "iqr_over_median": float((q3 - q1) / med) if med > 0 else None,
            "mean": float(x.mean()), "cv": float(x.std() / x.mean()) if x.mean() > 0 else None}


def psi6_sphere(unit_vectors, a_ref_deg=A_REF_DEG, cut=1.5):
    """S3 for the sphere: site-averaged hexatic order over SPHERICAL Delaunay
    bonds (= scipy.spatial.ConvexHull of the unit vectors), bearings measured in
    each site's local tangent frame, bonds longer than cut*a_ref dropped, sites
    with fewer than 3 remaining bonds dropped.  Registered substitution for
    Re6Zr's global |<exp(6 i theta_b)>|, which vanishes by symmetry on a sphere.
    """
    from scipy.spatial import ConvexHull
    u = np.asarray(unit_vectors, float)
    n = len(u)
    if n < 8:
        return float("nan"), 0, 0
    hull = ConvexHull(u)
    nb = [set() for _ in range(n)]
    for s in hull.simplices:
        for i in s:
            for j in s:
                if i != j:
                    nb[i].add(int(j))
    cut_cos = np.cos(np.radians(cut * a_ref_deg))
    zhat = np.array([0.0, 0.0, 1.0])
    xhat = np.array([1.0, 0.0, 0.0])
    vals, nbond = [], []
    for j in range(n):
        k = np.fromiter(nb[j], dtype=int)
        if k.size == 0:
            continue
        c = u[k] @ u[j]
        k = k[c >= cut_cos]
        if k.size < 3:
            continue
        ref = zhat if abs(u[j, 2]) < 0.99 else xhat
        e1 = ref - (ref @ u[j]) * u[j]
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(u[j], e1)
        t = u[k] - np.outer(u[k] @ u[j], u[j])
        beta = np.arctan2(t @ e2, t @ e1)
        vals.append(abs(np.mean(np.exp(6j * beta))))
        nbond.append(k.size)
    if not vals:
        return float("nan"), 0, 0
    return float(np.mean(vals)), len(vals), int(np.sum(nbond))


def q6_3d(points, a_ref=A_REF_MPC, cut=1.5):
    """S3 for 3-D: Steinhardt Q6 per site over neighbours within cut*a_ref,
    sites with fewer than 4 neighbours dropped; statistic = mean_j Q6_j."""
    from scipy.spatial import cKDTree
    from scipy.special import sph_harm_y
    p = np.asarray(points, float)
    if len(p) < 8:
        return float("nan"), 0
    tree = cKDTree(p)
    pairs = tree.query_ball_point(p, r=cut * a_ref)
    vals = []
    for j, idx in enumerate(pairs):
        idx = [i for i in idx if i != j]
        if len(idx) < 4:
            continue
        d = p[idx] - p[j]
        r = np.linalg.norm(d, axis=1)
        theta = np.arccos(np.clip(d[:, 2] / r, -1, 1))
        phi = np.arctan2(d[:, 1], d[:, 0])
        q = np.array([np.mean(sph_harm_y(6, m, theta, phi)) for m in range(-6, 7)])
        vals.append(np.sqrt(4 * np.pi / 13 * np.sum(np.abs(q) ** 2)))
    if not vals:
        return float("nan"), 0
    return float(np.mean(vals)), len(vals)


def empirical_two_sided_p(null, value):
    """Registered convention (X2's): p = min(1, 2*min(le, ge)) with
    le = (1 + #{null <= v})/(N+1), ge = (1 + #{null >= v})/(N+1)."""
    a = np.asarray(null, float)
    a = a[np.isfinite(a)]
    n = a.size
    if n == 0 or not np.isfinite(value):
        return None
    lo = (1 + int((a <= value).sum())) / (n + 1.0)
    hi = (1 + int((a >= value).sum())) / (n + 1.0)
    return float(min(1.0, 2 * min(lo, hi)))


def cmb_point_statistics(pix, tab_unused=None):
    """S1, S2, S3 and the absolute floor for a CMB pixel list, on 3-D chords."""
    import healpy as hp
    u = np.array(hp.pix2vec(NSIDE, pix)).T if len(pix) else np.zeros((0, 3))
    out = {"S2_count": int(len(pix))}
    d0, betti = alpha_h0_deaths(u, MAX_ALPHA_SQ_CMB)
    ss = spread_stats(d0)
    out["h0"] = ss
    out["S1_iqr_over_median"] = ss.get("iqr_over_median")
    out["betti_at_truncation"] = betti
    if d0.size:
        dg = chord2deg(d0)
        out["absolute_floor"] = {
            "median_death_deg": float(np.median(dg)),
            "q1_death_deg": float(np.percentile(dg, 25)),
            "q3_death_deg": float(np.percentile(dg, 75)),
            "twice_median_death_deg": float(2 * np.median(dg)),
            "pixel_scale_deg": PIXEL_SCALE_DEG,
            "n_bars_below_pixel_scale": int((dg < PIXEL_SCALE_DEG).sum()),
            "frac_bars_below_pixel_scale": float((dg < PIXEL_SCALE_DEG).mean()),
            "max_death_deg": float(dg.max()),
        }
        out["truncation"] = {
            "r_trunc_chord": float(R_TRUNC_CMB_CHORD),
            "r_trunc_deg": float(chord2deg(R_TRUNC_CMB_CHORD)),
            "frac_h0_deaths_at_or_above_r_trunc": float((d0 >= R_TRUNC_CMB_CHORD * (1 - 1e-9)).mean()),
        }
    s3, nsite, nbond = psi6_sphere(u)
    out["S3_psi6_site_mean"] = s3
    out["S3_n_sites_kept"] = nsite
    out["S3_n_bonds"] = nbond
    return out
