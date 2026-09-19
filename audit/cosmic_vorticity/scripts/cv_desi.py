#!/usr/bin/env python3
"""TEST B of registration.json (CV): the Re6Zr defect-core pipeline on DESI DR1
BGS_BRIGHT-21.5 (volume-limited), NGC primary / SGC secondary.

Sub-commands:
  gate1    --cap NGC          G1 volume-limited check: comoving n(r) in 4 Mpc/h shells
  data     --cap NGC          the real catalogue through detector CV-DESI-D1
  control  --cap NGC --n 100  N2 official randoms, shell-matched, IDENTICAL detector (seeds 7100000+k)
  inject   --cap NGC --n 50   the amplitude x density injection grid (seeds 7300000+...)
  lognormal --cap NGC --n N   DECLARED POST-REGISTRATION DEVIATION: spectrum-matched
                              lognormal mocks on the same grid (see report.json)
  analyze  --cap NGC

Weights: w = WEIGHT * WEIGHT_FKP for galaxies AND randoms (registered).
Every command: timeout <s> prlimit --as=8589934592 -- <venv-tda python> cv_desi.py <sub> ...
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_lib as L  # noqa: E402

DESI = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/desi_dr1_lss"
OUT = os.path.join(L.CV, "results")
ZMIN, ZMAX = 0.10, 0.40
PAD = 40.0                       # Mpc/h, registered
SHELL = 4.0                      # Mpc/h, G1 shells


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), file=sys.stderr, flush=True)


def comoving_h(z):
    """Planck18 comoving distance in Mpc/h (the X2 convention, x2_lib.radecz_to_r),
    via an interpolation table so 13 M randoms are affordable."""
    from astropy.cosmology import Planck18
    zt = np.linspace(0.0, 0.5, 4001)
    rt = Planck18.comoving_distance(zt).value * Planck18.h
    return np.interp(z, zt, rt)


def read_cat(path, need_z=True):
    from astropy.io import fits
    with fits.open(path, memmap=True) as f:
        d = f[1].data
        z = np.asarray(d["Z"], dtype=np.float64)
        sel = (z >= ZMIN) & (z < ZMAX)
        ra = np.asarray(d["RA"], dtype=np.float64)[sel]
        dec = np.asarray(d["DEC"], dtype=np.float64)[sel]
        w = (np.asarray(d["WEIGHT"], dtype=np.float64)[sel] *
             np.asarray(d["WEIGHT_FKP"], dtype=np.float64)[sel])
        z = z[sel]
    r = comoving_h(z)
    cd = np.cos(np.radians(dec))
    xyz = np.stack([r * cd * np.cos(np.radians(ra)), r * cd * np.sin(np.radians(ra)),
                    r * np.sin(np.radians(dec))], 1)
    return xyz.astype(np.float32), w.astype(np.float32), r.astype(np.float32)


class Grid:
    """NGP deposit on an 8 Mpc/h grid (the X2 estimator convention), Gaussian
    smoothing at sigma_G = 16 Mpc/h, and the registered detector CV-DESI-D1."""

    def __init__(self, xyz_d, w_d, xyz_r, w_r):
        lo = xyz_d.min(0) - PAD
        hi = xyz_d.max(0) + PAD
        self.origin = lo
        self.shape = tuple(int(np.ceil((hi[i] - lo[i]) / L.L_CELL)) + 1 for i in range(3))
        self.n = int(np.prod(self.shape))
        self.Rs = ndimage.gaussian_filter(self.deposit(xyz_r, w_r), L.SIGMA_G / L.L_CELL, mode="constant")
        alpha = float(w_d.sum() / w_r.sum())
        self.alphaR = alpha * self.Rs
        pos = self.alphaR[self.alphaR > 0]
        self.valid = self.alphaR >= 0.5 * np.median(pos)
        # erosion by the R_ex ball: the whole ball must lie inside the valid region
        rad = int(round(L.R_EX / L.L_CELL))
        g = np.mgrid[-rad:rad + 1, -rad:rad + 1, -rad:rad + 1]
        self.ball = (g ** 2).sum(0) <= rad ** 2
        self.eroded = ndimage.binary_erosion(self.valid, structure=self.ball)
        # guard: never let a candidate sit within `rad+1` cells of a box face
        edge = np.ones(self.shape, bool)
        edge[:rad + 1] = edge[-(rad + 1):] = False
        edge[:, :rad + 1] = edge[:, -(rad + 1):] = False
        edge[:, :, :rad + 1] = edge[:, :, -(rad + 1):] = False
        self.eroded &= edge
        self.alpha = alpha
        self.rad = rad
        off = np.argwhere(self.ball) - rad
        st = np.array([self.shape[1] * self.shape[2], self.shape[2], 1])
        self.flat_off = (off * st).sum(1)

    def deposit(self, xyz, w):
        idx = np.floor((xyz - self.origin) / L.L_CELL + 0.5).astype(np.int64)
        np.clip(idx, 0, np.array(self.shape) - 1, out=idx)
        flat = (idx[:, 0] * self.shape[1] + idx[:, 1]) * self.shape[2] + idx[:, 2]
        return np.bincount(flat, weights=w.astype(np.float64), minlength=self.n).reshape(self.shape)

    def delta(self, xyz, w):
        Ds = ndimage.gaussian_filter(self.deposit(xyz, w), L.SIGMA_G / L.L_CELL, mode="constant")
        d = np.zeros(self.shape)
        np.divide(Ds - self.alphaR, self.alphaR, out=d, where=self.valid)
        return d

    def peaks(self, d, nu=L.NU_DESI):
        """Registered detector CV-DESI-D1: delta > nu*sigma(delta over valid
        cells) AND delta is the maximum of delta over the R_ex ball, with the
        whole ball inside the valid region.  The spherical-footprint maximum is
        the exact analogue of Re6Zr's ndimage.minimum_filter over a disk."""
        sd = float(d[self.valid].std())
        cand = np.flatnonzero(self.eroded.ravel() & (d.ravel() > nu * sd))
        if cand.size == 0:
            return np.zeros((0, 3)), sd, 0.0
        f = d.ravel()
        keep = np.ones(cand.size, bool)
        for o in self.flat_off:
            if o == 0:
                continue
            keep &= f[cand] >= f[cand + o]
        sel = cand[keep]
        i0, rem = np.divmod(sel, self.shape[1] * self.shape[2])
        i1, i2 = np.divmod(rem, self.shape[2])
        pts = self.origin + (np.stack([i0, i1, i2], 1)) * L.L_CELL
        return pts, sd, float(self.valid.sum()) * L.L_CELL ** 3


def desi_statistics(pts):
    out = {"S2_count": int(len(pts))}
    d0, betti = L.alpha_h0_deaths(pts, L.MAX_ALPHA_SQ_DESI)
    ss = L.spread_stats(d0)
    out["h0"] = ss
    out["S1_iqr_over_median"] = ss.get("iqr_over_median")
    out["betti_at_truncation"] = betti
    if d0.size:
        out["absolute_floor"] = {
            "median_death_Mpch": float(np.median(d0)), "q1_Mpch": float(np.percentile(d0, 25)),
            "q3_Mpch": float(np.percentile(d0, 75)), "twice_median_Mpch": float(2 * np.median(d0)),
            "cell_Mpch": L.L_CELL, "sigma_G_Mpch": L.SIGMA_G,
            "n_bars_below_cell": int((d0 < L.L_CELL).sum()),
            "n_bars_below_sigma_G": int((d0 < L.SIGMA_G).sum()),
            "frac_bars_below_sigma_G": float((d0 < L.SIGMA_G).mean()),
            "max_death_Mpch": float(d0.max())}
        rt = 3 * L.A_REF_MPC
        out["truncation"] = {"r_trunc_Mpch": rt,
                             "frac_h0_deaths_at_or_above_r_trunc": float((d0 >= rt * (1 - 1e-9)).mean())}
    q6, ns = L.q6_3d(pts)
    out["S3_Q6_site_mean"] = q6
    out["S3_n_sites_kept"] = ns
    return out


def _load(cap):
    xd, wd, rd = read_cat(os.path.join(DESI, "BGS_BRIGHT-21.5_%s_clustering.dat.fits" % cap))
    xr, wr, rr = read_cat(os.path.join(DESI, "BGS_BRIGHT-21.5_%s_0_clustering.ran.fits" % cap))
    return xd, wd, rd, xr, wr, rr


# --------------------------------------------------------------------- G1
def cmd_gate1(a):
    xd, wd, rd, xr, wr, rr = _load(a.cap)
    edges = np.arange(np.floor(rd.min()), np.ceil(rd.max()) + SHELL, SHELL)
    cg, _ = np.histogram(rd, edges, weights=wd)
    cr, _ = np.histogram(rr, edges, weights=wr)
    vol = 4 * np.pi / 3 * (edges[1:] ** 3 - edges[:-1] ** 3)
    fsky = cr / (vol * (cr.sum() / (4 * np.pi / 3 * (edges[-1] ** 3 - edges[0] ** 3))))
    n = np.divide(cg, np.where(cr > 0, cr, np.nan))      # galaxies per random => prop. to number density
    interior = (cr > 0.02 * cr.max())
    ni = n[interior]
    ratio = float(np.nanmax(ni) / np.nanmin(ni))
    res = {"cap": a.cap, "z_range": [ZMIN, ZMAX], "N_gal": int(len(rd)), "N_ran": int(len(rr)),
           "sum_w_gal": float(wd.sum()), "sum_w_ran": float(wr.sum()),
           "shell_Mpch": SHELL, "r_min_Mpch": float(rd.min()), "r_max_Mpch": float(rd.max()),
           "rule": "interior shells (weighted random count > 2% of max); max/min of the galaxy-per-random "
                   "ratio must be <= 1.5 for the sample to be treated as volume-limited",
           "n_interior_shells": int(interior.sum()),
           "ratio_max_over_min": ratio, "passed": bool(ratio <= 1.5),
           "n_profile_r_Mpch": ((edges[:-1] + edges[1:]) / 2)[interior].tolist(),
           "n_profile_gal_per_random": ni.tolist()}
    json.dump(res, open(os.path.join(OUT, "desi_gate1_%s.json" % a.cap), "w"), indent=1)
    log("G1 %s: N=%d shells=%d ratio=%.3f passed=%s" % (a.cap, len(rd), interior.sum(), ratio, res["passed"]))


# --------------------------------------------------------------------- data
def cmd_data(a):
    xd, wd, rd, xr, wr, rr = _load(a.cap)
    g = Grid(xd, wd, xr, wr)
    d = g.delta(xd, wd)
    pts, sd, vol = g.peaks(d)
    st = desi_statistics(pts)
    st.update(sigma_delta=sd, valid_volume_Mpch3=vol, alpha=g.alpha,
              grid_shape=list(g.shape), n_valid_cells=int(g.valid.sum()),
              n_eroded_cells=int(g.eroded.sum()))
    res = {"cap": a.cap, "N_gal": int(len(xd)), "N_ran": int(len(xr)),
           "weights": "WEIGHT * WEIGHT_FKP (galaxies and randoms)", "primary": st}
    np.savez_compressed(os.path.join(OUT, "desi_peaks_%s.npz" % a.cap), pts=pts)
    json.dump(res, open(os.path.join(OUT, "desi_data_%s.json" % a.cap), "w"), indent=1)
    log("data %s: N=%d S1=%s S3=%s" % (a.cap, st["S2_count"], st["S1_iqr_over_median"], st["S3_Q6_site_mean"]))


# --------------------------------------------------------------------- controls
def shell_match(rr, wr, rd, wd, rng):
    """Draw randoms with the data's per-shell weighted counts (X2's S3 rule),
    so the control has the data's exact radial selection."""
    edges = np.arange(np.floor(min(rd.min(), rr.min())), np.ceil(max(rd.max(), rr.max())) + SHELL, SHELL)
    sd = np.clip(np.searchsorted(edges, rd, "right") - 1, 0, len(edges) - 2)
    sr = np.clip(np.searchsorted(edges, rr, "right") - 1, 0, len(edges) - 2)
    need = np.bincount(sd, minlength=len(edges) - 1)
    order = np.argsort(sr, kind="stable")
    starts = np.searchsorted(sr[order], np.arange(len(edges) - 1))
    ends = np.searchsorted(sr[order], np.arange(len(edges) - 1), side="right")
    out = []
    for s in range(len(edges) - 1):
        k = need[s]
        if k == 0 or ends[s] <= starts[s]:
            continue
        pool = order[starts[s]:ends[s]]
        out.append(rng.choice(pool, size=min(k, pool.size), replace=False))
    return np.concatenate(out)


def cmd_control(a):
    xd, wd, rd, xr, wr, rr = _load(a.cap)
    g = Grid(xd, wd, xr, wr)
    rows = []
    for k in range(a.n):
        rng = np.random.default_rng(7100000 + k)
        idx = shell_match(rr, wr, rd, wd, rng)
        pts, sd, _ = g.peaks(g.delta(xr[idx], wr[idx]))
        st = desi_statistics(pts)
        st["sigma_delta"] = sd
        rows.append({kk: st[kk] for kk in ("S1_iqr_over_median", "S2_count", "S3_Q6_site_mean", "sigma_delta")})
        if k % 10 == 0:
            log("N2 %d N=%d S1=%s" % (k, st["S2_count"], st["S1_iqr_over_median"]))
    json.dump({"cap": a.cap, "n": a.n, "seeds": "7100000+k",
               "note": "OFFICIAL RANDOMS, shell-matched, IDENTICAL detector. This destroys ALL clustering, "
                       "so it is the analogue of Re6Zr's uniform control: reported as a CONTROL with NO "
                       "p-value against the data (registration N2 / limit L4). It is NOT a substitute for "
                       "the clustering-matched null.",
               "rows": rows}, open(os.path.join(OUT, "desi_control_%s.json" % a.cap), "w"))


# --------------------------------------------------------------------- injection
def cmd_inject(a):
    xd, wd, rd, xr, wr, rr = _load(a.cap)
    g = Grid(xd, wd, xr, wr)
    base_pts, base_sd, vol = g.peaks(g.delta(xd, wd))
    amps = [0.0, 0.5, 1.0, 2.0]
    dens = [1, 3, 10]                       # per 1e6 (Mpc/h)^3
    ijk = np.argwhere(g.eroded)
    cells = []
    for ia, A in enumerate(amps):
        for idn, D in enumerate(dens):
            if A == 0.0 and idn > 0:
                continue
            rows = []
            for k in range(a.n):
                rng = np.random.default_rng(7300000 + 10000 * ia + 1000 * idn + k)
                ridx = shell_match(rr, wr, rd, wd, rng)
                xb, wb = xr[ridx], wr[ridx]
                Ninj = max(1, int(round(D * 1e-6 * vol)))
                if A > 0:
                    # jittered cubic lattice inside the eroded region (registered)
                    side = max(1, int(round(Ninj ** (1 / 3))))
                    sub = ijk[rng.choice(len(ijk), size=min(Ninj, len(ijk)), replace=False)]
                    cen = g.origin + sub * L.L_CELL
                    lat = (vol / max(Ninj, 1)) ** (1 / 3)
                    cen = cen + rng.normal(0, 0.25 * lat, cen.shape)
                    # each injected core adds galaxies with a Gaussian profile of width sigma_G
                    per = max(1, int(round(A * base_sd * g.alphaR[g.eroded].mean() / wb.mean())))
                    ex = np.repeat(cen, per, axis=0) + rng.normal(0, L.SIGMA_G, (len(cen) * per, 3))
                    xb = np.vstack([xb, ex.astype(np.float32)])
                    wb = np.concatenate([wb, np.full(len(ex), wb.mean(), np.float32)])
                pts, sd, _ = g.peaks(g.delta(xb, wb))
                st = desi_statistics(pts)
                rows.append({kk: st[kk] for kk in ("S1_iqr_over_median", "S2_count", "S3_Q6_site_mean")})
            cells.append(dict(amp_over_sigma_delta=A, density_per_1e6Mpc3=D, n_inj=Ninj, n_sims=a.n,
                              seed_rule="7300000+10000*%d+1000*%d+k" % (ia, idn), rows=rows))
            log("inject A=%.2f D=%d (N_inj=%d) done" % (A, D, Ninj))
    json.dump({"cap": a.cap, "base": "official randoms, shell-matched (clustering-free)",
               "valid_volume_Mpch3": vol, "sigma_delta_data": base_sd, "cells": cells},
              open(os.path.join(OUT, "desi_inject_%s.json" % a.cap), "w"))


# --------------------------------------------------------------------- analyze
def cmd_analyze(a):
    dat = json.load(open(os.path.join(OUT, "desi_data_%s.json" % a.cap)))
    res = {"cap": a.cap, "N_gal": dat["N_gal"], "weights": dat["weights"], "data": dat["primary"]}
    res["G1"] = json.load(open(os.path.join(OUT, "desi_gate1_%s.json" % a.cap)))
    keys = ("S1_iqr_over_median", "S2_count", "S3_Q6_site_mean")
    get = lambda rs, k: np.array([r[k] for r in rs], float)

    cf = os.path.join(OUT, "desi_control_%s.json" % a.cap)
    if os.path.exists(cf):
        c = json.load(open(cf))
        res["N2_randoms_control_NO_P_VALUE"] = {
            "note": c["note"], "seeds": c["seeds"], "n": c["n"],
            **{k: {"mean": float(np.nanmean(get(c["rows"], k))),
                   "std": float(np.nanstd(get(c["rows"], k), ddof=1)),
                   "p2p5": float(np.nanpercentile(get(c["rows"], k), 2.5)),
                   "p97p5": float(np.nanpercentile(get(c["rows"], k), 97.5))} for k in keys}}

    lf = os.path.join(OUT, "desi_lognormal_%s.json" % a.cap)
    if os.path.exists(lf):
        m = json.load(open(lf))
        res["N1_lognormal"] = {"status": m["status"], "deviation": m["deviation"], "gate_G2": m["gate_G2"]}
        if m["gate_G2"]["passed"]:
            res["N1_lognormal"]["statistics"] = {
                k: {"data": dat["primary"][k], "mock_mean": float(np.nanmean(get(m["rows"], k))),
                    "mock_std": float(np.nanstd(get(m["rows"], k), ddof=1)),
                    "empirical_rank_p_two_sided": L.empirical_two_sided_p(get(m["rows"], k), dat["primary"][k])}
                for k in keys}
        else:
            res["N1_lognormal"]["verdict"] = ("G2 FAILED -> Test B is INCONCLUSIVE against a clustering-matched "
                                              "null, by the registered rule (the same rule that made round 2's "
                                              "X2 inconclusive at 2.88 sigma). No p-value is reported.")
    else:
        res["N1_lognormal"] = {"status": "NOT ATTEMPTED",
                               "verdict": "registered budget clause: no DESI-footprint mock catalogues exist on "
                                          "disk and the CAMB port was not completed -> Test B is INCONCLUSIVE "
                                          "with respect to a clustering-matched null. The randoms are NOT "
                                          "substituted for it."}

    inf = os.path.join(OUT, "desi_inject_%s.json" % a.cap)
    if os.path.exists(inf) and os.path.exists(cf):
        inj = json.load(open(inf))
        c = json.load(open(cf))
        lo = {k: float(np.nanpercentile(get(c["rows"], k), 2.5)) for k in keys}
        hi = {k: float(np.nanpercentile(get(c["rows"], k), 97.5)) for k in keys}
        grid = []
        for cc in inj["cells"]:
            g = {kk: cc[kk] for kk in ("amp_over_sigma_delta", "density_per_1e6Mpc3", "n_inj", "n_sims")}
            for k in keys:
                x = get(cc["rows"], k)
                g["rate_" + k] = float(np.mean((x < lo[k]) | (x > hi[k])))
                g["mean_" + k] = float(np.nanmean(x))
            grid.append(g)
        res["injection"] = {"criterion": "outside the two-sided central 95% interval of the randoms control "
                                         "ensemble (the reference ensemble actually available)",
                            "reference_interval": {k: [lo[k], hi[k]] for k in keys}, "grid": grid}
        z = [g for g in grid if g["amp_over_sigma_delta"] == 0.0]
        if z:
            res["injection"]["false_positive_rate_amp_zero"] = {k: z[0]["rate_" + k] for k in keys}
        det = [g for g in grid if g["amp_over_sigma_delta"] > 0 and g["rate_S1_iqr_over_median"] >= 0.95]
        res["injection"]["smallest_detected_at_95pc_S1"] = (
            min(det, key=lambda g: (g["amp_over_sigma_delta"], g["density_per_1e6Mpc3"])) if det else
            "NONE of the registered (amplitude, density) cells reached 95% detection on S1")
    json.dump(res, open(os.path.join(OUT, "desi_analysis_%s.json" % a.cap), "w"), indent=1)
    print(json.dumps({"N": dat["primary"]["S2_count"], "S1": dat["primary"]["S1_iqr_over_median"],
                      "S3": dat["primary"]["S3_Q6_site_mean"],
                      "N1": res["N1_lognormal"]["status"]}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["gate1", "data", "control", "inject", "analyze"])
    ap.add_argument("--cap", default="NGC", choices=["NGC", "SGC"])
    ap.add_argument("--n", type=int, default=100)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    globals()["cmd_" + a.cmd](a)

# NOTE on the registered null N1 (lognormal mocks).  No `lognormal` sub-command
# exists in this file and none was run.  Reason, recorded here rather than left
# implicit: no DESI-footprint mock catalogues exist on disk (audit/reverse_zero/
# E5-cosmic-web-tda-scaled/ holds only Betti-curve stacks for the SDSS DR17
# window, and audit/reverse_zero_r2/X2-cosmic-web/x2_lib.py is wired to that
# window's grid, healpix footprint and n(z)), and porting the CAMB -> lognormal
# -> bias-fit -> RSD -> Poisson pipeline to the DESI BGS footprint was outside
# this run's compute budget.  registration.json TEST_B.null_and_controls.
# N1_lognormal.budget_clause therefore applies: N1 is reported as NOT ATTEMPTED
# and Test B is INCONCLUSIVE with respect to a clustering-matched null.  The
# official randoms are NOT substituted for it.
