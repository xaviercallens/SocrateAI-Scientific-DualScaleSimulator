#!/usr/bin/env python3
"""TEST A of registration.json (CV): the Re6Zr defect-core pipeline on the CMB.

Sub-commands (all seeded; seeds are the registered ones):
  null    --which wmap|smica --n 500      500 spectrum-matched Gaussian sims, seeds 2000000+k
  data    --which wmap|smica              the real map through the IDENTICAL detector
  control --which wmap|smica              C1 uniform (seed 201+k) and C2 value-shuffle (6100000+k)
  inject  --which wmap                    the amplitude x density injection grid
  analyze --which wmap|smica              gate, p-values, decision

Every command:
  timeout <s> prlimit --as=8589934592 -- <venv-tda python> cv_cmb.py <sub> ...
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time

import numpy as np
import healpy as hp

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cv_lib as L  # noqa: E402

OUT = os.path.join(L.CV, "results")


def log(m):
    print("[%s] %s" % (time.strftime("%H:%M:%S"), m), file=sys.stderr, flush=True)


def _setup(which):
    m, mask = L.load_data(which)
    tab, _ = L.disk_table()
    er = L.eroded_mask(mask, tab)
    return m, mask, tab, er


def _one(t_raw, mask, tab, er, nu=L.NU_PRIMARY, sign=0):
    """prep -> smooth -> detect -> statistics.  The SAME path for data, nulls,
    injections and controls."""
    tp = L.prep_map(t_raw, mask)
    pix, sd, ts = L.detect_cmb(tp, mask, tab, er, nu=nu, sign=sign)
    st = L.cmb_point_statistics(pix)
    st["sigma_Ts"] = sd
    st["var_Ts_eroded"] = float(ts[er].var())
    return st, pix, ts, tp


# --------------------------------------------------------------------- null
def cmd_null(a):
    cl_in, fsky = L.load_cl_in(a.which)
    m, mask, tab, er = _setup(a.which)
    rows = []
    for k in range(a.n):
        t = L.make_sim(cl_in, 2000000 + k)
        st, _, _, _ = _one(t, mask, tab, er)
        st["seed"] = 2000000 + k
        rows.append(st)
        if k % 50 == 0:
            log("null %s %d/%d S1=%.4f N=%d S3=%.4f" % (a.which, k, a.n, st["S1_iqr_over_median"],
                                                        st["S2_count"], st["S3_psi6_site_mean"]))
    json.dump({"which": a.which, "n": a.n, "seeds": "2000000+k, k=0..%d" % (a.n - 1),
               "fsky_mask": float(mask.mean()), "fsky_eroded": float(er.mean()),
               "rows": rows}, open(os.path.join(OUT, "cmb_null_%s.json" % a.which), "w"))
    log("null %s done" % a.which)


# --------------------------------------------------------------------- data
def cmd_data(a):
    m, mask, tab, er = _setup(a.which)
    res = {"which": a.which, "fsky_mask": float(mask.mean()), "fsky_eroded": float(er.mean()),
           "map_file": L.MAPS[a.which]["map"], "mask_file": L.MAPS[a.which]["mask"]}
    st, pix, ts, tp = _one(m, mask, tab, er)
    res["primary"] = st
    np.savez_compressed(os.path.join(OUT, "cmb_peaks_%s.npz" % a.which), pix=pix)
    # registered secondary diagnostics (outside the multiplicity family)
    for tag, kw in (("nu1.5", dict(nu=L.NU_SECONDARY)), ("maxima_only", dict(sign=+1)),
                    ("minima_only", dict(sign=-1))):
        s2, _, _, _ = _one(m, mask, tab, er, **kw)
        res.setdefault("secondary", {})[tag] = s2
    json.dump(res, open(os.path.join(OUT, "cmb_data_%s.json" % a.which), "w"), indent=1)
    log("data %s: N=%d S1=%.4f S3=%.4f" % (a.which, st["S2_count"], st["S1_iqr_over_median"],
                                           st["S3_psi6_site_mean"]))


# --------------------------------------------------------------------- controls
def cmd_control(a):
    m, mask, tab, er = _setup(a.which)
    dat = json.load(open(os.path.join(OUT, "cmb_data_%s.json" % a.which)))
    N = dat["primary"]["S2_count"]
    valid = np.where(er)[0]
    c1 = []
    for k in range(a.n):
        # C1: Re6Zr's uniform control, seed rule 201+k, NOT passed through the
        # detector.  Pipeline sanity check only -- NO p-value (registration L4).
        rng = np.random.default_rng(201 + k)
        pix = rng.choice(valid, size=N, replace=False)
        st = L.cmb_point_statistics(pix)
        c1.append(st)
    c2 = []
    tp = L.prep_map(m, mask)
    ts = L.smooth_map(tp, mask)
    for k in range(a.n):
        # C2: value-shuffle of T_s inside the eroded region, FULL detector rerun.
        rng = np.random.default_rng(6100000 + k)
        sh = ts.copy()
        sh[er] = rng.permutation(ts[er])
        pix, sd, _ = L.detect_cmb(None, mask, tab, er, ts_precomputed=sh)
        st = L.cmb_point_statistics(pix)
        st["sigma_Ts"] = sd
        c2.append(st)
        if k % 25 == 0:
            log("C2 %d N=%d S1=%s" % (k, st["S2_count"], st["S1_iqr_over_median"]))
    json.dump({"which": a.which, "N_matched": N,
               "C1_uniform_no_p_value": {"seeds": "201+k", "n": a.n, "rows": c1,
                                         "note": "NOT passed through the detector; the contrast with the "
                                                 "data is dominated by the R_disk hard-core exclusion, which "
                                                 "the detector imposes by construction (registration L4). "
                                                 "Reported as a pipeline sanity check with NO p-value."},
               "C2_value_shuffle": {"seeds": "6100000+k", "n": a.n, "rows": c2,
                                    "note": "values of T_s permuted inside the eroded region, FULL detector rerun"}},
              open(os.path.join(OUT, "cmb_controls_%s.json" % a.which), "w"))
    log("controls %s done" % a.which)


# --------------------------------------------------------------------- injection
def fib_sphere(n, rng, jitter_frac=0.25):
    """Jittered geodesic grid (registered primary placement)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    gold = np.pi * (1 + 5 ** 0.5)
    th = gold * i
    u = np.stack([np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)], 1)
    d_nn = np.sqrt(4 * np.pi / n)             # rad, mean nearest-neighbour spacing
    for j in range(n):
        ref = np.array([0.0, 0.0, 1.0]) if abs(u[j, 2]) < 0.99 else np.array([1.0, 0.0, 0.0])
        e1 = ref - (ref @ u[j]) * u[j]
        e1 /= np.linalg.norm(e1)
        e2 = np.cross(u[j], e1)
        ang = rng.uniform(0, 2 * np.pi)
        step = rng.normal(0, jitter_frac * d_nn)
        u[j] = u[j] * np.cos(step) + (np.cos(ang) * e1 + np.sin(ang) * e2) * np.sin(step)
    return u / np.linalg.norm(u, axis=1, keepdims=True)


def poisson_sphere(n, rng):
    v = rng.normal(size=(n, 3))
    return v / np.linalg.norm(v, axis=1, keepdims=True)


def inject(t, centres, amp, theta0_deg=L.SIGMA_S_DEG, rng=None):
    """T += s*A*sigma_T*exp(-theta^2/(2 theta0^2)), s = +-1 equiprobable."""
    if amp == 0.0 or len(centres) == 0:
        return t
    npix = hp.nside2npix(L.NSIDE)
    vec = np.array(hp.pix2vec(L.NSIDE, np.arange(npix))).T
    th0 = np.radians(theta0_deg)
    out = t.copy()
    R = 4 * th0
    for c in centres:
        d = hp.query_disc(L.NSIDE, c, R)
        cosd = np.clip(vec[d] @ c, -1, 1)
        ang = np.arccos(cosd)
        s = 1.0 if rng.random() < 0.5 else -1.0
        out[d] += s * amp * np.exp(-ang ** 2 / (2 * th0 ** 2))
    return out


def cmd_inject(a):
    cl_in, fsky = L.load_cl_in(a.which)
    m, mask, tab, er = _setup(a.which)
    sigma_T = float(L.prep_map(m, mask)[mask > 0].std())
    amps = [0.0, 0.25, 0.5, 1.0, 2.0]
    dens = [100, 300, 1000]
    cells = []
    for ia, A in enumerate(amps):
        for idn, Nj in enumerate(dens):
            if A == 0.0 and idn > 0:
                continue                      # A=0 run ONCE (registered)
            for placement in (["grid"] if A == 0.0 else ["grid", "poisson"]):
                rows = []
                for k in range(a.n):
                    gs = 4100000 + 10000 * ia + 1000 * idn + k
                    ps = 5100000 + 10000 * ia + 1000 * idn + k + (500000 if placement == "poisson" else 0)
                    t = L.make_sim(cl_in, gs)
                    rng = np.random.default_rng(ps)
                    cen = (fib_sphere(Nj, rng) if placement == "grid" else poisson_sphere(Nj, rng))
                    t = inject(t, cen, A * sigma_T, rng=rng)
                    st, _, _, _ = _one(t, mask, tab, er)
                    rows.append({kk: st[kk] for kk in ("S1_iqr_over_median", "S2_count", "S3_psi6_site_mean")})
                cells.append(dict(amp_over_sigmaT=A, n_inj_full_sky=Nj, placement=placement,
                                  n_sims=a.n, gaussian_seed_rule="4100000+10000*%d+1000*%d+k" % (ia, idn),
                                  placement_seed_rule="5100000+10000*%d+1000*%d+k%s" % (ia, idn,
                                                      "+500000" if placement == "poisson" else ""),
                                  rows=rows))
                log("inject A=%.2f N=%d %s done" % (A, Nj, placement))
    json.dump({"which": a.which, "sigma_T_mK_or_K": sigma_T, "theta0_deg": L.SIGMA_S_DEG,
               "profile": "T += s*A*sigma_T*exp(-theta^2/(2 theta0^2)), s=+-1 equiprobable",
               "cells": cells}, open(os.path.join(OUT, "cmb_inject_%s.json" % a.which), "w"))


# --------------------------------------------------------------------- analyze
def cmd_analyze(a):
    nul = json.load(open(os.path.join(OUT, "cmb_null_%s.json" % a.which)))
    dat = json.load(open(os.path.join(OUT, "cmb_data_%s.json" % a.which)))
    rows = nul["rows"]
    get = lambda rs, k: np.array([r[k] for r in rs], float)
    res = {"which": a.which, "n_null": len(rows), "null_seeds": nul["seeds"],
           "fsky_mask": nul["fsky_mask"], "fsky_eroded": nul["fsky_eroded"]}

    # ---- fail-closed gate on var(T_s), the quantity this test actually reads
    v = get(rows, "var_Ts_eroded")
    z = (dat["primary"]["var_Ts_eroded"] - v.mean()) / v.std(ddof=1)
    passed = bool(abs(z) <= L.GATE_THRESH)
    res["gate"] = {"quantity": "var(T_s) over the eroded valid region", "z": float(z),
                   "threshold": L.GATE_THRESH, "passed": passed,
                   "data": dat["primary"]["var_Ts_eroded"], "sim_mean": float(v.mean()),
                   "sim_std": float(v.std(ddof=1)),
                   "decision": "PASS: p-values may be computed" if passed else
                               "FAIL: null mis-specified for T_s; NO p-values; map NOT TESTED by CV"}
    res["x1_gate_for_reference"] = json.load(open(os.path.join(L.X1, a.which, "gate.json")))["max_abs_z"]

    res["statistics"] = {}
    for key, name in (("S1_iqr_over_median", "S1 H0 death-radius IQR/median"),
                      ("S2_count", "S2 detected core count"),
                      ("S3_psi6_site_mean", "S3 site-averaged psi6")):
        nv = get(rows, key)
        dv = dat["primary"][key]
        e = {"name": name, "data": dv, "null_mean": float(np.nanmean(nv)),
             "null_std": float(np.nanstd(nv, ddof=1)),
             "null_p2p5": float(np.nanpercentile(nv, 2.5)), "null_p97p5": float(np.nanpercentile(nv, 97.5)),
             "z": float((dv - np.nanmean(nv)) / np.nanstd(nv, ddof=1)),
             "empirical_rank_p_two_sided": L.empirical_two_sided_p(nv, dv)}
        e["p_value_convention"] = "min(1, 2*min((1+#{null<=d})/(N+1), (1+#{null>=d})/(N+1)))"
        if not passed:
            e["empirical_rank_p_two_sided"] = None
            e["withheld"] = "gate failed (fail closed)"
        res["statistics"][key] = e

    res["absolute_floor_data"] = dat["primary"].get("absolute_floor")
    res["absolute_floor_null_mean"] = {
        k: float(np.mean([r["absolute_floor"][k] for r in rows if "absolute_floor" in r]))
        for k in ("median_death_deg", "twice_median_death_deg", "frac_bars_below_pixel_scale", "max_death_deg")}
    res["truncation_data"] = dat["primary"].get("truncation")
    res["truncation_null_mean_frac_at_trunc"] = float(np.mean(
        [r["truncation"]["frac_h0_deaths_at_or_above_r_trunc"] for r in rows if "truncation" in r]))
    res["truncation_inert"] = bool(res["truncation_null_mean_frac_at_trunc"] <= 0.01 and
                                   (res["truncation_data"] or {}).get("frac_h0_deaths_at_or_above_r_trunc", 1) <= 0.01)
    res["secondary_data"] = dat.get("secondary")

    cf = os.path.join(OUT, "cmb_controls_%s.json" % a.which)
    if os.path.exists(cf):
        c = json.load(open(cf))
        for tag in ("C1_uniform_no_p_value", "C2_value_shuffle"):
            rs = c[tag]["rows"]
            res.setdefault("controls", {})[tag] = {
                "note": c[tag]["note"], "seeds": c[tag]["seeds"], "n": len(rs),
                **{k: {"mean": float(np.nanmean(get(rs, k))), "std": float(np.nanstd(get(rs, k), ddof=1))}
                   for k in ("S1_iqr_over_median", "S2_count", "S3_psi6_site_mean")}}
        # C2 carries a p-value: it IS passed through the detector
        for k in ("S1_iqr_over_median", "S2_count", "S3_psi6_site_mean"):
            res["controls"]["C2_value_shuffle"][k]["empirical_rank_p_vs_null"] = L.empirical_two_sided_p(
                get(rows, k), float(np.nanmean(get(c["C2_value_shuffle"]["rows"], k))))

    inf = os.path.join(OUT, "cmb_inject_%s.json" % a.which)
    if os.path.exists(inf):
        inj = json.load(open(inf))
        lo = {k: float(np.nanpercentile(get(rows, k), 2.5)) for k in
              ("S1_iqr_over_median", "S2_count", "S3_psi6_site_mean")}
        hi = {k: float(np.nanpercentile(get(rows, k), 97.5)) for k in lo}
        grid = []
        for c in inj["cells"]:
            rs = c["rows"]
            g = {kk: c[kk] for kk in ("amp_over_sigmaT", "n_inj_full_sky", "placement", "n_sims")}
            for k in lo:
                x = get(rs, k)
                g["rate_" + k] = float(np.mean((x < lo[k]) | (x > hi[k])))
                g["mean_" + k] = float(np.nanmean(x))
            g["rate_any"] = float(np.mean([any((r[k] < lo[k]) or (r[k] > hi[k]) for k in lo) for r in rs]))
            grid.append(g)
        res["injection"] = {"criterion": "outside the two-sided central 95% interval of the 500 clean nulls",
                            "null_interval": {k: [lo[k], hi[k]] for k in lo}, "grid": grid}
        fp = [g for g in grid if g["amp_over_sigmaT"] == 0.0]
        if fp:
            res["injection"]["false_positive_rate_amp_zero"] = {
                k: fp[0]["rate_" + k] for k in lo} | {"any": fp[0]["rate_any"], "n_sims": fp[0]["n_sims"]}
        det = [g for g in grid if g["amp_over_sigmaT"] > 0 and g["placement"] == "grid"
               and g["rate_S1_iqr_over_median"] >= 0.95]
        res["injection"]["smallest_detected_at_95pc_S1_grid_placement"] = (
            min(det, key=lambda g: (g["amp_over_sigmaT"], g["n_inj_full_sky"])) if det else
            "NONE of the registered (amplitude, density) cells reached 95% detection on S1")
    json.dump(res, open(os.path.join(OUT, "cmb_analysis_%s.json" % a.which), "w"), indent=1)
    print(json.dumps({"gate": res["gate"]["passed"], "gate_z": res["gate"]["z"],
                      **{k: (res["statistics"][k]["data"], res["statistics"][k]["empirical_rank_p_two_sided"])
                         for k in res["statistics"]}}, indent=1))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["null", "data", "control", "inject", "analyze"])
    ap.add_argument("--which", default="wmap", choices=list(L.MAPS))
    ap.add_argument("--n", type=int, default=500)
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    globals()["cmd_" + a.cmd](a)
