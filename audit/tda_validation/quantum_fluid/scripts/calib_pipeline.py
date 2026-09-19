"""Calibration of the two pipeline paths on synthetic inputs whose answer is
known by construction (run BEFORE writing expectations.json, so that the
expected TDA signatures are measured values, not hand arithmetic).

(1) alpha path (cosmic_web_tda_scaled.alpha_persistence) on
    - a perfect triangular lattice, spacing a = 1 (rows 30 x 30),
    - the same lattice with Gaussian jitter sigma = 0.05a, 0.10a (seed 1, 2),
    - a uniform random (Poisson, fixed N) point set with the SAME number of
      points in the SAME bounding box (seeds 11..15).
(2) lower-star path (cmb_tda.betti_curves_from_topology)
    - torus topology Euler check V - E + F = 0 for L = 8, 32, 64,
    - field on a 64x64 torus = 1 - sum of K narrow Gaussian dips at known
      separated positions: b0 at a threshold below the background must equal K;
      b0 of a constant+tiny-noise field at a threshold above everything must be 1
      and b1 must be 2 (torus H1).

Command:
  prlimit --as=8589934592 -- .venv-tda/bin/python calib_pipeline.py
Output: ../results/calib_pipeline.json
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qf_common as qc  # noqa: E402


def tri_lattice(nx, ny, a=1.0):
    pts = []
    for j in range(ny):
        for i in range(nx):
            pts.append([a * (i + 0.5 * (j % 2)), a * j * np.sqrt(3) / 2])
    return np.array(pts)


def alpha_summary(pts, label, max_r):
    web = qc.web()
    info, by_dim = web.alpha_persistence(pts, max_r ** 2, label, None)
    d0 = qc.h0_finite_deaths(by_dim)
    h1 = np.array(by_dim[1]) if by_dim[1] else np.empty((0, 2))
    tb0 = web.top_bars(by_dim, 0, k=3, r_trunc=max_r)
    out = {
        "label": label, "n_points": int(len(pts)),
        "h0_death_stats": qc.spread_stats(d0),
        "h1_n_bars": int(len(h1)),
        "h1_birth_stats": qc.spread_stats(h1[:, 0]) if len(h1) else None,
        "h1_death_stats": qc.spread_stats(h1[:, 1]) if len(h1) else None,
        "h1_persistence_stats": qc.spread_stats(h1[:, 1] - h1[:, 0]) if len(h1) else None,
        "h1_median_death_over_median_h0_death": float(np.median(h1[:, 1]) / np.median(d0)) if len(h1) else None,
        "top3_h0_bars": tb0.tolist(),
        "betti_at_truncation": info["betti_numbers_at_truncation"],
    }
    return out


def main():
    res = {"alpha_path": {}, "lower_star_path": {}}
    nx = ny = 30
    base = tri_lattice(nx, ny)
    res["alpha_path"]["perfect_triangular_a1"] = alpha_summary(base, "tri_a1", 3.0)
    for s, seed in [(0.05, 1), (0.10, 2)]:
        rng = np.random.default_rng(seed)
        res["alpha_path"][f"jitter_{s}_seed{seed}"] = alpha_summary(base + rng.normal(0, s, base.shape), f"tri_jit{s}", 3.0)
    lo, hi = base.min(0), base.max(0)
    pois = []
    for seed in range(11, 16):
        rng = np.random.default_rng(seed)
        p = lo + rng.random(base.shape) * (hi - lo)
        pois.append(alpha_summary(p, f"poisson_seed{seed}", 3.0))
    res["alpha_path"]["poisson_same_N_same_box_seeds11_15"] = pois
    res["alpha_path"]["poisson_iqr_over_median_mean"] = float(np.mean([p["h0_death_stats"]["iqr_over_median"] for p in pois]))
    res["alpha_path"]["poisson_median_h0_death_mean"] = float(np.mean([p["h0_death_stats"]["median"] for p in pois]))
    res["alpha_path"]["reference_values"] = {"a_over_2": 0.5, "a_over_sqrt3": float(1 / np.sqrt(3)), "ratio_2_over_sqrt3": float(2 / np.sqrt(3))}

    cmb = qc.cmb()
    eul = {}
    for L in (8, 32, 64):
        u, e, t = qc.torus_topology(L)
        eul[str(L)] = {"V": int(u.size), "E": int(e.shape[0]), "F": int(t.shape[0]), "chi": int(u.size - e.shape[0] + t.shape[0])}
    res["lower_star_path"]["torus_euler"] = eul
    L = 64
    u, e, t = qc.torus_topology(L)
    yy, xx = np.mgrid[0:L, 0:L]
    rng = np.random.default_rng(3)
    K = 7
    centers = [(8, 8), (8, 40), (30, 20), (50, 50), (40, 5), (20, 55), (55, 30)]
    f = np.ones((L, L))
    for (cx, cy) in centers:
        dx = qc.wrap(2 * np.pi * (xx - cx) / L) * L / (2 * np.pi)
        dy = qc.wrap(2 * np.pi * (yy - cy) / L) * L / (2 * np.pi)
        f -= np.exp(-(dx ** 2 + dy ** 2) / (2 * 1.5 ** 2))
    f += 1e-3 * rng.normal(size=f.shape)
    flat = f.ravel()
    sigma = flat.std()
    thr = np.array([0.5, 2.0]) / sigma  # absolute thresholds 0.5 and 2.0 in field units
    b0, b1, chi = cmb.betti_curves_from_topology(flat, u, e, t, thr, sublevel=True)
    res["lower_star_path"]["gaussian_dips_torus64"] = {
        "K_dips": K, "abs_thresholds": [0.5, 2.0], "b0": b0.tolist(), "b1": b1.tolist(), "chi": chi.tolist(),
        "expected": {"b0_at_0.5": K, "b0_at_2.0": 1, "b1_at_2.0": 2, "chi_at_2.0": -1}}
    with open(os.path.join(qc.RESULTS, "calib_pipeline.json"), "w") as fh:
        json.dump(res, fh, indent=1)
    print(json.dumps(res, indent=1)[:4000])


if __name__ == "__main__":
    main()
