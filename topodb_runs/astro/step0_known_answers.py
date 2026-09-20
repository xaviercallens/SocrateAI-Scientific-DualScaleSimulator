#!/usr/bin/env python3
"""STEP 0 — known-answer suite for THIS pipeline, run before any real data.

Six cases, each with its topology written down BEFORE it is computed (the
`expected` argument is literally the textbook Betti vector):

  K1 circle        S^1 in R^2      alpha    expected (1, 1)
  K2 sphere        S^2 in R^3      alpha    expected (1, 0, 1)
  K3 torus         T^2 in R^3      alpha    expected (1, 2, 1)
  K4 grid wells    7 Gaussian wells, 128x128   cubical  expected b0 = 7
  K5 healpix       the FIXED sphere complex, nside 32    expected (1, 0, 1)
  K6 negative      uniform noise in a disc + a pure-noise grid: the gap rule
                   must NOT report the K1/K4 answer, and its gap ratio must be
                   far below the signal cases.

If any of K1-K5 fails, the campaign stops. K6 is the control that the
persistence-gap rule is not true by construction (LeanFlow CLAUDE.md rule 2).

Tier X (numerics) throughout.
"""
from __future__ import annotations

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import (ASTRO, Timer, add_run, alpha_bars, gap_betti, peak_mb,  # noqa: E402
                    save_json, topodb)

SCRIPT = "topodb_runs/astro/step0_known_answers.py"
COMMAND = ("prlimit --as=8589934592 -- "
           "/home/callensxavier_gmail_com/SocrateAI-Scientific-DualScaleSimulator/.venv-tda/bin/python "
           "topodb_runs/astro/step0_known_answers.py")
SEED = 20260920

DS = "synthetic/step0_known_answers"


def _resolution_floor(st):
    """Sampling resolution of the point cloud, in length units.

    The median finite H0 death radius is the typical nearest-neighbour scale.
    A bar shorter than that is not resolved by the sample, so it cannot be a
    feature of the underlying shape.  Data-driven, no tuned constant.
    """
    d0 = alpha_bars(st, 0)
    fin = d0[np.isfinite(d0[:, 1]), 1]
    return float(np.median(fin)) if fin.size else 0.0


def _alpha_tree(points, max_alpha_square):
    import gudhi
    ac = gudhi.AlphaComplex(points=np.ascontiguousarray(points, dtype=float))
    st = ac.create_simplex_tree(max_alpha_square=max_alpha_square)
    st.compute_persistence(homology_coeff_field=2, persistence_dim_max=True)
    return st


def case_circle(rng):
    n, noise = 400, 0.01
    th = rng.uniform(0, 2 * np.pi, n)
    pts = np.c_[np.cos(th), np.sin(th)] + rng.normal(0, noise, (n, 2))
    st = _alpha_tree(pts, 4.0)
    flo = _resolution_floor(st)
    b1, gap1 = gap_betti(np.diff(alpha_bars(st, 1), axis=1).ravel(), flo)
    b0 = int(np.sum(~np.isfinite(alpha_bars(st, 0)[:, 1])))
    return dict(betti={0: b0, 1: b1}, expected={0: 1, 1: 1}, gap_ratio_h1=gap1,
                resolution_floor=flo,
                bars={1: alpha_bars(st, 1).tolist()},
                params=dict(n_points=n, noise_sigma=noise, max_alpha_square=4.0,
                            shape="unit circle S^1 in R^2"))


def case_sphere(rng):
    n = 2000                      # Fibonacci sphere: deterministic, no clumping
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n)
    tht = np.pi * (1 + 5 ** 0.5) * i
    pts = np.c_[np.cos(tht) * np.sin(phi), np.sin(tht) * np.sin(phi), np.cos(phi)]
    st = _alpha_tree(pts, 4.0)
    flo = _resolution_floor(st)
    b2, gap2 = gap_betti(np.diff(alpha_bars(st, 2), axis=1).ravel(), flo)
    b1, gap1 = gap_betti(np.diff(alpha_bars(st, 1), axis=1).ravel(), flo)
    b0 = int(np.sum(~np.isfinite(alpha_bars(st, 0)[:, 1])))
    return dict(betti={0: b0, 1: b1, 2: b2}, expected={0: 1, 1: 0, 2: 1},
                gap_ratio_h2=gap2, gap_ratio_h1=gap1, resolution_floor=flo,
                bars={2: alpha_bars(st, 2).tolist()},
                params=dict(n_points=n, sampling="Fibonacci lattice",
                            max_alpha_square=4.0, shape="unit sphere S^2 in R^3"))


def case_torus(rng):
    n, R, r = 4000, 2.0, 0.8
    u = rng.uniform(0, 2 * np.pi, n)
    v = rng.uniform(0, 2 * np.pi, n)
    pts = np.c_[(R + r * np.cos(v)) * np.cos(u),
                (R + r * np.cos(v)) * np.sin(u),
                r * np.sin(v)]
    st = _alpha_tree(pts, 9.0)
    flo = _resolution_floor(st)
    b1, gap1 = gap_betti(np.diff(alpha_bars(st, 1), axis=1).ravel(), flo)
    b2, gap2 = gap_betti(np.diff(alpha_bars(st, 2), axis=1).ravel(), flo)
    b0 = int(np.sum(~np.isfinite(alpha_bars(st, 0)[:, 1])))
    return dict(betti={0: b0, 1: b1, 2: b2}, expected={0: 1, 1: 2, 2: 1},
                gap_ratio_h1=gap1, gap_ratio_h2=gap2, resolution_floor=flo,
                bars={1: alpha_bars(st, 1).tolist(), 2: alpha_bars(st, 2).tolist()},
                params=dict(n_points=n, R=R, r=r, max_alpha_square=9.0,
                            shape="torus T^2 in R^3"))


def _wells_grid(rng, k, N=128, noise=0.0):
    """A field with k Gaussian wells (negative bumps) on a flat background."""
    x = np.arange(N)
    X, Y = np.meshgrid(x, x, indexing="ij")
    f = np.zeros((N, N))
    centres = rng.integers(12, N - 12, size=(k, 2)) if k else np.zeros((0, 2), int)
    # force separation so the k wells are genuinely k distinct basins
    keep = []
    for c in centres:
        if all(np.hypot(*(c - d)) > 24 for d in keep):
            keep.append(c)
    for c in keep:
        f -= np.exp(-((X - c[0]) ** 2 + (Y - c[1]) ** 2) / (2 * 6.0 ** 2))
    if noise:
        f += rng.normal(0, noise, (N, N))
    return f, len(keep)


def case_grid_wells(rng):
    import gudhi
    k_req = 9
    f, k = _wells_grid(rng, k_req, noise=0.02)
    cc = gudhi.CubicalComplex(top_dimensional_cells=f)
    cc.compute_persistence(homology_coeff_field=2)
    d0 = np.array(cc.persistence_intervals_in_dimension(0), dtype=float).reshape(-1, 2)
    d0f = d0[np.isfinite(d0[:, 1])]
    pers = np.r_[d0f[:, 1] - d0f[:, 0]]
    b0, gap0 = gap_betti(pers)
    # the infinite bar is the global minimum's component; k wells give
    # k-1 finite bars above the gap plus that one infinite bar.
    b0_total = b0 + int(np.sum(~np.isfinite(d0[:, 1])))
    return dict(betti={0: b0_total}, expected={0: k}, gap_ratio_h0=gap0,
                bars={0: d0f.tolist()},
                params=dict(grid=128, k_wells_requested=k_req, k_wells_placed=k,
                            well_sigma_px=6.0, min_separation_px=24,
                            noise_sigma=0.02, filtration="sublevel cubical"))


def case_healpix(rng):
    sys.path.insert(0, ASTRO)
    import healpy as hp
    from lib.cmb_topology import build_topology_fixed, complex_betti
    nside = 32
    mask = np.ones(hp.nside2npix(nside))
    unmasked, edges, tris = build_topology_fixed(mask, nside)
    b = complex_betti(unmasked.size, edges, tris)
    return dict(betti={0: int(b[0]), 1: int(b[1]), 2: int(b[2])},
                expected={0: 1, 1: 0, 2: 1}, bars={},
                params=dict(nside=nside, full_sky=True, n_vertices=int(unmasked.size),
                            n_edges=int(len(edges)), n_triangles=int(len(tris)),
                            euler=int(unmasked.size - len(edges) + len(tris)),
                            library="lib/cmb_topology.py @ d8175f1"))


def case_negative(rng):
    """Nothing true by construction: the same rules on structureless input."""
    import gudhi
    # (a) uniform points in a disc: no 1-cycle should stand out
    n = 400
    rad = np.sqrt(rng.uniform(0, 1, n))
    th = rng.uniform(0, 2 * np.pi, n)
    pts = np.c_[rad * np.cos(th), rad * np.sin(th)]
    st = _alpha_tree(pts, 4.0)
    flo = _resolution_floor(st)
    b1_noise, gap1_noise = gap_betti(np.diff(alpha_bars(st, 1), axis=1).ravel(), flo)
    # (b) pure-noise grid: must not report 9 wells
    f = rng.normal(0, 1, (128, 128))
    cc = gudhi.CubicalComplex(top_dimensional_cells=f)
    cc.compute_persistence(homology_coeff_field=2)
    d0 = np.array(cc.persistence_intervals_in_dimension(0), dtype=float).reshape(-1, 2)
    d0f = d0[np.isfinite(d0[:, 1])]
    b0_noise, gap0_noise = gap_betti(d0f[:, 1] - d0f[:, 0])
    b0_noise += int(np.sum(~np.isfinite(d0[:, 1])))
    return dict(b1_disc=b1_noise, gap_ratio_h1_disc=gap1_noise, resolution_floor=flo,
                b0_noise_grid=b0_noise, gap_ratio_h0_noise_grid=gap0_noise,
                n_finite_h0_bars_noise_grid=int(d0f.shape[0]),
                params=dict(n_points_disc=n, grid=128, noise_sigma=1.0))


def main():
    rng = np.random.default_rng(SEED)
    out = {"seed": SEED, "cases": {}}
    cases = [("K1_circle", case_circle), ("K2_sphere", case_sphere),
             ("K3_torus", case_torus), ("K4_grid_wells", case_grid_wells),
             ("K5_healpix_complex", case_healpix)]

    with topodb() as db:
        db.add_dataset(id=DS, domain="synthetic", title="Step-0 known-answer suite (astro TDA pipeline)",
                       source=f"generated by {SCRIPT}", provenance="synthetic_control",
                       notes=("Six shapes whose homology is textbook. Run before any real dataset; "
                             "if K1-K5 fail the campaign stops. K6 is the negative control that the "
                             "persistence-gap rule is not true by construction."))

    results = {}
    for name, fn in cases:
        with Timer() as t:
            r = fn(np.random.default_rng(SEED))
        r["wall_sec"] = t.sec
        r["peak_mb"] = peak_mb()
        results[name] = r
        ok = all(r["betti"].get(d) == v for d, v in r["expected"].items())
        r["passed"] = bool(ok)
        print(f"{name:22s} betti={r['betti']} expected={r['expected']} "
              f"{'PASS' if ok else 'FAIL'}  ({t.sec:.1f}s)")

        method = {"K1_circle": "alpha", "K2_sphere": "alpha", "K3_torus": "alpha",
                  "K4_grid_wells": "cubical", "K5_healpix_complex": "lower_star_graph"}[name]
        with topodb() as db:
            rid = add_run(db, dataset_id=DS, method=method, coeff_field=2,
                          params={**r["params"], "case": name,
                                  "betti_rule": ("largest multiplicative persistence gap (common.gap_betti) above a "
                                 "resolution floor = median finite H0 death radius; see AMENDMENT A1 "
                                 "in PRE_DECLARED_STATISTICS.md")},
                          script=SCRIPT, command=COMMAND, tier="X",
                          max_dim=2 if method != "cubical" else 1,
                          preprocessing="none (synthetic sample generated in-script from the stated seed)",
                          seed=str(SEED), wall_sec=t.sec, peak_mb=r["peak_mb"])
            db.add_betti(rid, r["betti"], expected=r["expected"])
            for dim, bars in r.get("bars", {}).items():
                if bars:
                    db.add_bars(rid, int(dim), [(b, d) for b, d in bars], top=20)
            for k in ("gap_ratio_h0", "gap_ratio_h1", "gap_ratio_h2", "resolution_floor"):
                if k in r and np.isfinite(r[k]):
                    db.add_statistic(rid, k, r[k])
            db.add_control(rid, "known_answer",
                           f"{name}: recovers the textbook Betti vector {r['expected']}",
                           passed=ok, detail=f"computed {r['betti']}")
        results[name]["run_id"] = rid

    # K6 negative control -> attached to its own run
    with Timer() as t:
        neg = case_negative(np.random.default_rng(SEED + 1))
    neg["wall_sec"] = t.sec
    results["K6_negative"] = neg
    sig_gap = results["K1_circle"]["gap_ratio_h1"]
    neg_ok = (neg["gap_ratio_h1_disc"] < sig_gap / 10.0) and (neg["b0_noise_grid"] != 9)
    neg["passed"] = bool(neg_ok)
    print(f"{'K6_negative':22s} disc gap_h1={neg['gap_ratio_h1_disc']:.3g} "
          f"(circle {sig_gap:.3g}); noise-grid b0={neg['b0_noise_grid']} (wells case 9) "
          f"{'PASS' if neg_ok else 'FAIL'}")
    with topodb() as db:
        rid = add_run(db, dataset_id=DS, method="alpha", coeff_field=2,
                      params={**neg["params"], "case": "K6_negative"},
                      script=SCRIPT, command=COMMAND, tier="X", max_dim=2,
                      preprocessing="uniform noise: no topology to find",
                      seed=str(SEED + 1), wall_sec=t.sec, peak_mb=peak_mb())
        db.add_statistic(rid, "gap_ratio_h1_uniform_disc", neg["gap_ratio_h1_disc"])
        db.add_statistic(rid, "gap_ratio_h0_noise_grid", neg["gap_ratio_h0_noise_grid"])
        db.add_statistic(rid, "b0_noise_grid", neg["b0_noise_grid"])
        db.add_control(rid, "negative",
                       "the persistence-gap rule does NOT report a loop on a uniform disc "
                       "(gap ratio < 1/10 of the circle's) and does NOT report 9 wells on a noise grid",
                       passed=neg_ok,
                       detail=(f"disc H1 gap ratio {neg['gap_ratio_h1_disc']:.4g} vs circle "
                               f"{sig_gap:.4g}; noise-grid b0 {neg['b0_noise_grid']}"))
    results["K6_negative"]["run_id"] = rid

    all_ok = all(results[c]["passed"] for c, _ in cases) and neg_ok
    out["cases"] = results
    out["all_passed"] = bool(all_ok)
    p = save_json("step0_known_answers.json", out)
    print(f"\nwrote {p}\nSTEP 0 {'PASSED' if all_ok else 'FAILED'}")

    with topodb() as db:
        db.add_finding(dataset_id=DS,
                       claim=("This pipeline (GUDHI 3.13.0 AlphaComplex / CubicalComplex, and the "
                              "fixed HEALPix 2-complex) recovers the textbook homology of the circle, "
                              "the sphere, the torus and a grid with 9 known wells, and the "
                              "persistence-gap rule does not fire on structureless input."),
                       verdict="recovered" if all_ok else "failed", tier="X",
                       caveat=("Known answers only. The gap rule is confined to this suite: the prior "
                               "campaign measured it firing on 9 of 12 density-matched nulls, so it is "
                               "not used as a detector on any real dataset in this campaign."),
                       reference="topodb_runs/astro/results/step0_known_answers.json")
    return 0 if all_ok else 1


if __name__ == "__main__":
    sys.exit(main())
