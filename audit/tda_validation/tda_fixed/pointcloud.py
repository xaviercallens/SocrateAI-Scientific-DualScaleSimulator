"""Point-cloud (alpha complex) persistence.

NO DEFECT WAS FOUND IN THIS PATH.  The known-answer suite on branch
loop/tda-simple exercised ``cosmic_web_tda_scaled.alpha_persistence`` in cases
P1, P2, P3a, P3b, P4, P10a, P10b, N1, N2 and it recovered the expected
topology; the two defects it found (see ``fix_expectations.json``) are both in
``cmb_tda.py``.  This module exists only so ``run_validation.py`` has ONE
import surface for all three paths.  It is a re-implementation of the same
alpha-complex call, not a fix, and ``run_validation.py`` checks that it gives
bars identical to the original pipeline function on the P2 cloud and says
which function produced the reported numbers.

Tier of every number this module produces: X (numerics).
"""
from __future__ import annotations

import time

import numpy as np

__all__ = ["alpha_persistence", "betti_curve", "sphere_cloud_P2"]


def alpha_persistence(xyz, max_alpha_sq=float("inf"), label="cloud"):
    """Alpha-complex persistence over Z/2, bars in RADIUS units.

    Mirrors ``audit/reverse_zero/E5-cosmic-web-tda-scaled/
    cosmic_web_tda_scaled.py::alpha_persistence`` (same gudhi calls, same
    sqrt of the squared alpha values, same coefficient field), minus the
    on-disk side effects.
    """
    import gudhi

    t0 = time.time()
    ac = gudhi.AlphaComplex(points=xyz)
    st = ac.create_simplex_tree(max_alpha_square=max_alpha_sq)
    diag = st.persistence(homology_coeff_field=2, min_persistence=0.0)
    dt = time.time() - t0
    by_dim = {0: [], 1: [], 2: []}
    for dim, (b, d) in diag:
        if dim in by_dim:
            by_dim[dim].append([float(np.sqrt(max(b, 0.0))),
                                float(np.sqrt(d)) if np.isfinite(d) else float("inf")])
    return {
        "label": label,
        "n_points": int(np.asarray(xyz).shape[0]),
        "runtime_sec": dt,
        "betti_numbers_at_truncation": list(st.betti_numbers()),
        "function": "tda_fixed.pointcloud.alpha_persistence (gudhi.AlphaComplex, Z/2)",
        "tier": "X",
    }, by_dim


def betti_curve(by_dim, dim, r_grid, r_trunc):
    """beta_dim(r) = #{bars (b, d): b <= r < d}, infinite deaths alive throughout."""
    bars = np.array(by_dim[dim], dtype=float) if by_dim.get(dim) else np.empty((0, 2))
    curve = np.zeros(len(r_grid))
    if bars.size == 0:
        return curve
    b = bars[:, 0]
    d = np.where(np.isfinite(bars[:, 1]), bars[:, 1], r_trunc + 1e9)
    for i, r in enumerate(r_grid):
        curve[i] = np.sum((b <= r) & (r < d))
    return curve


def sphere_cloud_P2():
    """The suite's P2 cloud, reproduced exactly: 5000 normalised Gaussian
    vectors (np.random.default_rng(12)) + isotropic noise sigma 0.02.
    Seed and construction from audit/tda_validation/simple_suite/
    expectations.json case P2."""
    rng = np.random.default_rng(12)
    v = rng.normal(size=(5000, 3))
    v = v / np.linalg.norm(v, axis=1, keepdims=True)
    return v + rng.normal(0, 0.02, (5000, 3))
