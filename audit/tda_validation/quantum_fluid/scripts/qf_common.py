"""Shared helpers for the quantum-fluid TDA validation.

The pipeline functions under test are IMPORTED (never reimplemented) from
  audit/reverse_zero/E5-cmb-tda/cmb_tda.py                     -> betti_curves_from_topology
  audit/reverse_zero/E5-cosmic-web-tda-scaled/cosmic_web_tda_scaled.py -> alpha_persistence, betti_curve, top_bars
via importlib.util.spec_from_file_location. Both modules guard their main()
behind `if __name__ == "__main__"`, so importing them runs only module-level
constants/imports (checked by reading both files).

This file only builds INPUTS in the formats those functions expect
(vertex list `unmasked`, local-index edge pairs, local-index triangle
triples) and small helper statistics that are not part of the pipeline.
"""
import importlib.util
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
WT_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", ".."))
CMB_PATH = os.path.join(WT_ROOT, "audit/reverse_zero/E5-cmb-tda/cmb_tda.py")
WEB_PATH = os.path.join(WT_ROOT, "audit/reverse_zero/E5-cosmic-web-tda-scaled/cosmic_web_tda_scaled.py")
DATA_ROOT = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/quantum_fluid"
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_cmb = None
_web = None


def cmb():
    global _cmb
    if _cmb is None:
        _cmb = _load("cmb_tda_under_test", CMB_PATH)
    return _cmb


def web():
    global _web
    if _web is None:
        _web = _load("cosmic_web_tda_under_test", WEB_PATH)
    return _web


def torus_topology(L):
    """Periodic L x L grid (flat index i = y*L + x), Freudenthal triangulation:
    edges right, up, and one diagonal (x,y)-(x+1,y+1); two triangles per square.
    Returns (unmasked, edges_arr, tris_arr) in the format of
    cmb_tda.build_topology (local indices == flat indices here, all vertices
    kept). Euler characteristic of the torus: V - E + F = 0 (asserted)."""
    idx = np.arange(L * L).reshape(L, L)  # idx[y, x]
    r = np.roll(idx, -1, axis=1)          # (x+1, y)
    u = np.roll(idx, -1, axis=0)          # (x, y+1)
    d = np.roll(r, -1, axis=0)            # (x+1, y+1)
    a = idx.ravel()
    e = np.concatenate([np.stack([a, r.ravel()], 1), np.stack([a, u.ravel()], 1), np.stack([a, d.ravel()], 1)])
    e = np.sort(e, axis=1)
    t1 = np.sort(np.stack([a, r.ravel(), d.ravel()], 1), axis=1)
    t2 = np.sort(np.stack([a, u.ravel(), d.ravel()], 1), axis=1)
    tris = np.concatenate([t1, t2])
    V, E, F = L * L, e.shape[0], tris.shape[0]
    assert len({tuple(x) for x in e}) == E, "duplicate edges"
    assert len({tuple(x) for x in tris}) == F, "duplicate triangles"
    assert V - E + F == 0, (V, E, F)
    return a.astype(np.int64), e.astype(np.int64), tris.astype(np.int64)


def masked_grid_topology(mask2d):
    """Open (non-periodic) grid restricted to mask2d (bool [ny, nx]); flat index
    i = y*nx + x. Edges right/up/diagonal and the two Freudenthal triangles,
    kept only if all vertices are in the mask. Returns (unmasked, edges_local,
    tris_local) with LOCAL indices 0..N-1 into `unmasked`, exactly like
    cmb_tda.build_topology."""
    ny, nx = mask2d.shape
    flat = np.arange(ny * nx).reshape(ny, nx)
    unmasked = flat[mask2d]
    loc = -np.ones(ny * nx, dtype=np.int64)
    loc[unmasked] = np.arange(unmasked.size)
    m = mask2d
    A = flat[:-1, :-1]; R = flat[:-1, 1:]; U = flat[1:, :-1]; D = flat[1:, 1:]
    mA = m[:-1, :-1]; mR = m[:-1, 1:]; mU = m[1:, :-1]; mD = m[1:, 1:]
    edges = []
    # horizontal edges over full grid
    edges.append(np.stack([flat[:, :-1][m[:, :-1] & m[:, 1:]], flat[:, 1:][m[:, :-1] & m[:, 1:]]], 1))
    edges.append(np.stack([flat[:-1, :][m[:-1, :] & m[1:, :]], flat[1:, :][m[:-1, :] & m[1:, :]]], 1))
    k = mA & mD
    edges.append(np.stack([A[k], D[k]], 1))
    e = loc[np.concatenate(edges)]
    k1 = mA & mR & mD
    k2 = mA & mU & mD
    t = loc[np.concatenate([np.stack([A[k1], R[k1], D[k1]], 1), np.stack([A[k2], U[k2], D[k2]], 1)])]
    e = np.sort(e, 1); t = np.sort(t, 1)
    return unmasked.astype(np.int64), e.astype(np.int64), t.astype(np.int64)


def wrap(a):
    """wrap angle differences to (-pi, pi]"""
    return np.angle(np.exp(1j * a))


def plaquette_winding(theta):
    """theta [L, L] periodic (index [y, x]). Winding number of each plaquette
    with corners (x,y),(x+1,y),(x+1,y+1),(x,y+1), counter-clockwise, using
    bond differences wrapped to (-pi, pi]. Returns int array [L, L] with
    plaquette (y, x) centred at (x+0.5, y+0.5)."""
    t00 = theta
    t10 = np.roll(theta, -1, axis=1)
    t11 = np.roll(t10, -1, axis=0)
    t01 = np.roll(theta, -1, axis=0)
    s = wrap(t10 - t00) + wrap(t11 - t10) + wrap(t01 - t11) + wrap(t00 - t01)
    return np.rint(s / (2 * np.pi)).astype(int)


def psi6_local(points, box=None):
    """Local hexatic order |psi6_j| = |mean_k exp(6 i theta_jk)| over Delaunay
    neighbours (scipy.spatial.Delaunay; independent of GUDHI). Returns array
    per point. Non-periodic."""
    from scipy.spatial import Delaunay
    tri = Delaunay(points)
    indptr, nb = tri.vertex_neighbor_vertices
    out = np.zeros(len(points), dtype=complex)
    for j in range(len(points)):
        n = nb[indptr[j]:indptr[j + 1]]
        dv = points[n] - points[j]
        ang = np.arctan2(dv[:, 1], dv[:, 0])
        out[j] = np.mean(np.exp(6j * ang))
    return out, tri


def psi6_global_bonds(points):
    """Global bond-orientational order |<exp(6 i theta_b)>| over all Delaunay
    bonds (one reading of the Duhan et al. 2025 definition)."""
    from scipy.spatial import Delaunay
    tri = Delaunay(points)
    indptr, nb = tri.vertex_neighbor_vertices
    angs = []
    for j in range(len(points)):
        for k in nb[indptr[j]:indptr[j + 1]]:
            if k > j:
                dv = points[k] - points[j]
                angs.append(np.arctan2(dv[1], dv[0]))
    angs = np.array(angs)
    return float(np.abs(np.mean(np.exp(6j * angs)))), angs.size


def h0_finite_deaths(by_dim):
    bars = np.array(by_dim[0]) if by_dim[0] else np.empty((0, 2))
    if bars.size == 0:
        return np.empty(0)
    d = bars[:, 1]
    return d[np.isfinite(d)]


def spread_stats(x):
    x = np.asarray(x, float)
    if x.size == 0:
        return {"n": 0}
    q1, med, q3 = np.percentile(x, [25, 50, 75])
    return {"n": int(x.size), "median": float(med), "q1": float(q1), "q3": float(q3),
            "iqr_over_median": float((q3 - q1) / med) if med > 0 else None,
            "mean": float(x.mean()), "cv": float(x.std() / x.mean()) if x.mean() > 0 else None}
