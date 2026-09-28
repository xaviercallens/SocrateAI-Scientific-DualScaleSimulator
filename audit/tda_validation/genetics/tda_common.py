"""Shared helpers for the genetics validation of the DualScaleSimulator TDA pipeline.

The pipeline functions under test are IMPORTED (not reimplemented) from
audit/reverse_zero/E5-cosmic-web-tda-scaled/cosmic_web_tda_scaled.py.
That module only runs its analysis under `if __name__ == "__main__"`, so a
plain import does not execute main(); we still load it by path with importlib
under a private module name so nothing else in sys.path is shadowed.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")  # machine is oversubscribed; 1 BLAS thread is ~10x faster here

import hashlib
import importlib.util
import json
import os
import time

import numpy as np
import gudhi

HERE = os.path.dirname(os.path.abspath(__file__))
WT_ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
PIPE_PATH = os.path.join(WT_ROOT, "audit", "reverse_zero", "E5-cosmic-web-tda-scaled", "cosmic_web_tda_scaled.py")
DATA_ROOT = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/genetics"


def _load_pipeline():
    spec = importlib.util.spec_from_file_location("cwtda_under_test", PIPE_PATH)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)  # module body only defines constants/functions
    return mod


PIPE = _load_pipeline()
alpha_persistence = PIPE.alpha_persistence
top_bars = PIPE.top_bars
betti_curve = PIPE.betti_curve
euler_curve = PIPE.euler_curve


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def pipeline_sha256():
    return sha256(PIPE_PATH)


def alpha_h1_summary(xyz, label="x", k=5):
    """Run the pipeline's alpha_persistence on a 3-D cloud (full Delaunay,
    max_alpha_sq=inf) and summarise H1 with top_bars, passing r_trunc
    explicitly (the module default 50.0 is a Mpc/h cosmology value)."""
    xyz = np.ascontiguousarray(xyz, dtype=float)
    res, by_dim = alpha_persistence(xyz, float("inf"), label, None)
    h1 = np.array(by_dim[1]) if by_dim[1] else np.empty((0, 2))
    n_inf = int(np.sum(~np.isfinite(h1[:, 1]))) if h1.size else 0
    finite_deaths = [d for dim in (0, 1, 2) for (_, d) in by_dim[dim] if np.isfinite(d)]
    r_trunc = float(max(finite_deaths)) if finite_deaths else 1.0
    tb = top_bars(by_dim, 1, k=k, r_trunc=r_trunc)
    pers = (tb[:, 1] - tb[:, 0]) if tb.size else np.zeros(0)
    p1 = float(pers[0]) if pers.size > 0 else 0.0
    p2 = float(pers[1]) if pers.size > 1 else 0.0
    c = xyz - xyz.mean(0)
    r_rms = float(np.sqrt(np.mean(np.sum(c ** 2, 1))))
    return {
        "n_points": int(xyz.shape[0]),
        "n_h1_bars": int(len(by_dim[1])),
        "n_h1_infinite_death": n_inf,
        "r_trunc_used": r_trunc,
        "top_h1_bars_birth_death": tb.tolist(),
        "top_h1_persistence": pers.tolist(),
        "P1": p1, "P2": p2,
        "P1_over_P2": (p1 / p2) if p2 > 0 else float("inf"),
        "R_rms": r_rms,
        "S": p1 / r_rms if r_rms > 0 else 0.0,
        "dominant_birth_over_death": float(tb[0, 0] / tb[0, 1]) if tb.size and tb[0, 1] > 0 else None,
        "betti_numbers_full_complex": res["betti_numbers_at_truncation"],
        "runtime_sec": res["runtime_sec"],
    }


def rips_h1_summary(points=None, distance_matrix=None, k=5):
    """Plain gudhi RipsComplex cross-check (full dimension / full metric),
    1-skeleton edge-collapsed then expanded to dimension 2."""
    t0 = time.time()
    if distance_matrix is not None:
        D = np.asarray(distance_matrix, dtype=float)
        rc = gudhi.RipsComplex(distance_matrix=D)
    else:
        P = np.asarray(points, dtype=float)
        rc = gudhi.RipsComplex(points=P)
        from scipy.spatial.distance import pdist, squareform
        D = squareform(pdist(P))
    st = rc.create_simplex_tree(max_dimension=1)
    st.collapse_edges()
    st.expansion(2)
    st.compute_persistence(homology_coeff_field=2)
    h1 = st.persistence_intervals_in_dimension(1)
    h1 = np.array(h1) if len(h1) else np.empty((0, 2))
    Dn = D.copy()
    np.fill_diagonal(Dn, np.inf)
    nn_med = float(np.median(Dn.min(1)))
    if h1.size:
        fin = np.isfinite(h1[:, 1])
        dmax = float(np.max(D))
        deaths = np.where(fin, h1[:, 1], dmax)
        pers = deaths - h1[:, 0]
        o = np.argsort(pers)[::-1]
        tb = np.column_stack([h1[o, 0], deaths[o]])[:k]
        pers_sorted = pers[o]
    else:
        tb = np.empty((0, 2))
        pers_sorted = np.zeros(0)
        fin = np.zeros(0, bool)
    p1 = float(pers_sorted[0]) if pers_sorted.size > 0 else 0.0
    p2 = float(pers_sorted[1]) if pers_sorted.size > 1 else 0.0
    return {
        "n_points": int(D.shape[0]),
        "n_h1_bars": int(h1.shape[0]),
        "n_h1_infinite_death": int(np.sum(~fin)) if h1.size else 0,
        "top_h1_bars_birth_death": tb.tolist(),
        "P1": p1, "P2": p2,
        "P1_over_P2": (p1 / p2) if p2 > 0 else float("inf"),
        "median_nn_distance": nn_med,
        "S_rips": p1 / nn_med if nn_med > 0 else 0.0,
        "all_h1_persistence_sorted": pers_sorted.tolist(),
        "runtime_sec": time.time() - t0,
    }


def pval_ge(null, obs):
    null = np.asarray(null, dtype=float)
    return float((1 + np.sum(null >= obs)) / (1 + null.size))


def circ_corr(a, b):
    """Jammalamadaka-SenGupta circular-circular correlation."""
    a = np.asarray(a, float)
    b = np.asarray(b, float)
    ma = np.angle(np.mean(np.exp(1j * a)))
    mb = np.angle(np.mean(np.exp(1j * b)))
    sa = np.sin(a - ma)
    sb = np.sin(b - mb)
    return float(np.sum(sa * sb) / np.sqrt(np.sum(sa ** 2) * np.sum(sb ** 2)))


def pca3(Z):
    Zc = Z - Z.mean(0)
    U, s, Vt = np.linalg.svd(Zc, full_matrices=False)
    X = U[:, :3] * s[:3]
    return X, (s ** 2 / np.sum(s ** 2))


def classical_mds(D, k=3):
    D = np.asarray(D, float)
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D ** 2) @ J
    w, V = np.linalg.eigh(B)
    o = np.argsort(w)[::-1]
    w, V = w[o], V[:, o]
    pos = np.clip(w[:k], 0, None)
    return V[:, :k] * np.sqrt(pos), w


def mad_low(x, nmad=3.0):
    med = np.median(x)
    mad = 1.4826 * np.median(np.abs(x - med))
    return med - nmad * mad, med + nmad * mad


def dump(obj, path):
    def conv(o):
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        raise TypeError(type(o))
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, default=conv)


class Budget(Exception):
    pass


_T_START = time.time()


def resumable_seed_loop(cache_path, n, fn, budget_sec):
    """Evaluate fn(seed) for seed = 0..n-1, caching finished values in
    cache_path (.json list) so a run that hits the wall-clock budget can be
    re-invoked with identical arguments and continues exactly where it
    stopped (same seeds, same results). Raises Budget if not finished."""
    vals = []
    if os.path.exists(cache_path):
        with open(cache_path) as f:
            vals = json.load(f)
    last_save = time.time()
    while len(vals) < n:
        vals.append(fn(len(vals)))
        if time.time() - last_save > 20 or len(vals) == n:
            with open(cache_path, "w") as f:
                json.dump(vals, f)
            last_save = time.time()
        if time.time() - _T_START > budget_sec and len(vals) < n:
            with open(cache_path, "w") as f:
                json.dump(vals, f)
            raise Budget(f"{cache_path}: {len(vals)}/{n}")
    return vals
