"""
lib_ph.py -- point samplers, embeddings and the pre-registered PH recovery criterion
for Part B (tier X). Criterion as in expectations.json partB.criterion_recovered:
  beta_k(eps) = #{H_k bars with birth <= eps < death}; RECOVERED iff some window
  [e1, e2] with e2 >= R*e1 (R = 1.5) has beta_k(eps) == expected_k for all tested k.
"""
from __future__ import annotations

import numpy as np

R_WINDOW = 1.5


# ---------------------------------------------------------------- samplers
def sample_torus(n, N, rng):
    th = rng.uniform(0, 2 * np.pi, size=(N, n))
    X = np.empty((N, 2 * n))
    X[:, 0::2] = np.cos(th)
    X[:, 1::2] = np.sin(th)
    return X, th


def orbifold_embed(th):
    """Z2-invariant map T^4 -> R^14: (cos th_i)_i, (sin th_i sin th_j)_{i<=j}.
    Invariant under th -> -th; separates th from th' unless th' = +-th (see injectivity check)."""
    c, s = np.cos(th), np.sin(th)
    cols = [c[:, i] for i in range(4)]
    for i in range(4):
        for j in range(i, 4):
            cols.append(s[:, i] * s[:, j])
    return np.stack(cols, axis=1)


def quotient_dist(a, b):
    """distance on T^4/+-1 between angle vectors (flat metric)."""
    def td(x, y):
        d = np.abs((x - y + np.pi) % (2 * np.pi) - np.pi)
        return np.sqrt((d ** 2).sum(-1))
    return np.minimum(td(a, b), td(a, -b))


def sample_fermat_k3(M, rng):
    """M draws of random complex Gaussian (x,y,z); all 4 roots w of w^4 = -(x^4+y^4+z^4).
    Returns 4M points of CP^3 as unit vectors in C^4 and max residual |sum z_i^4|/|z|^4."""
    xyz = rng.normal(size=(M, 3)) + 1j * rng.normal(size=(M, 3))
    s = -(xyz ** 4).sum(1)
    r = np.abs(s) ** 0.25
    ang = np.angle(s) / 4
    pts = []
    for m in range(4):
        w = r * np.exp(1j * (ang + m * np.pi / 2))
        pts.append(np.concatenate([xyz, w[:, None]], axis=1))
    Z = np.concatenate(pts, axis=0)
    Z = Z / np.linalg.norm(Z, axis=1, keepdims=True)
    resid = np.abs((Z ** 4).sum(1)).max()
    return Z, float(resid)


def projector_embed(Z):
    """|z><z| (|z|=1) as a point of R^16 with Frobenius-isometric coordinates:
    diagonal entries, sqrt(2)*Re and sqrt(2)*Im of the 6 upper off-diagonal entries."""
    P = Z[:, :, None] * Z.conj()[:, None, :]
    cols = [P[:, i, i].real for i in range(4)]
    for i in range(4):
        for j in range(i + 1, 4):
            cols.append(np.sqrt(2) * P[:, i, j].real)
            cols.append(np.sqrt(2) * P[:, i, j].imag)
    return np.stack(cols, axis=1)


def sample_null_cube(N, D, diam, rng):
    """Uniform points in a 4-cube of the given diameter, placed by a random isometry in R^D."""
    side = diam / 2.0
    Y = rng.uniform(0, side, size=(N, 4))
    Q, _ = np.linalg.qr(rng.normal(size=(D, 4)))
    return Y @ Q.T


def farthest_point_subsample(X, n, rng):
    """max-min (farthest point) subsample of n points, first point random."""
    N = len(X)
    idx = [int(rng.integers(N))]
    d = np.linalg.norm(X - X[idx[0]], axis=1)
    for _ in range(n - 1):
        i = int(np.argmax(d))
        idx.append(i)
        d = np.minimum(d, np.linalg.norm(X - X[i], axis=1))
    return X[idx], float(d.max())


# ---------------------------------------------------------------- criterion
def betti_curve(bars, t):
    b = np.array([x[0] for x in bars]) if bars else np.zeros(0)
    d = np.array([x[1] for x in bars]) if bars else np.zeros(0)
    b.sort(); d.sort()
    return np.searchsorted(b, t, side="right") - np.searchsorted(d, t, side="right")


def recovery(bars_by_dim, expected, tau):
    """bars_by_dim[k] = list of (birth, death) (death may be inf); expected = list per dim.
    Returns dict with best window ratio for the full expected vector and per-dim gap ratios."""
    ev = {0.0, float(tau)}
    for k in range(len(expected)):
        for b, d in bars_by_dim.get(k, []):
            if b < tau:
                ev.add(float(b))
            if d < tau:
                ev.add(float(d))
    t = np.array(sorted(ev))
    ok = np.ones(len(t), dtype=bool)
    curves = {}
    for k, e in enumerate(expected):
        c = betti_curve(bars_by_dim.get(k, []), t)
        curves[k] = c
        ok &= (c == e)
    # intervals [t_i, t_{i+1}) with ok[i]; merge runs
    best = (0.0, None, None)
    i = 0
    while i < len(t) - 1:
        if ok[i]:
            j = i
            while j < len(t) - 1 and ok[j]:
                j += 1
            a, b = t[i], t[j]
            ratio = (b / a) if a > 0 else float("inf")
            if ratio > best[0]:
                best = (float(ratio), float(a), float(b))
            i = j
        else:
            i += 1
    gaps = {}
    for k, e in enumerate(expected):
        pers = sorted(((min(d, tau) - b) for b, d in bars_by_dim.get(k, []) if b < tau), reverse=True)
        if k == 0:
            pers = pers[1:]  # drop the essential H0 bar
            e = e - 1
        if e == 0:
            gaps[k] = {"top_noise_persistence": float(pers[0]) if pers else 0.0}
        else:
            gaps[k] = {"pers_expected_th": float(pers[e - 1]) if len(pers) >= e else None,
                       "pers_next": float(pers[e]) if len(pers) > e else 0.0,
                       "ratio": (float(pers[e - 1] / pers[e]) if len(pers) > e and pers[e] > 0 else
                                 (float("inf") if len(pers) >= e else 0.0))}
    # max beta per dim over eps in (0,tau) for diagnostics, and beta at a few scales
    return {"best_window_ratio": best[0], "best_window": [best[1], best[2]],
            "recovered": best[0] >= R_WINDOW, "gap": gaps,
            "max_beta_k": {k: int(curves[k][1:].max()) if len(t) > 1 else 0 for k in curves}}


def bars_from_st(st, maxk):
    out = {}
    for k in range(maxk + 1):
        out[k] = [(float(b), float(d)) for b, d in st.persistence_intervals_in_dimension(k)]
    return out
