"""Shared primitives for the TopoDB quantum-fluid breadth run.

Every length returned by this module is a LENGTH, never a squared filtration
value: GUDHI's AlphaComplex stores squared circumradii and `assert_alpha_convention`
checks that on an equilateral triangle at import time of the caller.

Used by: step0_re6zr.py, compute_stm_fields.py, compute_xy.py, compute_gpe.py,
compute_breadth.py, compute_disorder.py.
"""
from __future__ import annotations

import hashlib
import json
import os
import re
import resource
import sys
import time

import numpy as np
from scipy import ndimage, spatial

PHI0 = 2.067833848e-15  # Wb
DATA_ROOT = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/quantum_fluid"
STM_ROOT = os.path.join(DATA_ROOT, "raw/zenodo_14780459/extracted/STM data/Fig 1 and 2")
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.abspath(os.path.join(HERE, "..", "results"))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))

PRLIMIT = "prlimit --as=8589934592 -- "


# --------------------------------------------------------------------- meta
def command() -> str:
    """The exact command, including the memory-cap wrapper it was run under."""
    return PRLIMIT + sys.executable + " " + " ".join(sys.argv)


class Meter:
    def __enter__(self):
        self.t0 = time.perf_counter()
        return self

    def __exit__(self, *a):
        self.wall_sec = time.perf_counter() - self.t0
        self.peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def sha256_file(path: str) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(name: str, obj) -> str:
    os.makedirs(RESULTS, exist_ok=True)
    p = os.path.join(RESULTS, name)
    with open(p, "w") as fh:
        json.dump(obj, fh, indent=1, sort_keys=True, default=_jd)
    return p


def _jd(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(type(o))


# ------------------------------------------------------------------- gudhi
def alpha_bars(points, max_alpha_sq=None, max_dim=2):
    """Alpha persistence. Returns {dim: ndarray (n,2)} of (birth, death) RADII (lengths).

    death = inf is kept as np.inf. GUDHI filtration values are squared radii; we sqrt.
    """
    import gudhi

    pts = np.ascontiguousarray(np.asarray(points, dtype=float))
    ac = gudhi.AlphaComplex(points=pts)
    st = ac.create_simplex_tree() if max_alpha_sq is None else ac.create_simplex_tree(max_alpha_square=float(max_alpha_sq))
    st.compute_persistence(homology_coeff_field=2)
    out = {}
    for d in range(max_dim + 1):
        iv = st.persistence_intervals_in_dimension(d)
        if len(iv) == 0:
            out[d] = np.empty((0, 2))
            continue
        iv = np.asarray(iv, dtype=float)
        # numerical negatives from roundoff at 0
        iv[:, 0] = np.sqrt(np.clip(iv[:, 0], 0, None))
        fin = np.isfinite(iv[:, 1])
        iv[fin, 1] = np.sqrt(np.clip(iv[fin, 1], 0, None))
        out[d] = iv
    return out


def assert_alpha_convention(tol=1e-9):
    """Equilateral triangle of side s: the H1 class of the 3 vertices must die at
    the circumradius s/sqrt(3). Catches a missing sqrt, which would otherwise
    produce plausible-looking but wrong lengths everywhere."""
    s = 7.0
    tri = np.array([[0.0, 0.0], [s, 0.0], [s / 2, s * np.sqrt(3) / 2]])
    # a bare triangle has no H1 (it is filled); use a hexagon hole instead:
    ang = np.arange(6) * np.pi / 3
    R = 5.0
    hexa = np.c_[R * np.cos(ang), R * np.sin(ang)]
    bars = alpha_bars(hexa)
    h1 = bars[1]
    assert len(h1) == 1, f"hexagon should give exactly one H1 class, got {len(h1)}"
    got = h1[0, 1]
    assert abs(got - R) < 1e-6, f"hexagon H1 must die at circumradius {R}, got {got}"
    # H0 merge radius of a single edge of length s: half the length
    bars2 = alpha_bars(np.array([[0.0, 0.0], [s, 0.0], [0.0, 3 * s]]))
    d0 = np.sort(bars2[0][np.isfinite(bars2[0][:, 1]), 1])
    assert abs(d0[0] - s / 2) < 1e-6, f"H0 first merge must be s/2 = {s/2}, got {d0[0]}"
    del tri
    return True


def h0_deaths(bars) -> np.ndarray:
    a = bars[0]
    return a[np.isfinite(a[:, 1]), 1]


def h1_deaths(bars) -> np.ndarray:
    a = bars[1]
    return a[np.isfinite(a[:, 1]), 1] if len(a) else np.empty(0)


def betti_from_bars(bars, thresh=None) -> dict:
    """Betti numbers. thresh=None -> essential classes (infinite bars) only.
    thresh=float -> classes alive at that radius."""
    out = {}
    for d, iv in bars.items():
        if len(iv) == 0:
            out[d] = 0
            continue
        if thresh is None:
            out[d] = int(np.sum(~np.isfinite(iv[:, 1])))
        else:
            out[d] = int(np.sum((iv[:, 0] <= thresh) & ((~np.isfinite(iv[:, 1])) | (iv[:, 1] > thresh))))
    return out


def top_bars(bars, dim, n=10, r_trunc=None):
    iv = bars.get(dim, np.empty((0, 2)))
    if len(iv) == 0:
        return []
    rows = []
    for b, d in iv:
        if not np.isfinite(d):
            d = None if r_trunc is None else float(r_trunc)
        rows.append((float(b), None if d is None else float(d)))
    rows.sort(key=lambda bd: (np.inf if bd[1] is None else bd[1] - bd[0]), reverse=True)
    return rows[:n]


# -------------------------------------------------------------- statistics
def spread_stats(x) -> dict:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return {"n": 0}
    q1, med, q3 = np.percentile(x, [25, 50, 75])
    return {"n": int(x.size), "median": float(med), "q1": float(q1), "q3": float(q3),
            "iqr_over_median": float((q3 - q1) / med) if med else None,
            "mean": float(x.mean()), "cv": float(x.std() / x.mean()) if x.mean() else None}


def psi6_global(pts) -> tuple:
    """|<exp(6 i theta_b)>| over Delaunay bonds. Global (not per-particle) so that
    it is a single orientational-order number comparable across datasets."""
    pts = np.asarray(pts, dtype=float)
    if len(pts) < 4:
        return float("nan"), 0
    try:
        tri = spatial.Delaunay(pts)
    except Exception:
        return float("nan"), 0
    edges = set()
    for s in tri.simplices:
        for i in range(3):
            a, b = sorted((int(s[i]), int(s[(i + 1) % 3])))
            edges.add((a, b))
    e = np.array(sorted(edges))
    if len(e) == 0:
        return float("nan"), 0
    d = pts[e[:, 1]] - pts[e[:, 0]]
    th = np.arctan2(d[:, 1], d[:, 0])
    return float(abs(np.mean(np.exp(6j * th)))), int(len(e))


def nn_distance(pts, box=None, margin=0.0) -> np.ndarray:
    """Nearest-neighbour distances. If margin > 0 and box = (xmin,xmax,ymin,ymax),
    only points at least `margin` from the border are used as sources (the
    neighbour may lie outside), which removes the border bias."""
    pts = np.asarray(pts, dtype=float)
    if len(pts) < 2:
        return np.empty(0)
    tree = spatial.cKDTree(pts)
    d, _ = tree.query(pts, k=2)
    d = d[:, 1]
    if margin > 0 and box is not None:
        xmin, xmax, ymin, ymax = box
        keep = ((pts[:, 0] >= xmin + margin) & (pts[:, 0] <= xmax - margin) &
                (pts[:, 1] >= ymin + margin) & (pts[:, 1] <= ymax - margin))
        d = d[keep]
    return d


def a_density(n, area) -> float:
    """Triangular-lattice constant at areal density n/area."""
    return float(np.sqrt(2.0 / (np.sqrt(3.0) * (n / area))))


def rank_p(observed, null_samples, side="two") -> float:
    """Rank p-value with the 1/(n+1) floor made explicit by the caller's n_null."""
    null = np.asarray(null_samples, dtype=float)
    null = null[np.isfinite(null)]
    n = null.size
    if n == 0:
        return float("nan")
    ge = float(np.sum(null >= observed))
    le = float(np.sum(null <= observed))
    if side == "greater":
        return (ge + 1) / (n + 1)
    if side == "less":
        return (le + 1) / (n + 1)
    return min(1.0, 2.0 * min((ge + 1) / (n + 1), (le + 1) / (n + 1)))


# ----------------------------------------------------------------- S(q)/FFT
def a_fft(img, px_nm, a_expect_nm, pad=4, qlo_frac=0.4, qhi_frac=2.0):
    """Detector-free, TDA-free lattice constant from the first Bragg ring.

    Triangular lattice: q1 = 4 pi / (sqrt(3) a)  ->  a = 4 pi / (sqrt(3) q1).
    Plane-subtract, Hann window, zero-pad by `pad`, azimuthally average |FFT|^2,
    take the maximum in [qlo_frac, qhi_frac] * q1_expect and refine it with a
    3-point parabolic fit on the radial profile.
    Returns (a_nm, q1, profile_q, profile_S).
    """
    im = plane_subtract(img)
    ny, nx = im.shape
    w = np.outer(np.hanning(ny), np.hanning(nx))
    im = (im - im.mean()) * w
    Ny, Nx = ny * pad, nx * pad
    F = np.fft.fftshift(np.abs(np.fft.fft2(im, s=(Ny, Nx))) ** 2)
    qx = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(Nx, d=px_nm))
    qy = 2 * np.pi * np.fft.fftshift(np.fft.fftfreq(Ny, d=px_nm))
    QX, QY = np.meshgrid(qx, qy)
    Q = np.sqrt(QX ** 2 + QY ** 2)
    dq = 2 * np.pi / (max(Nx, Ny) * px_nm)
    nb = int(Q.max() / dq)
    idx = np.clip((Q / dq).astype(int), 0, nb - 1)
    S = np.bincount(idx.ravel(), weights=F.ravel(), minlength=nb)
    cnt = np.bincount(idx.ravel(), minlength=nb)
    good = cnt > 0
    S = np.where(good, S / np.maximum(cnt, 1), 0.0)
    qc = (np.arange(nb) + 0.5) * dq
    q1e = 4 * np.pi / (np.sqrt(3) * a_expect_nm)
    sel = (qc >= qlo_frac * q1e) & (qc <= qhi_frac * q1e) & good
    if not sel.any():
        return float("nan"), float("nan"), qc, S
    j = int(np.argmax(np.where(sel, S, -np.inf)))
    q1 = qc[j]
    if 0 < j < nb - 1 and good[j - 1] and good[j + 1]:
        y0, y1, y2 = S[j - 1], S[j], S[j + 1]
        den = (y0 - 2 * y1 + y2)
        if den != 0:
            q1 = qc[j] + 0.5 * dq * (y0 - y2) / den
    return float(4 * np.pi / (np.sqrt(3) * q1)), float(q1), qc, S


def plane_subtract(img):
    img = np.asarray(img, dtype=float)
    ny, nx = img.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    ok = np.isfinite(img)
    A = np.c_[xx[ok], yy[ok], np.ones(int(ok.sum()))]
    coef, *_ = np.linalg.lstsq(A, img[ok], rcond=None)
    out = img - (coef[0] * xx + coef[1] * yy + coef[2])
    if not ok.all():
        out[~ok] = float(np.nanmean(out[ok]))
    return out


# -------------------------------------------------------------------- .sxm
def read_sxm(path):
    """Nanonis .sxm: ASCII header to ':SCANIT_END:', 0x1A 0x04 marker, then
    big-endian float32 images, channels in DATA_INFO ('both' = fwd then bwd)."""
    b = open(path, "rb").read()
    i = b.find(b":SCANIT_END:")
    hdr = b[:i].decode("latin1")
    j = b.find(b"\x1a\x04", i)
    data = b[j + 2:]
    tags = {}
    for m in re.finditer(r":([^:\n]+):\n(.*?)(?=\n:[^:\n]+:\n|\Z)", hdr, re.S):
        tags[m.group(1)] = m.group(2)
    nx, ny = [int(v) for v in tags["SCAN_PIXELS"].split()]
    rx, ry = [float(v) for v in tags["SCAN_RANGE"].split()]
    lines = [ln.strip().split("\t") for ln in tags["DATA_INFO"].strip().split("\n")]
    head, rows = lines[0], lines[1:]
    chans = []
    for r in rows:
        d = dict(zip(head, r))
        for dd in (["fwd", "bwd"] if d["Direction"] == "both" else ["fwd"]):
            chans.append((d["Name"], dd))
    arr = np.frombuffer(data[: len(chans) * nx * ny * 4], dtype=">f4").reshape(len(chans), ny, nx)
    return ({c: arr[k].astype(np.float64) for k, c in enumerate(chans)},
            dict(nx=nx, ny=ny, range_x_m=rx, range_y_m=ry))


def detect_minima(img, px_nm, a_nm):
    """Vortex cores = local conductance minima. Plane-subtract, Gaussian low-pass
    with real-space sigma = a/6, keep pixels equal to the minimum of a disk of
    radius a/3. Same detector as the earlier validation run (deliberately: the
    Step-0 gate does not rely on it)."""
    im = plane_subtract(img)
    ny, nx = im.shape
    s_px = (a_nm / 6.0) / px_nm
    kx = 2 * np.pi * np.fft.fftfreq(nx)
    ky = 2 * np.pi * np.fft.fftfreq(ny)
    KX, KY = np.meshgrid(kx, ky)
    f = np.real(np.fft.ifft2(np.fft.fft2(im) * np.exp(-0.5 * s_px ** 2 * (KX ** 2 + KY ** 2))))
    rad = (a_nm / 3.0) / px_nm
    R = int(np.ceil(rad))
    fy, fx = np.mgrid[-R:R + 1, -R:R + 1]
    fp = (fx ** 2 + fy ** 2) <= rad ** 2
    mn = ndimage.minimum_filter(f, footprint=fp, mode="nearest")
    my, mx = np.nonzero(f == mn)
    return np.stack([mx * px_nm, my * px_nm], 1)


def a_tri(B_T) -> float:
    """Triangular vortex-lattice constant in nm at induction B (T)."""
    return float(1.075 * np.sqrt(PHI0 / B_T) * 1e9)


# ----------------------------------------------------------- phase winding
def winding_vortices(theta, periodic=True):
    """Plaquette phase circulation on a 2-D phase field. Returns (pos, charge)
    with pos in (col, row) plaquette-centre units. theta in radians."""
    th = np.asarray(theta, dtype=float)
    if periodic:
        t00 = th
        t10 = np.roll(th, -1, axis=1)
        t11 = np.roll(np.roll(th, -1, axis=1), -1, axis=0)
        t01 = np.roll(th, -1, axis=0)
    else:
        t00 = th[:-1, :-1]
        t10 = th[:-1, 1:]
        t11 = th[1:, 1:]
        t01 = th[1:, :-1]

    def wrap(d):
        return (d + np.pi) % (2 * np.pi) - np.pi

    circ = wrap(t10 - t00) + wrap(t11 - t10) + wrap(t01 - t11) + wrap(t00 - t01)
    q = np.rint(circ / (2 * np.pi)).astype(int)
    ry, rx = np.nonzero(q != 0)
    return np.stack([rx + 0.5, ry + 0.5], 1).astype(float), q[ry, rx]


# --------------------------------------------------------------- cubical PH
def cubical_bars(field, sublevel=True):
    """Cubical persistence of a 2-D field. Returns {dim: (n,2)} in FIELD units.
    sublevel=False computes the superlevel filtration by negating the field
    (values are returned negated back, so a bar's endpoints are field values)."""
    import gudhi

    f = np.asarray(field, dtype=float)
    g = f if sublevel else -f
    cc = gudhi.CubicalComplex(top_dimensional_cells=g)
    cc.compute_persistence(homology_coeff_field=2)
    out = {}
    for d in (0, 1):
        iv = cc.persistence_intervals_in_dimension(d)
        iv = np.asarray(iv, dtype=float).reshape(-1, 2)
        if not sublevel and len(iv):
            iv = -iv[:, ::-1]
        out[d] = iv
    return out
