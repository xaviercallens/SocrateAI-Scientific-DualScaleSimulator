#!/usr/bin/env python3
"""
Discrete-symmetry ("crystallinity") lens for the CMB and the cosmic web.

WHAT THIS IS.  A search template, tier X (numerics), for residual DISCRETE
angular anisotropy: does a field on S^2 carry more power in the subspace
invariant under a finite rotation group G than an isotropic Gaussian field with
the same power spectrum and the same mask does?  The finite groups tried first
are taken from LeanMaster Stream 8 and from this project's D_4 work, as a
motivated shortlist of WHICH groups to try, not as a derivation of a signal.

WHAT THIS IS NOT.  There is no derivation from K3 x T^2 to any CMB observable.
LeanMaster says so itself, verbatim (v3.28.0, commit 64f905f,
docs/STREAM8_WHICH_K3.md:115):

    "(iii) Observables: none -- `N = 4`, non-chiral."

Nothing in this file is a prediction of K3 x T^2 and nothing here is proved.
The full design note and pre-registration is lens_spec.json in this directory,
committed before this file was run.

THE GROUP THAT ACTS ON THE SKY.  (Z_2)^4 |x A_4, order 192, does not act on S^2:
the (Z_2)^4 factor is 2-torsion translation on the abelian surface.  Only A_4
can.  It is reached from the 24 Hurwitz units (the binary tetrahedral group 2T)
through the standard double cover

    phi : Sp(1) -> SO(3),   phi(q) : v |-> q v qbar   (v purely imaginary,
                                                       basis (i, j, k)),

whose kernel is {+-1}.  So |phi(2T)| = 24/2 = 12, NOT 24.  The control group of
the same order is therefore C_12, not C_24.

STAGES (each is a separate process; see report.json for the exact commands):

    --stage groups      tier-B group / representation checks, no data
    --stage operators   build and cache the Wigner-D projector operators
    --stage nulls       Gaussian null ensemble for a given map's spectrum
    --stage injection   known-answer sensitivity + group specificity
    --stage data        the ONE pre-registered application to a real map
    --stage tda         lower-star Betti curves of m_G and r_G (secondary)
    --stage pointcloud  point-set path, synthetic known-answer test
"""
import argparse
import hashlib
import importlib.util
import json
import os
import platform
import resource
import subprocess
import sys
import time
from fractions import Fraction
from itertools import product

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REVERSE_WT = "/mnt/disks/disk-socrateai-local-1/dualscale-wt-reverse"
CMB_TDA_REL = "audit/reverse_zero/E5-cmb-tda/cmb_tda.py"

# ---------------------------------------------------------------------------
# frozen analysis constants (pre-registered in lens_spec.json; do not edit)
# ---------------------------------------------------------------------------
NSIDE_WORK = 128
LMAX = 64
BANDS = [(2, 8), (9, 16), (17, 32), (33, 64)]
N_ORIENT = 192
ORIENT_SEED = 20260919
N_NULL = 200
NULL_BASE_SEED = 900000
INJ_BASE_SEED = 700000
INJ_PATTERN_SEED = 31415
AMPLITUDES = [0.0, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20]
N_INJ = 100
ALPHA = 0.05
N_TESTS = 4
ALPHA_BONF = ALPHA / N_TESTS

WMAP_MAP = os.path.join(REVERSE_WT, "data/real2/cmb/wmap_ilc_9yr_v5.fits")
WMAP_MASK = os.path.join(
    REVERSE_WT, "data/real2/cmb/wmap_temperature_kq85_analysis_mask_r9_9yr_v5.fits")
PLANCK_MAP = ("/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/planck_maps/"
              "COM_CMB_IQU-smica_2048_R3.00_full.fits")
PLANCK_MASK = ("/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/planck_maps/"
               "COM_Mask_CMB-common-Mask-Int_2048_R3.00.fits")


def log(msg):
    print("[%7.1fs] %s" % (time.time() - _T0, msg), flush=True)


_T0 = time.time()


def peak_rss_mb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def sha256_of(path, chunk=1 << 22):
    """Chunked: the Planck SMICA map is 2.0 GB and a single .read() would add a
    2 GB bytes object on top of the downgraded map and the ~150 MB operator
    dict, under the 8 GiB `prlimit --as` cap."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(chunk), b""):
            h.update(blk)
    return h.hexdigest()


def jdump(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1, sort_keys=False, default=_default)
    log("wrote %s (%d bytes)" % (path, os.path.getsize(path)))


def _default(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(repr(o))


# ===========================================================================
# 1.  Quaternions, the Hurwitz units, D_4, and the homomorphism to SO(3)
# ===========================================================================
def quat_mul(a, b):
    """Hamilton product, (w, x, y, z)."""
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return (aw * bw - ax * bx - ay * by - az * bz,
            aw * bx + ax * bw + ay * bz - az * by,
            aw * by - ax * bz + ay * bw + az * bx,
            aw * bz + ax * by - ay * bx + az * bw)


def hurwitz_units():
    """The 24 units of the Hurwitz quaternion order = the 24 norm-1 vectors of
    the D_4 lattice = the binary tetrahedral group 2T.  Exact rational
    arithmetic (Fraction), so closure is checked exactly, not to a tolerance.

    8 of Lipschitz type (+-1, +-i, +-j, +-k) and 16 of type
    (+-1 +- i +- j +- k)/2.
    """
    q = []
    one = Fraction(1)
    zero = Fraction(0)
    for axis in range(4):
        for s in (one, -one):
            v = [zero] * 4
            v[axis] = s
            q.append(tuple(v))
    half = Fraction(1, 2)
    for signs in product((half, -half), repeat=4):
        q.append(tuple(signs))
    assert len(q) == 24
    return q


def d4_roots():
    """The 24 norm-2 vectors of D_4 = {x in Z^4 : sum x_i even}: all
    (+-1, +-1, 0, 0) and permutations.  These are the roots; |D_4 roots| = 24.
    """
    roots = []
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [0, 0, 0, 0]
                    v[i] = si
                    v[j] = sj
                    roots.append(tuple(v))
    assert len(roots) == 24
    return roots


def aut_of_root_system(roots, gram_dim=4):
    """|Aut(L)| by exhaustive search: count the linear maps that permute the
    root set and preserve the inner product.  A root-system automorphism is
    determined by the image of a basis of roots, so we enumerate images of one
    fixed basis and check the whole set maps to itself.

    NOTE (first run returned 384 for D_4).  The map must NOT be required to have
    integer entries in the ambient Z^n coordinates: the extra (triality)
    automorphisms of D_4 have HALF-integer entries -- they are exactly the
    left/right multiplications by the Hurwitz units, which preserve D_4 but not
    Z^4.  Requiring A in M_n(Z) computes Aut(Z^n) instead, which is 384 for
    n = 4 and is why the first run returned that number.  The condition kept is
    the correct one: A preserves the Gram matrix of the basis and maps the root
    set onto itself (which forces it to preserve the lattice the roots span).
    """
    R = np.array(roots, dtype=np.int64)
    n = R.shape[1]
    G = R @ R.T
    # pick a basis of the lattice spanned by the roots, from the roots
    basis_idx = []
    M = np.zeros((0, n), dtype=np.int64)
    for i in range(len(R)):
        M2 = np.vstack([M, R[i]])
        if np.linalg.matrix_rank(M2.astype(float)) > M.shape[0]:
            M = M2
            basis_idx.append(i)
        if M.shape[0] == gram_dim:
            break
    assert M.shape[0] == gram_dim, "roots do not span rank %d" % gram_dim
    B = M                                   # (r, n) basis of roots
    Bg = B @ B.T                            # gram of the basis
    rootset = set(map(tuple, R.tolist()))
    # candidate images: for each basis root, any root with the same norm
    cands = [[j for j in range(len(R)) if R[j] @ R[j] == B[k] @ B[k]]
             for k in range(B.shape[0])]
    count = 0
    Binv = np.linalg.pinv(B.astype(float))
    for choice in product(*cands):
        C = R[list(choice)]                 # (r, n) images
        if not np.array_equal(C @ C.T, Bg):
            continue
        # the unique linear map with B -> C on the span (may be half-integral)
        A = (Binv @ C.astype(float))        # acts on row vectors as x -> x @ A
        img = R @ A
        imgi = np.rint(img)
        if not np.allclose(img, imgi, atol=1e-9):
            continue
        if set(map(tuple, imgi.astype(np.int64).tolist())) == rootset:
            count += 1
    return count


def quat_to_so3(q):
    """phi(q) : v |-> q v qbar in the basis (i, j, k).  q must be a unit
    quaternion (w, x, y, z)."""
    w, x, y, z = [float(c) for c in q]
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - z * w), 2 * (x * z + y * w)],
        [2 * (x * y + z * w), 1 - 2 * (x * x + z * z), 2 * (y * z - x * w)],
        [2 * (x * z - y * w), 2 * (y * z + x * w), 1 - 2 * (x * x + y * y)],
    ])


def rot_angle(R):
    """Rotation angle in [0, pi] of an SO(3) matrix."""
    c = (np.trace(R) - 1.0) / 2.0
    return float(np.arccos(np.clip(c, -1.0, 1.0)))


def element_order(R, max_order=64):
    P = np.eye(3)
    for k in range(1, max_order + 1):
        P = P @ R
        if np.allclose(P, np.eye(3), atol=1e-9):
            return k
    return -1


def dedup_matrices(mats, tol=1e-9):
    out = []
    for M in mats:
        if not any(np.allclose(M, N, atol=tol) for N in out):
            out.append(M)
    return out


def group_A4():
    """The image phi(2T) in SO(3): the tetrahedral rotation group, order 12,
    isomorphic to A_4."""
    qs = hurwitz_units()
    mats = [quat_to_so3(q) for q in qs]
    return dedup_matrices(mats), mats, qs


def group_C12():
    """Control: the cyclic group of order 12 about the z axis.  Same order as
    the target, different group."""
    out = []
    for k in range(12):
        t = 2.0 * np.pi * k / 12.0
        out.append(np.array([[np.cos(t), -np.sin(t), 0.0],
                             [np.sin(t), np.cos(t), 0.0],
                             [0.0, 0.0, 1.0]]))
    return out


GROUPS = {"A4": lambda: group_A4()[0], "C12": group_C12}


def check_closure(mats, tol=1e-9):
    n = len(mats)
    for A in mats:
        for B in mats:
            P = A @ B
            if not any(np.allclose(P, M, atol=tol) for M in mats):
                return False
    return True


# ===========================================================================
# 2.  Wigner-D by eigendecomposition of J_y  (no factorials, stable to l=64)
# ===========================================================================
def jy_matrix(l):
    """J_y in the |l m> basis, m = -l..l at index m+l.  Hermitian, tridiagonal.

        <l,m+1|J_y|l,m> = -i/2 * sqrt(l(l+1) - m(m+1))
    """
    n = 2 * l + 1
    J = np.zeros((n, n), dtype=np.complex128)
    for m in range(-l, l):
        c = np.sqrt(l * (l + 1) - m * (m + 1))
        J[m + 1 + l, m + l] = -0.5j * c
        J[m + l, m + 1 + l] = +0.5j * c
    return J


class WignerCache:
    """Per-l eigendecomposition of J_y, plus the invariant-subspace basis of a
    finite group."""

    def __init__(self, lmax):
        self.lmax = lmax
        self.V = []
        self.lam = []
        for l in range(lmax + 1):
            w, v = np.linalg.eigh(jy_matrix(l))
            self.lam.append(w)
            self.V.append(v)
        self.m = [np.arange(-l, l + 1) for l in range(lmax + 1)]

    def D(self, l, alpha, beta, gamma):
        """D^l_{m'm}(a,b,c) = e^{-i m' a} d^l_{m'm}(b) e^{-i m c}."""
        V, lam, m = self.V[l], self.lam[l], self.m[l]
        d = (V * np.exp(-1j * beta * lam)[None, :]) @ V.conj().T
        return np.exp(-1j * m * alpha)[:, None] * d * np.exp(-1j * m * gamma)[None, :]

    def UD(self, l, U, alpha, beta, gamma):
        """U^dagger D^l(a,b,c), computed without ever forming d^l explicitly
        (cost 2 * dim(U) * (2l+1)^2 instead of (2l+1)^3)."""
        V, lam, m = self.V[l], self.lam[l], self.m[l]
        if U.shape[1] == 0:
            return np.zeros((0, 2 * l + 1), dtype=np.complex128)
        W = (U.conj().T * np.exp(-1j * m * alpha)[None, :]) @ V
        return (W * np.exp(-1j * beta * lam)[None, :]) @ \
               (V.conj().T * np.exp(-1j * m * gamma)[None, :])


def euler_zyz(R):
    """(alpha, beta, gamma) with R = Rz(alpha) Ry(beta) Rz(gamma)."""
    c = np.clip(R[2, 2], -1.0, 1.0)
    beta = float(np.arccos(c))
    s = np.sqrt(max(0.0, 1.0 - c * c))
    if s < 1e-12:
        gamma = 0.0
        if c > 0:
            alpha = float(np.arctan2(R[1, 0], R[0, 0]))
        else:
            alpha = float(np.arctan2(-R[1, 0], -R[0, 0]))
    else:
        alpha = float(np.arctan2(R[1, 2], R[0, 2]))
        gamma = float(np.arctan2(R[2, 1], -R[2, 0]))
    return alpha, beta, gamma


def character(l, theta):
    """chi_l(theta) = sin((l+1/2)theta)/sin(theta/2), chi_l(0) = 2l+1."""
    if abs(theta) < 1e-12:
        return float(2 * l + 1)
    return float(np.sin((l + 0.5) * theta) / np.sin(theta / 2.0))


def invariant_dims_by_character(mats, lmax):
    thetas = [rot_angle(M) for M in mats]
    return [int(round(sum(character(l, t) for t in thetas) / len(mats)))
            for l in range(lmax + 1)]


def invariant_basis(wc, mats, lmax, tol=1e-8):
    """Orthonormal basis U_l of the G-invariant subspace at each l, from the
    projector P = (1/|G|) sum_g D(g).  Also returns the projector diagnostics."""
    eulers = [euler_zyz(M) for M in mats]
    Us, dims, diag = [], [], []
    for l in range(lmax + 1):
        n = 2 * l + 1
        P = np.zeros((n, n), dtype=np.complex128)
        for (a, b, c) in eulers:
            P += wc.D(l, a, b, c)
        P /= len(mats)
        idem = float(np.max(np.abs(P @ P - P)))
        herm = float(np.max(np.abs(P - P.conj().T)))
        tr = complex(np.trace(P))
        w, v = np.linalg.eigh((P + P.conj().T) / 2.0)
        keep = w > 0.5
        U = v[:, keep]
        Us.append(U)
        dims.append(int(keep.sum()))
        diag.append({"l": l, "idempotency_max_abs": idem,
                     "hermiticity_max_abs": herm,
                     "trace_real": float(tr.real), "trace_imag": float(tr.imag),
                     "dim": int(keep.sum())})
    return Us, dims, diag


# ===========================================================================
# 3.  alm handling (healpy real-field convention <-> full m = -l..l)
# ===========================================================================
def alm_to_full(alm, lmax):
    """healpy alm (m >= 0) -> list of complex vectors a_l[m+l], m = -l..l,
    using a_{l,-m} = (-1)^m conj(a_{l,m})."""
    import healpy as hp
    out = []
    for l in range(lmax + 1):
        v = np.zeros(2 * l + 1, dtype=np.complex128)
        for m in range(0, l + 1):
            x = alm[hp.Alm.getidx(lmax, l, m)]
            v[m + l] = x
            if m > 0:
                v[-m + l] = ((-1) ** m) * np.conj(x)
        out.append(v)
    return out


def full_to_alm(full, lmax):
    import healpy as hp
    alm = np.zeros(hp.Alm.getsize(lmax), dtype=np.complex128)
    for l in range(lmax + 1):
        for m in range(0, l + 1):
            alm[hp.Alm.getidx(lmax, l, m)] = full[l][m + l]
    return alm


def band_norms(full, lmax):
    return np.array([float(np.sum(np.abs(full[l]) ** 2))
                     for l in range(lmax + 1)])


# ===========================================================================
# 4.  Orientation grid and the crystallinity operators
# ===========================================================================
def orientation_grid(n=N_ORIENT, seed=ORIENT_SEED):
    """n uniformly distributed rotations, frozen by seed.  Element 0 is the
    identity.  The SAME grid is used for data and for every null."""
    rng = np.random.default_rng(seed)
    mats = [np.eye(3)]
    while len(mats) < n:
        q = rng.normal(size=4)
        nq = np.linalg.norm(q)
        if nq < 1e-8:
            continue
        q = q / nq
        mats.append(quat_to_so3((q[0], q[1], q[2], q[3])))
    return mats


def build_operators(group, lmax=LMAX, n_orient=N_ORIENT):
    """Q[k][l] = U_l^dagger D^l(R_k^{-1}), the operator whose squared norm on
    a_l is the G-invariant power of the map seen with G placed at R_k."""
    mats = GROUPS[group]()
    wc = WignerCache(lmax)
    Us, dims, diag = invariant_basis(wc, mats, lmax)
    dims_char = invariant_dims_by_character(mats, lmax)
    grid = orientation_grid(n_orient)
    Q = []
    for k, R in enumerate(grid):
        a, b, c = euler_zyz(R.T)            # R^{-1} = R^T for SO(3)
        Q.append([wc.UD(l, Us[l], a, b, c) for l in range(lmax + 1)])
        if (k + 1) % 48 == 0:
            log("  operators: %d/%d orientations" % (k + 1, n_orient))
    return {"group": group, "mats": mats, "Q": Q, "Us": Us, "dims": dims,
            "dims_char": dims_char, "diag": diag, "grid": grid, "wc": wc}


def crystallinity(ops, full, bands=BANDS):
    """C[k, b] = invariant power fraction in band b at orientation k."""
    lmax = len(full) - 1
    nrm = band_norms(full, lmax)
    nk = len(ops["Q"])
    out = np.zeros((nk, len(bands)))
    for k in range(nk):
        Qk = ops["Q"][k]
        inv = np.zeros(lmax + 1)
        for l in range(lmax + 1):
            if Qk[l].shape[0] == 0:
                continue
            inv[l] = float(np.sum(np.abs(Qk[l] @ full[l]) ** 2))
        for bi, (l0, l1) in enumerate(bands):
            den = nrm[l0:l1 + 1].sum()
            out[k, bi] = inv[l0:l1 + 1].sum() / den if den > 0 else 0.0
    return out


def cmax(C):
    """C_b^max: max over the frozen orientation grid, per band."""
    return C.max(axis=0)


def z_from(cm, mean, std):
    return (cm - mean) / np.where(std > 0, std, np.inf)


def Z_stat(cm, mean, std):
    return float(np.max(z_from(cm, mean, std)))


def rank_p(data_stat, null_stats):
    n = len(null_stats)
    return (1.0 + int(np.sum(np.asarray(null_stats) >= data_stat))) / (n + 1.0)


# ===========================================================================
# 5.  Sky I/O (identical code path for data and nulls)
# ===========================================================================
def load_map_and_mask(which):
    import healpy as hp
    if which == "wmap":
        mpath, kpath = WMAP_MAP, WMAP_MASK
    elif which == "planck":
        mpath, kpath = PLANCK_MAP, PLANCK_MASK
    else:
        raise SystemExit("unknown map %r" % which)
    m = hp.read_map(mpath, field=0)
    k = hp.read_map(kpath, field=0)
    ns_m, ns_k = hp.get_nside(m), hp.get_nside(k)
    m_dg = hp.ud_grade(m, NSIDE_WORK)
    k_frac = hp.ud_grade(k.astype(np.float64), NSIDE_WORK)
    k_dg = (k_frac >= 0.5).astype(np.uint8)
    return m_dg, k_dg, {"map_path": mpath, "mask_path": kpath,
                        "native_nside_map": int(ns_m),
                        "native_nside_mask": int(ns_k),
                        "downgrade": ("hp.ud_grade to NSIDE_WORK=%d; mask "
                                      "ud_graded as float then thresholded at "
                                      ">= 0.5 to a binary mask; no extra "
                                      "smoothing" % NSIDE_WORK),
                        "fsky_work": float(k_dg.mean()),
                        "map_sha256": None}


def map_to_full(m, mask, lmax=LMAX):
    """THE code path.  Mask multiply -> mean subtract on unmasked -> map2alm ->
    zero l = 0, 1 -> full-m vectors.  Data and every sim go through this."""
    import healpy as hp
    x = np.asarray(m, dtype=np.float64).copy()
    sel = mask > 0
    x[~sel] = 0.0
    x[sel] -= x[sel].mean()
    x[~sel] = 0.0
    alm = hp.map2alm(x, lmax=lmax, iter=0)
    full = alm_to_full(alm, lmax)
    full[0][:] = 0.0
    if lmax >= 1:
        full[1][:] = 0.0
    return full


def estimate_cl(m, mask, lmax_cl):
    """pseudo-C_l / f_sky from the data's own masked map, with l = 0, 1 zeroed
    (following the v2 fix in cmb_tda.estimate_cl)."""
    import healpy as hp
    x = np.asarray(m, dtype=np.float64).copy()
    sel = mask > 0
    x[~sel] = 0.0
    x[sel] -= x[sel].mean()
    x[~sel] = 0.0
    fsky = float(mask.mean())
    cl = hp.anafast(x, lmax=lmax_cl)
    cl = cl / max(fsky, 1e-12)
    cl[0] = 0.0
    if len(cl) > 1:
        cl[1] = 0.0
    return cl


def gaussian_sim(cl, seed, nside=NSIDE_WORK):
    import healpy as hp
    np.random.seed(seed)
    return hp.synfast(cl, nside, lmax=len(cl) - 1, pixwin=False, verbose=False) \
        if "verbose" in hp.synfast.__code__.co_varnames \
        else hp.synfast(cl, nside, lmax=len(cl) - 1, pixwin=False)


# ===========================================================================
# 6.  Injection of a G-symmetric pattern
# ===========================================================================
def symmetric_pattern(ops, ells, seed, per_l_equal=False):
    """A random unit vector inside the G-invariant subspace at the given ells,
    expressed as full-m vectors (zero elsewhere).  Real-field reality is
    imposed by symmetrising; the result is re-projected so it stays exactly
    inside the invariant subspace."""
    rng = np.random.default_rng(seed)
    lmax = ops["wc"].lmax
    full = [np.zeros(2 * l + 1, dtype=np.complex128) for l in range(lmax + 1)]
    for l in ells:
        U = ops["Us"][l]
        if U.shape[1] == 0:
            continue
        c = rng.normal(size=U.shape[1]) + 1j * rng.normal(size=U.shape[1])
        v = U @ c
        # impose the real-field condition a_{l,-m} = (-1)^m conj(a_{l,m})
        mm = np.arange(-l, l + 1)
        v_conj = ((-1.0) ** mm) * np.conj(v[::-1])
        v = 0.5 * (v + v_conj)
        v = U @ (U.conj().T @ v)            # back into the invariant subspace
        if per_l_equal:
            nv = np.linalg.norm(v)
            if nv > 0:
                v = v / nv
        full[l] = v
    tot = np.sqrt(sum(float(np.sum(np.abs(full[l]) ** 2)) for l in range(lmax + 1)))
    if tot > 0:
        for l in range(lmax + 1):
            full[l] = full[l] / tot
    return full


def rotate_full(wc, full, R):
    a, b, c = euler_zyz(R)
    out = []
    for l in range(wc.lmax + 1):
        if not np.any(full[l]):
            out.append(full[l].copy())
        else:
            out.append(wc.D(l, a, b, c) @ full[l])
    return out


def random_rotation(rng):
    q = rng.normal(size=4)
    q = q / np.linalg.norm(q)
    return quat_to_so3((q[0], q[1], q[2], q[3]))


def inject(full_host, pattern, ells, f):
    """Add the pattern scaled so injected power = f * host power in the
    injection ells."""
    host_p = sum(float(np.sum(np.abs(full_host[l]) ** 2)) for l in ells)
    pat_p = sum(float(np.sum(np.abs(pattern[l]) ** 2)) for l in ells)
    if pat_p <= 0 or f <= 0:
        return [v.copy() for v in full_host]
    s = np.sqrt(f * host_p / pat_p)
    out = [v.copy() for v in full_host]
    for l in ells:
        out[l] = out[l] + s * pattern[l]
    return out


# ===========================================================================
# 7.  TDA path
# ===========================================================================
def import_cmb_tda():
    path = os.path.join(REVERSE_WT, CMB_TDA_REL)
    spec = importlib.util.spec_from_file_location("cmb_tda_for_lens", path)
    mod = importlib.util.module_from_spec(spec)
    assert mod.__name__ != "__main__"
    spec.loader.exec_module(mod)
    return mod, sha256_of(path)


def betti_curves_of(field, mask, nside, nu_grid, cmbmod):
    """b0, b1 via cmb_tda.build_topology + betti_curves_from_topology.

    b2 is NEVER used: simple_suite/report.json records that build_topology's
    complex has a hollow-tetrahedron defect giving spurious H2 classes, with
    "b0 and b1 are unaffected".  The curve b0 - b1 is NOT an Euler
    characteristic of that complex and is not called one here.
    """
    unmasked, edges, tris = cmbmod.build_topology(mask, nside)
    b0, b1, _ = cmbmod.betti_curves_from_topology(field, unmasked, edges, tris,
                                                  nu_grid, sublevel=True)
    return np.asarray(b0), np.asarray(b1)


def cubical_gnomonic(field, mask, nside, npix_side=256, reso_arcmin=None):
    """Independent b0/b1 cross-check: gudhi CubicalComplex lower-star
    persistence on a gnomonic projection centred on the north galactic pole."""
    import healpy as hp
    import gudhi
    if reso_arcmin is None:
        reso_arcmin = 60.0 * 60.0 / npix_side * 1.0
    x = np.asarray(field, dtype=np.float64).copy()
    x[mask <= 0] = np.nan
    proj = hp.gnomview(x, rot=(0, 90), xsize=npix_side, reso=reso_arcmin,
                       return_projected_map=True, no_plot=True)
    arr = np.ma.filled(np.asarray(proj, dtype=np.float64), np.nan)
    finite = np.isfinite(arr)
    arr = np.where(finite, arr, np.nanmax(arr[finite]) + 1.0)
    cc = gudhi.CubicalComplex(top_dimensional_cells=arr)
    cc.compute_persistence()
    b0 = cc.persistence_intervals_in_dimension(0)
    b1 = cc.persistence_intervals_in_dimension(1)
    return {"n_pixels": int(arr.size), "n_finite": int(finite.sum()),
            "n_H0_bars": int(len(b0)), "n_H1_bars": int(len(b1)),
            "reso_arcmin": float(reso_arcmin), "xsize": int(npix_side)}


# ===========================================================================
# 8.  Point-set path (3-D)
# ===========================================================================
def orbit_residuals(X, mats, kdt=None):
    """For each point x and each non-identity g in G, the distance from g.x to
    the nearest point of X.  A point set that is a union of G-orbits gives 0."""
    from scipy.spatial import cKDTree
    if kdt is None:
        kdt = cKDTree(X)
    res = []
    for M in mats:
        if np.allclose(M, np.eye(3), atol=1e-12):
            continue
        d, _ = kdt.query(X @ M.T, k=1)
        res.append(d)
    return np.concatenate(res)


def fold_points(X, mats):
    """G-fold: map every point to the representative of its orbit with the
    largest (z, y, x) key.  Persistence of the folded set is the point-set
    analogue of the residual r_G."""
    out = np.empty_like(X)
    for i, x in enumerate(X):
        orb = np.array([M @ x for M in mats])
        key = orb[:, 2] * 1e6 + orb[:, 1] * 1e3 + orb[:, 0]
        out[i] = orb[np.argmax(key)]
    return out


def alpha_betti(X, max_alpha_sq=np.inf):
    import gudhi
    ac = gudhi.AlphaComplex(points=X.tolist())
    st = ac.create_simplex_tree(max_alpha_square=max_alpha_sq)
    st.compute_persistence(homology_coeff_field=2)
    return {d: int(len(st.persistence_intervals_in_dimension(d)))
            for d in (0, 1, 2)}


# ===========================================================================
# STAGES
# ===========================================================================
def stage_groups(args):
    out = {"stage": "groups", "tier": "B (exact arithmetic with negative "
                                      "controls) for the group and character "
                                      "checks; X for the numerical projector "
                                      "diagnostics"}

    # --- Hurwitz units, exact closure -------------------------------------
    qs = hurwitz_units()
    qset = set(qs)
    closed = all(quat_mul(a, b) in qset for a in qs for b in qs)
    orders_q = {}
    for q in qs:
        p, k = q, 1
        while p != (Fraction(1), Fraction(0), Fraction(0), Fraction(0)) and k < 32:
            p = quat_mul(p, q)
            k += 1
        orders_q[k] = orders_q.get(k, 0) + 1
    out["hurwitz_units"] = {
        "count": len(qs), "closed_exactly_under_quaternion_product": bool(closed),
        "element_order_multiset": {str(k): v for k, v in sorted(orders_q.items())},
        "note": "the 24 units of the Hurwitz order = the 24 norm-1 vectors of "
                "D_4 = the binary tetrahedral group 2T"}

    # --- D_4 roots and |Aut(D_4)| with negative controls -------------------
    roots = d4_roots()
    t0 = time.time()
    aut_d4 = aut_of_root_system(roots, 4)
    t_d4 = time.time() - t0
    # negative controls: A_3 = D_3 (|Aut| = 48 = 2 * |W(A_3)| = 2*24) and the
    # 8 norm-1 vectors of Z^4 (|Aut(Z^4)| = 2^4 * 4! = 384)
    a3_roots = [tuple(list(r)[:3]) for r in roots if r[3] == 0]
    aut_a3 = aut_of_root_system(a3_roots, 3)
    z4 = []
    for i in range(4):
        for s in (1, -1):
            v = [0, 0, 0, 0]
            v[i] = s
            z4.append(tuple(v))
    aut_z4 = aut_of_root_system(z4, 4)
    out["D4_lattice"] = {
        "n_norm2_roots": len(roots), "n_norm1_units": 24,
        "Aut_order_exhaustive": int(aut_d4),
        "expected_from_project_track_F": 1152,
        "matches": bool(aut_d4 == 1152),
        "seconds": round(t_d4, 2),
        "negative_controls": {
            "A3_root_system": {"computed": int(aut_a3), "expected": 48,
                               "matches": bool(aut_a3 == 48)},
            "Z4_norm1_vectors_not_D4": {"computed": int(aut_z4), "expected": 384,
                                        "matches": bool(aut_z4 == 384)}},
        "control_note": "the same exhaustive routine returns 48 and 384 on the "
                        "controls, so 1152 is not produced by construction"}

    # --- the homomorphism to SO(3) ----------------------------------------
    mats, all_mats, _ = group_A4()
    preimages = []
    for M in mats:
        preimages.append(sum(1 for N in all_mats if np.allclose(M, N, atol=1e-9)))
    orders = {}
    for M in mats:
        o = element_order(M)
        orders[o] = orders.get(o, 0) + 1
    out["homomorphism"] = {
        "map": "phi(q) : v |-> q v qbar on purely imaginary quaternions, basis "
               "(i, j, k); kernel {+-1}",
        "n_quaternions_in": len(all_mats),
        "n_distinct_SO3_matrices_out": len(mats),
        "expected_out": 12,
        "double_cover_ok": bool(len(mats) == 12 and set(preimages) == {2}),
        "preimage_counts": sorted(set(preimages)),
        "closed_in_SO3": bool(check_closure(mats)),
        "element_order_multiset": {str(k): v for k, v in sorted(orders.items())},
        "expected_order_multiset": {"1": 1, "2": 3, "3": 8},
        "note": "order 12, NOT 24.  The control group is therefore C_12."}

    ctrl = group_C12()
    orders_c = {}
    for M in ctrl:
        o = element_order(M)
        orders_c[o] = orders_c.get(o, 0) + 1
    out["control_group"] = {
        "name": "C_12 (cyclic, about z)", "order": len(ctrl),
        "closed": bool(check_closure(ctrl)),
        "element_order_multiset": {str(k): v for k, v in sorted(orders_c.items())},
        "expected_order_multiset": {"1": 1, "2": 1, "3": 2, "4": 2, "6": 2, "12": 4}}

    # --- orbit structure on S^2 vs the Frame shapes ------------------------
    rng = np.random.default_rng(12345)
    sizes = {}
    for _ in range(2000):
        v = rng.normal(size=3)
        v /= np.linalg.norm(v)
        orb = [M @ v for M in mats]
        u = []
        for o in orb:
            if not any(np.allclose(o, w, atol=1e-8) for w in u):
                u.append(o)
        sizes[len(u)] = sizes.get(len(u), 0) + 1
    special = {}
    for M in mats:
        if element_order(M) == 1:
            continue
        w, v = np.linalg.eig(M)
        ax = np.real(v[:, np.argmin(np.abs(w - 1.0))])
        ax = ax / np.linalg.norm(ax)
        orb = [N @ ax for N in mats]
        u = []
        for o in orb:
            if not any(np.allclose(o, x, atol=1e-8) for x in u):
                u.append(o)
        special[len(u)] = special.get(len(u), 0) + 1
    out["orbit_structure_and_frame_shapes"] = {
        "generic_orbit_sizes_2000_random_directions": {str(k): v for k, v in sorted(sizes.items())},
        "orbit_sizes_at_rotation_axes": {str(k): v for k, v in sorted(special.items())},
        "frame_shapes_quoted": "LeanMaster v3.28.0 docs/VERIFIED_FOUNDATION.md:194-195: "
                               "\"The 192 elements of `(Z_2)^4 |x A_4` act as distinct "
                               "automorphisms of `H^2(Km A, Q)`. Their Frame shapes, from "
                               "Lefschetz numbers of powers, are `1^24, 1^8 2^8, 1^6 3^6, "
                               "1^4 2^2 4^4` (1, 27, 128, 36 elements): classes `1A, 2A, "
                               "3A, 4B`, all geometric.\"",
        "comparison_where_meaningful": {
            "1A_order_1": "matched: 1 identity in phi(2T)",
            "2A_order_2": "matched in element order: 3 elements of order 2 in phi(2T) "
                          "(the coordinate-axis pi rotations)",
            "3A_order_3": "matched in element order: 8 elements of order 3 in phi(2T) "
                          "(the body-diagonal 2pi/3 rotations)",
            "4B_order_4": "NOT MATCHED and cannot be.  phi(2T) = A_4 has no element of "
                          "order 4: the order-4 quaternion units +-i, +-j, +-k map under "
                          "phi to rotations by pi, i.e. order 2.  The order-4 structure "
                          "lives in the double cover 2T and in the (Z_2)^4 factor, "
                          "neither of which acts on S^2.  This is the boundary of where "
                          "the Frame-shape comparison is meaningful; it is stated, not "
                          "fudged."},
        "element_count_note": "the counts 1, 27, 128, 36 are of the order-192 group, not "
                              "of its A_4 quotient; phi(2T) has 1, 3, 8 elements of "
                              "orders 1, 2, 3.  No claim is made that these agree."}

    # --- representation theory: the pre-stated invariant dimensions ---------
    wc = WignerCache(LMAX)
    rep = {}
    for gname in ("A4", "C12"):
        gm = GROUPS[gname]()
        dchar = invariant_dims_by_character(gm, LMAX)
        Us, dnum, diag = invariant_basis(wc, gm, LMAX)
        rep[gname] = {
            "dims_by_character_l0_to_l64": dchar,
            "dims_by_projector_rank_l0_to_l64": dnum,
            "two_methods_agree": bool(dchar == dnum),
            "max_idempotency_residual": max(d["idempotency_max_abs"] for d in diag),
            "max_hermiticity_residual": max(d["hermiticity_max_abs"] for d in diag),
            "trace_equals_dim_max_error": max(abs(d["trace_real"] - d["dim"]) for d in diag),
        }
    rep["A4"]["prestated_l0_to_l6"] = [1, 0, 0, 1, 1, 0, 2]
    rep["A4"]["prestated_matches"] = bool(rep["A4"]["dims_by_character_l0_to_l64"][:7]
                                          == [1, 0, 0, 1, 1, 0, 2])
    rep["C12"]["prestated_rule"] = "d_l = 1 + 2*floor(l/12)"
    rep["C12"]["prestated_matches"] = bool(
        all(rep["C12"]["dims_by_character_l0_to_l64"][l] == 1 + 2 * (l // 12)
            for l in range(LMAX + 1)))
    out["representation_theory"] = rep

    # --- Wigner-D known-answer tests ---------------------------------------
    rng = np.random.default_rng(777)
    hom_err, unit_err, chi_err = 0.0, 0.0, 0.0
    for _ in range(20):
        R1, R2 = random_rotation(rng), random_rotation(rng)
        for l in (0, 1, 2, 3, 5, 8, 16, 33, 64):
            D1 = wc.D(l, *euler_zyz(R1))
            D2 = wc.D(l, *euler_zyz(R2))
            D12 = wc.D(l, *euler_zyz(R1 @ R2))
            hom_err = max(hom_err, float(np.max(np.abs(D1 @ D2 - D12))))
            unit_err = max(unit_err, float(np.max(np.abs(
                D1.conj().T @ D1 - np.eye(2 * l + 1)))))
            chi_err = max(chi_err, abs(complex(np.trace(D1)).real
                                       - character(l, rot_angle(R1))))
    out["wigner_known_answer"] = {
        "homomorphism_D(R1)D(R2)=D(R1R2)_max_abs_error": hom_err,
        "unitarity_max_abs_error": unit_err,
        "trace_D(R)=chi_l(theta_R)_max_abs_error": chi_err,
        "tolerance": 1e-8,
        "pass": bool(hom_err < 1e-8 and unit_err < 1e-8 and chi_err < 1e-8),
        "note": "trace(D^l(R)) = chi_l(theta_R) is convention-independent, so "
                "together with the homomorphism property it pins the "
                "representation to the rotation it is meant to represent."}

    # --- rotation known-answer against the pixel sphere --------------------
    try:
        import healpy as hp
        nside = 64
        rng2 = np.random.default_rng(2468)
        cl = np.zeros(9)
        cl[2:] = 1.0
        np.random.seed(2468)
        m = hp.synfast(cl, nside, lmax=8, pixwin=False)
        alm = hp.map2alm(m, lmax=8, iter=3)
        full = alm_to_full(alm, 8)
        R = random_rotation(rng2)
        rot_full = rotate_full(WignerCache(8), full, R)
        m_rot = hp.alm2map(full_to_alm(rot_full, 8), nside, lmax=8)
        th, ph = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)))
        n = np.array(hp.ang2vec(th, ph))
        ninv = n @ R                                  # R^{-1} n  (rows)
        m_ref64 = hp.get_interp_val(m, *hp.vec2ang(ninv))
        err64 = float(np.max(np.abs(m_rot - m_ref64)) / np.max(np.abs(m)))
        # the reference above is limited by HEALPix bilinear interpolation on an
        # nside-64 grid, not by the rotation.  Refine the REFERENCE (not the
        # tolerance): interpolate the same band-limited function on an nside-512
        # grid, where the bilinear error is ~64x smaller.
        m_hi = hp.alm2map(alm, 512, lmax=8)
        m_ref512 = hp.get_interp_val(m_hi, *hp.vec2ang(ninv))
        err512 = float(np.max(np.abs(m_rot - m_ref512)) / np.max(np.abs(m)))
        out["rotation_vs_pixel_sphere"] = {
            "test": "band-limited lmax=8 map, my D applied in harmonic space vs "
                    "bilinear interpolation of T at R^{-1}n",
            "tolerance": 5e-3,
            "max_rel_error_nside64_reference": err64,
            "max_rel_error_nside512_reference": err512,
            "max_rel_error": err512,
            "pass": bool(err512 < 5e-3),
            "note": "the tolerance 5e-3 was fixed a priori and is NOT relaxed. "
                    "The first attempt used an nside-64 interpolation reference "
                    "and gave 6.9e-3, i.e. at the level of HEALPix bilinear "
                    "interpolation error itself; the REFERENCE was then refined "
                    "to nside 512 (same band-limited function, ~64x smaller "
                    "interpolation error) and the test is judged at the same "
                    "5e-3. Both numbers are reported. The decisive convention "
                    "checks are the exact ones above: D(R1)D(R2) = D(R1R2) and "
                    "trace D^l(R) = chi_l(theta_R), both at ~1e-13; a wrong "
                    "convention would fail those by O(1)."}
    except Exception as e:                            # pragma: no cover
        out["rotation_vs_pixel_sphere"] = {"error": repr(e)}

    out["all_checks_pass"] = bool(
        out["hurwitz_units"]["closed_exactly_under_quaternion_product"]
        and out["D4_lattice"]["matches"]
        and out["D4_lattice"]["negative_controls"]["A3_root_system"]["matches"]
        and out["D4_lattice"]["negative_controls"]["Z4_norm1_vectors_not_D4"]["matches"]
        and out["homomorphism"]["double_cover_ok"]
        and out["homomorphism"]["closed_in_SO3"]
        and rep["A4"]["two_methods_agree"] and rep["A4"]["prestated_matches"]
        and rep["C12"]["two_methods_agree"] and rep["C12"]["prestated_matches"]
        and out["wigner_known_answer"]["pass"]
        and out["rotation_vs_pixel_sphere"].get("pass", False))
    out["peak_rss_mb"] = peak_rss_mb()
    jdump(out, os.path.join(HERE, "groups_check.json"))


def _ops_cache_path(group):
    return os.path.join(HERE, "operators_%s.npz" % group)


def stage_operators(args):
    for group in ("A4", "C12"):
        t0 = time.time()
        ops = build_operators(group)
        blob = {}
        for k in range(len(ops["Q"])):
            for l in range(LMAX + 1):
                if ops["Q"][k][l].shape[0]:
                    blob["Q_%d_%d" % (k, l)] = ops["Q"][k][l]
        blob["dims"] = np.array(ops["dims"])
        blob["dims_char"] = np.array(ops["dims_char"])
        np.savez_compressed(_ops_cache_path(group), **blob)
        log("%s: built in %.1f s, cached %s (%.1f MB), dims[0:7]=%s"
            % (group, time.time() - t0, _ops_cache_path(group),
               os.path.getsize(_ops_cache_path(group)) / 1e6, ops["dims"][:7]))


def load_ops(group):
    """Rebuild the operator structure, using the cache for Q if present."""
    ops = {"group": group}
    mats = GROUPS[group]()
    wc = WignerCache(LMAX)
    Us, dims, diag = invariant_basis(wc, mats, LMAX)
    ops.update({"mats": mats, "wc": wc, "Us": Us, "dims": dims,
                "dims_char": invariant_dims_by_character(mats, LMAX),
                "grid": orientation_grid()})
    p = _ops_cache_path(group)
    if os.path.exists(p):
        z = np.load(p)
        Q = []
        for k in range(N_ORIENT):
            row = []
            for l in range(LMAX + 1):
                key = "Q_%d_%d" % (k, l)
                row.append(z[key] if key in z
                           else np.zeros((0, 2 * l + 1), dtype=np.complex128))
            Q.append(row)
        ops["Q"] = Q
        log("loaded cached operators for %s" % group)
    else:
        ops["Q"] = build_operators(group)["Q"]
    return ops


def _null_path(which, group):
    return os.path.join(HERE, "null_%s_%s.npz" % (which, group))


def stage_nulls(args):
    import healpy as hp
    which = args.map
    m, mask, meta = load_map_and_mask(which)
    cl = estimate_cl(m, mask, 3 * NSIDE_WORK - 1)
    del m
    for group in ("A4", "C12"):
        ops = load_ops(group)
        CM = np.zeros((N_NULL, len(BANDS)))
        for i in range(N_NULL):
            sim = gaussian_sim(cl, NULL_BASE_SEED + i)
            full = map_to_full(sim, mask)
            CM[i] = cmax(crystallinity(ops, full))
            if (i + 1) % 25 == 0:
                log("%s %s: null %d/%d  rss=%.0f MB"
                    % (which, group, i + 1, N_NULL, peak_rss_mb()))
        np.savez_compressed(_null_path(which, group), cmax=CM, cl=cl)
        log("%s %s: nulls written, mean C_b^max = %s"
            % (which, group, np.round(CM.mean(axis=0), 6)))


def _null_Z(CM):
    """Leave-one-out standardised Z for each null realisation."""
    n = CM.shape[0]
    Z = np.zeros(n)
    for i in range(n):
        o = np.delete(CM, i, axis=0)
        Z[i] = Z_stat(CM[i], o.mean(axis=0), o.std(axis=0, ddof=1))
    return Z


def stage_injection(args):
    import healpy as hp
    which = args.map
    m, mask, meta = load_map_and_mask(which)
    cl = estimate_cl(m, mask, 3 * NSIDE_WORK - 1)
    del m
    ells = {"A4": [3, 4, 6], "C12": [12, 13]}
    ops = {g: load_ops(g) for g in ("A4", "C12")}
    nullZ, nullmu, nullsd = {}, {}, {}
    for g in ("A4", "C12"):
        CM = np.load(_null_path(which, g))["cmax"]
        nullZ[g] = _null_Z(CM)
        nullmu[g] = CM.mean(axis=0)
        nullsd[g] = CM.std(axis=0, ddof=1)
    patterns = {g: symmetric_pattern(ops[g], ells[g], INJ_PATTERN_SEED + i)
                for i, g in enumerate(("A4", "C12"))}
    results = []
    rng_master = np.random.default_rng(INJ_BASE_SEED)
    for inj_g in ("A4", "C12"):
        for f in AMPLITUDES:
            pw = {"A4": 0, "C12": 0}
            pvals = {"A4": [], "C12": []}
            for r in range(N_INJ):
                seed = INJ_BASE_SEED + 1000 * AMPLITUDES.index(f) \
                    + 100000 * (inj_g == "C12") + r
                sim = gaussian_sim(cl, seed)
                host = map_to_full(sim, mask)
                rng = np.random.default_rng(seed)
                Rr = random_rotation(rng)
                pat = rotate_full(ops[inj_g]["wc"], patterns[inj_g], Rr)
                fld = inject(host, pat, ells[inj_g], f)
                for test_g in ("A4", "C12"):
                    cm = cmax(crystallinity(ops[test_g], fld))
                    Z = Z_stat(cm, nullmu[test_g], nullsd[test_g])
                    p = rank_p(Z, nullZ[test_g])
                    pvals[test_g].append(p)
                    if p <= ALPHA:
                        pw[test_g] += 1
                if (r + 1) % 25 == 0:
                    log("inject %s f=%g: %d/%d  power(A4)=%d power(C12)=%d"
                        % (inj_g, f, r + 1, N_INJ, pw["A4"], pw["C12"]))
            results.append({
                "injected_group": inj_g, "injection_ells": ells[inj_g],
                "amplitude_fraction_of_band_power": f,
                "n_realisations": N_INJ, "alpha": ALPHA,
                "power_A4_statistic": pw["A4"] / N_INJ,
                "power_C12_statistic": pw["C12"] / N_INJ,
                "median_p_A4": float(np.median(pvals["A4"])),
                "median_p_C12": float(np.median(pvals["C12"]))})
            log("DONE inject %s f=%g -> power A4 %.2f, C12 %.2f"
                % (inj_g, f, pw["A4"] / N_INJ, pw["C12"] / N_INJ))
    out = {"stage": "injection", "tier": "X", "map_spectrum_from": which,
           "bands": BANDS, "lmax": LMAX, "n_orient": N_ORIENT,
           "n_null": N_NULL, "p_floor": 1.0 / (N_NULL + 1),
           "pattern_seed": INJ_PATTERN_SEED, "base_seed": INJ_BASE_SEED,
           "orientation_of_each_injection": "random, seeded per realisation, "
                                            "NOT drawn from the search grid",
           "results": results, "peak_rss_mb": peak_rss_mb()}
    # smallest amplitude with power >= 0.95 on the matched statistic
    sens = {}
    for g in ("A4", "C12"):
        hits = [r["amplitude_fraction_of_band_power"] for r in results
                if r["injected_group"] == g
                and r["power_%s_statistic" % g] >= 0.95]
        sens[g] = min(hits) if hits else None
    out["smallest_amplitude_with_power_0.95_at_alpha_0.05"] = sens
    jdump(out, os.path.join(HERE, "injection_power_%s.json" % which))


EXT_AMPLITUDES = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0]
BROAD_AMPLITUDES = [0.0, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2]


def stage_injection_ext(args):
    """Extended and broadband injection, ADDED AFTER the pre-registered ladder
    was run and BEFORE any real map was looked at.  Reason, recorded here:

      (i)  the pre-registered ladder (f <= 0.2, pattern confined to l = 3,4,6)
           produced no power at all -- power stayed at the alpha = 0.05 level
           everywhere -- so the required group-specificity demonstration was
           vacuous: showing that the CONTROL statistic does not fire is
           meaningless at an amplitude where the MATCHED statistic does not
           fire either.  Specificity has to be shown at an amplitude where the
           matched statistic reaches power >= 0.95.
      (ii) the pre-registered pattern draw put 0.007 / 0.214 / 0.779 of its
           power at l = 3 / 4 / 6, so a nominally three-multipole pattern was
           effectively single-multipole.  The extended runs normalise the
           pattern to EQUAL power per injection multipole.  That is a change to
           the SIGNAL MODEL, not to the test statistic; the statistic, the
           bands, the orientation grid, the nulls and the decision rule are
           exactly as pre-registered and are not touched.
      (iii) the task's own wording is "a sum of G-orbit-symmetrised spherical
           harmonics", i.e. broadband.  Confining the signal to 3 multipoles was
           a narrowing introduced here, and is the case where the statistic has
           least power.

    Nothing here is a data test and N_TESTS stays 4.
    """
    which = args.map
    m, mask, meta = load_map_and_mask(which)
    cl = estimate_cl(m, mask, 3 * NSIDE_WORK - 1)
    del m
    ops = {g: load_ops(g) for g in ("A4", "C12")}
    nullZ, nullmu, nullsd = {}, {}, {}
    for g in ("A4", "C12"):
        CM = np.load(_null_path(which, g))["cmax"]
        nullZ[g] = _null_Z(CM)
        nullmu[g] = CM.mean(axis=0)
        nullsd[g] = CM.std(axis=0, ddof=1)

    narrow = {"A4": [3, 4, 6], "C12": [12, 13]}
    broad = {"A4": list(range(2, LMAX + 1)), "C12": list(range(2, LMAX + 1))}
    scenarios = [("narrow_equal_power_per_l", narrow, EXT_AMPLITUDES),
                 ("broadband_l2_to_64", broad, BROAD_AMPLITUDES)]
    results = []
    for sc_name, ells_of, ladder in scenarios:
        pats = {g: symmetric_pattern(ops[g], ells_of[g], INJ_PATTERN_SEED + 7,
                                     per_l_equal=True) for g in ("A4", "C12")}
        for inj_g in ("A4", "C12"):
            for f in ladder:
                pw = {"A4": 0, "C12": 0}
                for r in range(N_INJ):
                    seed = (INJ_BASE_SEED + 500000 + 10000 * ladder.index(f)
                            + 200000 * (inj_g == "C12")
                            + 3000000 * scenarios.index(
                                (sc_name, ells_of, ladder)) + r)
                    sim = gaussian_sim(cl, seed)
                    host = map_to_full(sim, mask)
                    rng = np.random.default_rng(seed)
                    Rr = random_rotation(rng)
                    pat = rotate_full(ops[inj_g]["wc"], pats[inj_g], Rr)
                    fld = inject(host, pat, ells_of[inj_g], f)
                    for test_g in ("A4", "C12"):
                        cm = cmax(crystallinity(ops[test_g], fld))
                        Z = Z_stat(cm, nullmu[test_g], nullsd[test_g])
                        if rank_p(Z, nullZ[test_g]) <= ALPHA:
                            pw[test_g] += 1
                results.append({
                    "scenario": sc_name, "injected_group": inj_g,
                    "injection_ells": [ells_of[inj_g][0], ells_of[inj_g][-1]]
                    if len(ells_of[inj_g]) > 3 else ells_of[inj_g],
                    "n_injection_ells": len(ells_of[inj_g]),
                    "amplitude_fraction_of_injected_band_power": f,
                    "n_realisations": N_INJ, "alpha": ALPHA,
                    "power_A4_statistic": pw["A4"] / N_INJ,
                    "power_C12_statistic": pw["C12"] / N_INJ})
                log("DONE %s inject %s f=%g -> power A4 %.2f, C12 %.2f"
                    % (sc_name, inj_g, f, pw["A4"] / N_INJ, pw["C12"] / N_INJ))
    out = {"stage": "injection_ext", "tier": "X",
           "why_added": stage_injection_ext.__doc__,
           "added_before_any_real_map_was_read": True,
           "map_spectrum_from": which, "results": results}
    sens, spec = {}, {}
    for sc_name, _, _ in scenarios:
        for g in ("A4", "C12"):
            hits = [r["amplitude_fraction_of_injected_band_power"]
                    for r in results if r["scenario"] == sc_name
                    and r["injected_group"] == g
                    and r["power_%s_statistic" % g] >= 0.95]
            sens["%s/%s" % (sc_name, g)] = min(hits) if hits else None
            if hits:
                f0 = min(hits)
                other = "C12" if g == "A4" else "A4"
                row = [r for r in results if r["scenario"] == sc_name
                       and r["injected_group"] == g
                       and r["amplitude_fraction_of_injected_band_power"] == f0][0]
                spec["%s: inject %s at f=%g" % (sc_name, g, f0)] = {
                    "matched_statistic_%s_power" % g: row["power_%s_statistic" % g],
                    "control_statistic_%s_power" % other:
                        row["power_%s_statistic" % other],
                    "alpha": ALPHA}
    out["smallest_amplitude_with_power_0.95_at_alpha_0.05"] = sens
    out["group_specificity_at_a_detectable_amplitude"] = spec
    out["peak_rss_mb"] = peak_rss_mb()
    jdump(out, os.path.join(HERE, "injection_ext_%s.json" % which))


def stage_probe(args):
    """Fixed-orientation amplitude probe.  This exists so that the numbers cited
    in report.json for the linear-interference effect and for the uneven pattern
    draw come from a COMMITTED script and not from an inline `python -c`
    (LeanFlow CLAUDE.md section 6: no numbers without a script).

    It injects the PRE-REGISTERED pattern (uneven per-l draw, seed
    INJ_PATTERN_SEED) UNROTATED into one fixed Gaussian host (seed 1) and reads
    the crystallinity at grid element 0, which is the identity orientation, so
    the group is exactly where the pattern is.  Any failure to rise is then a
    property of the statistic, not a misalignment.
    """
    which = args.map
    m, mask, meta = load_map_and_mask(which)
    cl = estimate_cl(m, mask, 3 * NSIDE_WORK - 1)
    del m
    ops = load_ops("A4")
    ells = [3, 4, 6]
    pat = symmetric_pattern(ops, ells, INJ_PATTERN_SEED)
    pat_eq = symmetric_pattern(ops, ells, INJ_PATTERN_SEED, per_l_equal=True)
    host = map_to_full(gaussian_sim(cl, 1), mask)
    rows = []
    for f in [0.0, 0.05, 0.2, 1.0, 5.0, 50.0]:
        C = crystallinity(ops, inject(host, pat, ells, f))
        rows.append({"f": f, "C_at_identity_orientation": C[0].tolist(),
                     "C_b_max_over_grid": C.max(axis=0).tolist()})
        log("probe f=%g C[k=0]=%s" % (f, np.round(C[0], 4)))
    out = {"stage": "probe", "tier": "X", "why": stage_probe.__doc__,
           "host_sim_seed": 1, "pattern_seed": INJ_PATTERN_SEED,
           "injection_ells": ells,
           "pattern_power_per_l_preregistered_draw":
               {str(l): float(np.sum(np.abs(pat[l]) ** 2)) for l in ells},
           "pattern_power_per_l_equal_power_draw":
               {str(l): float(np.sum(np.abs(pat_eq[l]) ** 2)) for l in ells},
           "amplitude_scan_at_fixed_identity_orientation": rows,
           "analytic_null_level_band1": {
               "sum_d_l_over_band": int(sum(ops["dims"][2:9])),
               "sum_2l+1_over_band": int(sum(2 * l + 1 for l in range(2, 9))),
               "ratio": float(sum(ops["dims"][2:9])
                              / sum(2 * l + 1 for l in range(2, 9)))},
           "reading": "the machinery is correct -- at large f the crystallinity "
                      "at the identity orientation tends to 1 -- but at small f "
                      "it can FALL, because the cross term 2Re<P a, s p> "
                      "dominates s^2 and is destructive as often as "
                      "constructive.  That is the sensitivity floor of a "
                      "power-FRACTION statistic: it is set by the sample "
                      "variance of the invariant subspace's own random content."}
    out["peak_rss_mb"] = peak_rss_mb()
    jdump(out, os.path.join(HERE, "probe_%s.json" % which))


def stage_data(args):
    which = args.map
    m, mask, meta = load_map_and_mask(which)
    meta["map_sha256"] = sha256_of(meta["map_path"])
    meta["mask_sha256"] = sha256_of(meta["mask_path"])
    full = map_to_full(m, mask)
    res = {"stage": "data", "tier": "X", "map": which, "meta": meta,
           "bands": BANDS, "lmax": LMAX, "n_orient": N_ORIENT,
           "orientation_seed": ORIENT_SEED, "n_null": N_NULL,
           "null_base_seed": NULL_BASE_SEED,
           "p_floor": 1.0 / (N_NULL + 1),
           "alpha_family": ALPHA, "n_tests": N_TESTS,
           "alpha_bonferroni": ALPHA_BONF, "tests": {}}
    for group in ("A4", "C12"):
        ops = load_ops(group)
        CM = np.load(_null_path(which, group))["cmax"]
        mu, sd = CM.mean(axis=0), CM.std(axis=0, ddof=1)
        C = crystallinity(ops, full)
        cm = cmax(C)
        z = z_from(cm, mu, sd)
        Z = float(np.max(z))
        nullZ = _null_Z(CM)
        p = rank_p(Z, nullZ)
        res["tests"][group] = {
            "C_b_max_data": cm.tolist(),
            "C_b_max_null_mean": mu.tolist(),
            "C_b_max_null_std": sd.tolist(),
            "z_per_band": z.tolist(),
            "Z_max_over_bands": Z,
            "argmax_band": BANDS[int(np.argmax(z))],
            "best_orientation_index_per_band": C.argmax(axis=0).tolist(),
            "rank_p": p,
            "detection_at_bonferroni": bool(p <= ALPHA_BONF),
            "invariant_dims_l0_to_l16": ops["dims"][:17]}
        log("%s %s: Z=%.3f  p=%.4f  (Bonferroni threshold %.4f)"
            % (which, group, Z, p, ALPHA_BONF))
    res["peak_rss_mb"] = peak_rss_mb()
    jdump(res, os.path.join(HERE, "real_map_result_%s.json" % which))


def stage_tda(args):
    import healpy as hp
    which = args.map
    cmbmod, sha = import_cmb_tda()
    m, mask, meta = load_map_and_mask(which)
    full = map_to_full(m, mask)
    out = {"stage": "tda", "tier": "X", "map": which,
           "complex": "cmb_tda.build_topology (sha256 %s) -- b0 and b1 ONLY; "
                      "b2 is never used and b0-b1 is NOT called an Euler "
                      "characteristic, because simple_suite/report.json records "
                      "a hollow-tetrahedron defect in this complex that produces "
                      "spurious H2 classes while leaving b0 and b1 unaffected"
                      % sha,
           "tda_fixed_status": "audit/tda_validation/tda_fixed/ does not exist on "
                               "branch loop/tda-simple (git ls-tree -r --name-only "
                               "loop/tda-simple | grep -i tda_fixed -> empty); the "
                               "only directory there is simple_suite/",
           "curves": {}}
    nu = np.linspace(-4, 4, 41)
    for group in ("A4", "C12"):
        ops = load_ops(group)
        inv = []
        for l in range(LMAX + 1):
            U = ops["Us"][l]
            inv.append(U @ (U.conj().T @ full[l]) if U.shape[1] else
                       np.zeros_like(full[l]))
        res = [full[l] - inv[l] for l in range(LMAX + 1)]
        mg = hp.alm2map(full_to_alm(inv, LMAX), NSIDE_WORK, lmax=LMAX)
        rg = hp.alm2map(full_to_alm(res, LMAX), NSIDE_WORK, lmax=LMAX)
        t = hp.alm2map(full_to_alm(full, LMAX), NSIDE_WORK, lmax=LMAX)
        orth = float(abs(np.sum(mg[mask > 0] * rg[mask > 0]))
                     / (np.linalg.norm(mg[mask > 0]) * np.linalg.norm(rg[mask > 0])))
        b0m, b1m = betti_curves_of(mg, mask, NSIDE_WORK, nu, cmbmod)
        b0r, b1r = betti_curves_of(rg, mask, NSIDE_WORK, nu, cmbmod)
        b0t, b1t = betti_curves_of(t, mask, NSIDE_WORK, nu, cmbmod)
        out["curves"][group] = {
            "nu_grid": nu.tolist(),
            "b0_mG": b0m.tolist(), "b1_mG": b1m.tolist(),
            "b0_rG": b0r.tolist(), "b1_rG": b1r.tolist(),
            "b0_T": b0t.tolist(), "b1_T": b1t.tolist(),
            "power_fraction_in_mG": float(sum(np.sum(np.abs(inv[l]) ** 2)
                                              for l in range(LMAX + 1))
                                          / sum(np.sum(np.abs(full[l]) ** 2)
                                                for l in range(LMAX + 1))),
            "mG_rG_normalised_overlap_on_unmasked_pixels": orth,
            "max_b1_mG": int(b1m.max()), "max_b1_rG": int(b1r.max()),
            "max_b1_T": int(b1t.max())}
        try:
            out["curves"][group]["cubical_gnomonic_mG"] = \
                cubical_gnomonic(mg, mask, NSIDE_WORK)
            out["curves"][group]["cubical_gnomonic_rG"] = \
                cubical_gnomonic(rg, mask, NSIDE_WORK)
        except Exception as e:
            out["curves"][group]["cubical_gnomonic_error"] = repr(e)
        log("tda %s done" % group)
    out["gated"] = False
    out["note"] = "secondary, reported not gated (lens_spec.json section 11)"
    out["peak_rss_mb"] = peak_rss_mb()
    jdump(out, os.path.join(HERE, "tda_secondary_%s.json" % which))


def stage_pointcloud(args):
    """Known-answer test of the point-set path: a union of G-orbits must give
    ~0 orbit-distance residual and a G-folded set with |X|/|G| points; a random
    (control) set must give a large residual.  The 3-torus / cosmic-web
    application is NOT run here (see report.json)."""
    mats = GROUPS["A4"]()
    ctrl = GROUPS["C12"]()
    rng = np.random.default_rng(555)
    n_seed = 60
    seeds = rng.normal(size=(n_seed, 3))
    seeds /= np.linalg.norm(seeds, axis=1, keepdims=True)
    seeds *= rng.uniform(0.3, 1.0, size=(n_seed, 1))
    X_sym = np.vstack([np.array([M @ s for M in mats]) for s in seeds])
    X_rnd = rng.normal(size=X_sym.shape)
    X_rnd /= np.linalg.norm(X_rnd, axis=1, keepdims=True)
    X_rnd *= rng.uniform(0.3, 1.0, size=(X_rnd.shape[0], 1))
    out = {"stage": "pointcloud", "tier": "X",
           "construction": "%d random seed points, each replaced by its full "
                           "A_4 orbit -> %d points; control = the same number "
                           "of points drawn from the same radial/angular "
                           "distribution with no symmetry"
                           % (n_seed, X_sym.shape[0]),
           "seed": 555}
    for name, X in (("A4_orbit_union", X_sym), ("random_control", X_rnd)):
        r_t = orbit_residuals(X, mats)
        r_c = orbit_residuals(X, ctrl)
        fold = fold_points(X, mats)
        uniq = np.unique(np.round(fold, 6), axis=0)
        out[name] = {
            "n_points": int(X.shape[0]),
            "A4_orbit_residual_mean": float(r_t.mean()),
            "A4_orbit_residual_median": float(np.median(r_t)),
            "C12_orbit_residual_mean": float(r_c.mean()),
            "C12_orbit_residual_median": float(np.median(r_c)),
            "n_distinct_after_A4_folding": int(uniq.shape[0]),
            "expected_n_distinct_if_exact_orbits": int(X.shape[0] // len(mats)),
            "betti_of_folded_set": alpha_betti(uniq)}
    out["group_specific"] = {
        "target_residual_symmetric_over_random":
            out["A4_orbit_union"]["A4_orbit_residual_mean"]
            / max(out["random_control"]["A4_orbit_residual_mean"], 1e-12),
        "control_group_residual_on_symmetric_set_is_NOT_small":
            out["A4_orbit_union"]["C12_orbit_residual_mean"],
        "reading": "a A_4-symmetric point set has ~0 residual under A_4 and a "
                   "large residual under C_12; the statistic is group-specific"}
    out["not_run"] = ("the cosmic-web application (SDSS/DESI point set with the "
                      "CAMB lognormal mocks of audit/reverse_zero/"
                      "E5-cosmic-web-tda-scaled/ as null) was NOT run in this "
                      "session; only this synthetic known-answer test was")
    out["peak_rss_mb"] = peak_rss_mb()
    jdump(out, os.path.join(HERE, "pointcloud_kat.json"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True,
                    choices=["groups", "operators", "nulls", "injection",
                             "injection_ext", "probe", "data", "tda",
                             "pointcloud"])
    ap.add_argument("--map", default="wmap", choices=["wmap", "planck"])
    args = ap.parse_args()
    log("stage=%s map=%s python=%s host=%s"
        % (args.stage, args.map, sys.version.split()[0], platform.node()))
    {"groups": stage_groups, "operators": stage_operators,
     "nulls": stage_nulls, "injection": stage_injection,
     "injection_ext": stage_injection_ext, "probe": stage_probe,
     "data": stage_data,
     "tda": stage_tda, "pointcloud": stage_pointcloud}[args.stage](args)
    log("done, peak RSS %.0f MB" % peak_rss_mb())


if __name__ == "__main__":
    main()
