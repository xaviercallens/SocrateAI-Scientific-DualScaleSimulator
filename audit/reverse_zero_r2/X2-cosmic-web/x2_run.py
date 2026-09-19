#!/usr/bin/env python3
"""Ensemble runner.  cd audit/reverse_zero_r2/X2-cosmic-web && python x2_run.py MODE K0 K1 [--b B] [--workers 4]
MODE: biasgrid (k = 10*i + j, i<17: b=0.8+0.1 i, seed 6000000+100*i+j), fiducial (seed 7000000+k, topology),
poisson (7100000+k), x3 (7200000+k, b*sqrt3), control_void (7300000+k, topology), control_ring (7400000+k, topology).
One npz per mock in mocks/MODE_kkkkk.npz; existing files skipped (resumable)."""
import sys, os, time, pathlib, argparse
import numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import x2_lib as L
from scipy.stats import ks_2samp

OUT = L.HERE / "mocks"
_G = {}

def init():
    _G["c"] = L.Ctx(); _G["g"] = L.Gen(_G["c"])

def plant_voids(c, xyz, rng, radius=35.0, keep=0.05):
    ra, dec, r = L.xyz_to_rad(xyz)
    cen = []
    for t in range(4):
        idx = np.argwhere(c.tile_valid[t]); m = np.median(idx, 0)
        j = np.argmin(((idx - m) ** 2).sum(1)); cen.append(c.origin + (idx[j] + 0.5) * L.CELL)
    inside = np.zeros(len(xyz), bool)
    for cc in cen:
        inside |= ((xyz - cc) ** 2).sum(1) < radius ** 2
    kill = inside & (rng.uniform(size=len(xyz)) > keep)
    xk = xyz[~kill]
    # refill removed galaxies shell by shell, Poisson in volume outside the voids
    sh = np.clip(np.searchsorted(c.edges, np.sqrt((xyz[kill] ** 2).sum(1)), "right") - 1, 0, c.nsh - 1)
    need = np.bincount(sh, minlength=c.nsh)
    add = []
    while need.sum() > 0:
        p = c.draw_shell_points(need, rng)
        ok = np.ones(len(p), bool)
        for cc in cen:
            ok &= ((p - cc) ** 2).sum(1) >= radius ** 2
        p = p[ok]
        add.append(p)
        rr = np.sqrt((p ** 2).sum(1)); s2 = np.clip(np.searchsorted(c.edges, rr, "right") - 1, 0, c.nsh - 1)
        need = need - np.bincount(s2, minlength=c.nsh); need = np.maximum(need, 0)
    return np.vstack([xk] + add), int(kill.sum())

def plant_ring(c, xyz, rng, radius=55.0, tube=5.0, n_ring=9000):
    idx = np.argwhere(c.tile_valid[1]); m = np.median(idx, 0)
    j = np.argmin(((idx - m) ** 2).sum(1)); cen = c.origin + (idx[j] + 0.5) * L.CELL
    los = cen / np.linalg.norm(cen)
    e1 = np.cross(los, [0, 0, 1.0]); e1 /= np.linalg.norm(e1); e2 = np.cross(los, e1)
    ph = rng.uniform(0, 2 * np.pi, n_ring)
    ring = cen + radius * (np.cos(ph)[:, None] * e1 + np.sin(ph)[:, None] * e2) + rng.normal(0, tube, (n_ring, 3))
    ra, dec, r = L.xyz_to_rad(ring)
    ok = (ra >= 140) & (ra <= 220) & (dec >= 0) & (dec <= 50)
    ok &= c.maskpix[L.pix_of(np.where(ok, ra, 180.0), np.where(ok, dec, 25.0))]
    sh = np.searchsorted(c.edges, r, "right") - 1
    ok &= (sh >= 0) & (sh < c.nsh)
    ring, sh = ring[ok], sh[ok]
    rsh = np.searchsorted(c.edges, np.sqrt((xyz ** 2).sum(1)), "right") - 1
    drop = np.zeros(len(xyz), bool)
    for s in np.unique(sh):
        pool = np.nonzero((rsh == s) & ~drop)[0]
        n = min((sh == s).sum(), len(pool))
        drop[rng.choice(pool, n, replace=False)] = True
    return np.vstack([xyz[~drop], ring]), int(len(ring))

def one(job):
    mode, k, b = job
    f = OUT / ("%s_%05d.npz" % (mode, k))
    if f.exists():
        return str(f), 0.0
    c, g = _G["c"], _G["g"]
    t = time.time()
    topo = mode in ("fiducial", "control_void", "control_ring")
    extra = {}
    if mode == "biasgrid":
        i, j = divmod(k, 10); bb = round(0.8 + 0.1 * i, 6)
        seed = L.SEEDS["biasgrid"] + 100 * i + j; xyz = g.mock(bb, seed); extra["b"] = bb
    elif mode == "fiducial":
        seed = L.SEEDS["fiducial"] + k; xyz = g.mock(b, seed); extra["b"] = b
    elif mode == "x3":
        seed = L.SEEDS["x3"] + k; xyz = g.mock(b * np.sqrt(3.0), seed); extra["b"] = b * np.sqrt(3.0)
    elif mode == "poisson":
        seed = L.SEEDS["poisson"] + k; xyz = c.draw_shell_points(c.shell_cnt, np.random.default_rng(seed))
    elif mode == "control_void":
        seed = L.SEEDS["control_void"] + k; xyz = g.mock(b, seed)
        xyz, nk = plant_voids(c, xyz, np.random.default_rng(seed + 1)); extra["n_removed"] = nk
    elif mode == "control_ring":
        seed = L.SEEDS["control_ring"] + k; xyz = g.mock(b, seed)
        xyz, nr = plant_ring(c, xyz, np.random.default_rng(seed + 1)); extra["n_ring_in_footprint"] = nr
    else:
        raise SystemExit("mode")
    assert len(xyz) == c.N, (len(xyz), c.N)
    D = c.cgrid(xyz)
    xi = c.xi_from_grid(D)
    r = np.sqrt((xyz ** 2).sum(1))
    ksp = ks_2samp(r, c.r).pvalue
    res = dict(xi=xi, ks_r_p=ksp, short=getattr(g, "last_short", 0), attempts=getattr(g, "last_attempts", 0), seed=seed, **extra)
    if topo:
        res["betti"] = c.betti_curves(D)
    OUT.mkdir(exist_ok=True)
    np.savez(f, **res)
    return str(f), time.time() - t

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("mode"); ap.add_argument("k0", type=int); ap.add_argument("k1", type=int)
    ap.add_argument("--b", type=float, default=None); ap.add_argument("--workers", type=int, default=4)
    a = ap.parse_args()
    import multiprocessing as mp
    jobs = [(a.mode, k, a.b) for k in range(a.k0, a.k1) if not (OUT / ("%s_%05d.npz" % (a.mode, k))).exists()]
    print("todo", len(jobs), flush=True)
    t0 = time.time()
    with mp.get_context("fork").Pool(a.workers, initializer=init) as pool:
        for i, (f, dt) in enumerate(pool.imap_unordered(one, jobs)):
            if i % 10 == 0:
                print(i, f, round(dt, 1), round(time.time() - t0), flush=True)
