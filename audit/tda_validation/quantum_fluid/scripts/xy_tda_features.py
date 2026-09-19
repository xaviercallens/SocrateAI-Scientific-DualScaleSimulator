"""Test 1 TDA features from saved XY configurations, using the pipeline under test.

Per configuration (lower-star path, cmb_tda.betti_curves_from_topology):
  field theta' = wrap(theta - arg(M)), M = sum exp(i theta) (gauge fix, stated
  in expectations.json), on the Freudenthal torus (qf_common.torus_topology,
  V-E+F = 0 asserted). Sublevel Betti curves b0, b1 at 33 absolute thresholds
  linspace(-pi, pi, 33), passed to the function in its sigma units
  (nu = threshold / std(theta')), since the function divides by that std.
  Vortex count N_v from wrapped plaquette winding (independent observable).
Negative control (--shuffle-n K): for the first K configurations per T, the
  same field with sites permuted by numpy default_rng(12345 + config index).
Secondary (--alpha-n K, alpha path, cosmic_web_tda_scaled.alpha_persistence):
  combined vortex+antivortex positions (plaquette centres, non-periodic),
  max_alpha_sq = 100; f_short = fraction of finite H0 deaths <= 0.75; random
  null = same number of points placed without replacement on plaquette
  centres, numpy default_rng(777 + config index).

Usage:
  prlimit --as=8589934592 -- .venv-tda/bin/python xy_tda_features.py --L 64 --stride 2 \
      --shuffle-n 50 --alpha-n 100 [--temps 0.40,0.45]
Output: DATA_ROOT/xy_tda/L{L}/T{T:.3f}.npz
"""
import argparse
import glob
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import qf_common as qc  # noqa: E402

THR = np.linspace(-np.pi, np.pi, 33)


def betti(theta2d, topo, cmb):
    u, e, t = topo
    f = theta2d.ravel().astype(np.float64)
    sig = f.std()
    b0, b1, _ = cmb.betti_curves_from_topology(f, u, e, t, THR / sig, sublevel=True)
    return b0, b1


def alpha_fshort(pts, web):
    info, by_dim = web.alpha_persistence(pts, 100.0, "xy_vortex", None)
    d0 = qc.h0_finite_deaths(by_dim)
    if d0.size == 0:
        return np.nan, 0
    return float(np.mean(d0 <= 0.75)), int(d0.size)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--stride", type=int, default=2)
    ap.add_argument("--shuffle-n", type=int, default=0)
    ap.add_argument("--alpha-n", type=int, default=0)
    ap.add_argument("--temps", type=str, default=None)
    a = ap.parse_args()
    cmb, web = qc.cmb(), qc.web()
    L = a.L
    topo = qc.torus_topology(L)
    src = os.path.join(qc.DATA_ROOT, "xy", f"L{L}")
    dst = os.path.join(qc.DATA_ROOT, "xy_tda", f"L{L}")
    os.makedirs(dst, exist_ok=True)
    files = sorted(glob.glob(os.path.join(src, "T*.npz")))
    if a.temps:
        want = {f"T{float(x):.3f}.npz" for x in a.temps.split(",")}
        files = [f for f in files if os.path.basename(f) in want]
    for fn in files:
        out = os.path.join(dst, os.path.basename(fn))
        if os.path.exists(out):
            continue
        t0 = time.time()
        d = np.load(fn)
        T = float(d["T"])
        cfgs = d["configs"][:: a.stride]
        n = len(cfgs)
        B0 = np.zeros((n, THR.size), int); B1 = np.zeros((n, THR.size), int); NV = np.zeros(n, int)
        nsh = min(a.shuffle_n, n)
        S0 = np.zeros((nsh, THR.size), int); S1 = np.zeros((nsh, THR.size), int)
        na = min(a.alpha_n, n)
        FS = np.full(na, np.nan); FR = np.full(na, np.nan); NA = np.zeros(na, int)
        for i, th in enumerate(cfgs):
            th = th.astype(np.float64)
            M = np.exp(1j * th).sum()
            thp = qc.wrap(th - np.angle(M))
            q = qc.plaquette_winding(thp)
            assert q.sum() == 0
            NV[i] = np.abs(q).sum()
            B0[i], B1[i] = betti(thp, topo, cmb)
            if i < nsh:
                perm = np.random.default_rng(12345 + i).permutation(L * L)
                S0[i], S1[i] = betti(thp.ravel()[perm].reshape(L, L), topo, cmb)
            if i < na:
                yy, xx = np.nonzero(q)
                NA[i] = yy.size
                if yy.size >= 4:
                    pts = np.stack([xx + 0.5, yy + 0.5], 1).astype(float)
                    FS[i], _ = alpha_fshort(pts, web)
                    sel = np.random.default_rng(777 + i).choice(L * L, size=yy.size, replace=False)
                    rp = np.stack([sel % L + 0.5, sel // L + 0.5], 1).astype(float)
                    FR[i], _ = alpha_fshort(rp, web)
        np.savez_compressed(out, T=T, L=L, thr=THR, b0=B0, b1=B1, nv=NV, b0_shuf=S0, b1_shuf=S1,
                            fshort=FS, fshort_rand=FR, n_alpha_pts=NA, stride=a.stride)
        print(f"L={L} T={T:.3f} n={n} <nv>={NV.mean():.1f} b0(0)={B0[:,16].mean():.1f} "
              f"shuf b0(0)={S0[:,16].mean() if nsh else float('nan'):.1f} t={time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
