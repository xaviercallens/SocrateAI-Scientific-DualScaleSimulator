"""2-D XY model Monte Carlo (Wolff single-cluster, numba), test 1.

H = -J sum_<ij> cos(theta_i - theta_j), J = k_B = 1, periodic L x L.
Wolff 1989 embedding for O(2): pick a random reflection direction r = (cos phi,
sin phi); flipping a spin reflects it across the line perpendicular to r:
theta -> 2 phi + pi - theta. Bond i-j (both neighbours) is added to the
cluster with probability p = 1 - exp(min(0, -2 beta (s_i.r)(s_j.r))) where
s.r is evaluated BEFORE flipping. A sweep-equivalent = clusters until >= L^2
spins have been flipped (cumulative).

Observables per measurement (independent of TDA):
  e      energy per site
  nv     number of vortices+antivortices (plaquette winding, wrapped bond
         differences), and total winding (must be exactly 0 on a torus)
  hx, hy sum over x (y) bonds of cos(dtheta) ; sx, sy sum of sin(dtheta)
  Helicity modulus (per direction mu):
     Ups_mu = (1/L^2) [ <sum cos> - beta <(sum sin)^2> ]
  (second derivative of the free energy wrt a uniform twist, standard form;
  Upsilon -> 1 as T -> 0 and -> 0 at high T).
  mx, my magnetisation components.

Configurations are saved as float32 angles for later TDA (run with the
.venv-tda python, which has gudhi; this script needs only numpy+numba and is
run with /mnt/disks/disk-socrateai-local-1/venv-tdaval/bin/python).

Usage:
  prlimit --as=8589934592 -- venv-tdaval/bin/python xy_mc.py --L 64 --temps 0.40:1.60:0.05 \
      --ntherm 2000 --nmeas 400 --gap 2 --save-configs 1
Seeds: 20260919 + 1000*L + round(1000*T).
Output: DATA_ROOT/xy/L{L}/T{T:.3f}.npz (observables time series + configs)
"""
import argparse
import os
import time

import numpy as np
from numba import njit

DATA_ROOT = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/quantum_fluid"


@njit(cache=True)
def wolff_cluster(theta, L, beta, stack, incl):
    phi = np.random.random() * 2.0 * np.pi
    rx = np.cos(phi); ry = np.sin(phi)
    N = L * L
    s0 = np.random.randint(N)
    for k in range(N):
        incl[k] = False
    top = 0
    stack[top] = s0; top += 1
    incl[s0] = True
    size = 0
    while top > 0:
        top -= 1
        i = stack[top]
        # projection before flip
        pi_ = np.cos(theta[i]) * rx + np.sin(theta[i]) * ry
        theta[i] = 2.0 * phi + np.pi - theta[i]
        size += 1
        x = i % L; y = i // L
        for d in range(4):
            if d == 0:
                j = y * L + (x + 1) % L
            elif d == 1:
                j = y * L + (x - 1 + L) % L
            elif d == 2:
                j = ((y + 1) % L) * L + x
            else:
                j = ((y - 1 + L) % L) * L + x
            if not incl[j]:
                pj = np.cos(theta[j]) * rx + np.sin(theta[j]) * ry
                arg = -2.0 * beta * pi_ * pj
                if arg < 0.0:
                    p = 1.0 - np.exp(arg)
                    if np.random.random() < p:
                        incl[j] = True
                        stack[top] = j; top += 1
    return size


@njit(cache=True)
def sweep_eq(theta, L, beta, stack, incl):
    flipped = 0
    ncl = 0
    while flipped < L * L:
        flipped += wolff_cluster(theta, L, beta, stack, incl)
        ncl += 1
    return ncl


@njit(cache=True)
def seed_numba(s):
    np.random.seed(s)


def wrap(a):
    return (a + np.pi) % (2 * np.pi) - np.pi


def measure(theta2d):
    t = theta2d
    dx = np.roll(t, -1, axis=1) - t
    dy = np.roll(t, -1, axis=0) - t
    cx, cy = np.cos(dx).sum(), np.cos(dy).sum()
    sx, sy = np.sin(dx).sum(), np.sin(dy).sum()
    L = t.shape[0]
    e = -(cx + cy) / (L * L)
    t10 = np.roll(t, -1, axis=1); t11 = np.roll(t10, -1, axis=0); t01 = np.roll(t, -1, axis=0)
    w = wrap(t10 - t) + wrap(t11 - t10) + wrap(t01 - t11) + wrap(t - t01)
    q = np.rint(w / (2 * np.pi)).astype(np.int64)
    return e, cx, cy, sx, sy, int(np.abs(q).sum()), int(q.sum()), np.cos(t).sum() / (L * L), np.sin(t).sum() / (L * L)


def iat(x, cmax=200):
    x = np.asarray(x, float) - np.mean(x)
    v = x.var()
    if v == 0:
        return 0.5
    tau = 0.5
    for k in range(1, min(cmax, len(x) // 4)):
        c = np.mean(x[:-k] * x[k:]) / v
        if c <= 0:
            break
        tau += c
    return float(tau)


def run_one(L, T, ntherm, nmeas, gap, save_configs, outdir):
    seed = 20260919 + 1000 * L + int(round(1000 * T))
    rng = np.random.default_rng(seed)
    seed_numba(seed % (2 ** 31))
    theta = rng.uniform(-np.pi, np.pi, L * L)
    beta = 1.0 / T
    stack = np.zeros(L * L, dtype=np.int64)
    incl = np.zeros(L * L, dtype=np.bool_)
    t0 = time.time()
    for _ in range(ntherm):
        sweep_eq(theta, L, beta, stack, incl)
    obs = np.zeros((nmeas, 9))
    cfgs = np.zeros((nmeas, L, L), dtype=np.float32) if save_configs else None
    for m in range(nmeas):
        for _ in range(gap):
            sweep_eq(theta, L, beta, stack, incl)
        theta = wrap(theta)
        t2 = theta.reshape(L, L)
        obs[m] = measure(t2)
        assert obs[m, 6] == 0, "total winding on the torus must be exactly 0"
        if save_configs:
            cfgs[m] = t2.astype(np.float32)
    e, cx, cy, sx, sy, nv = obs[:, 0], obs[:, 1], obs[:, 2], obs[:, 3], obs[:, 4], obs[:, 5]
    N = L * L
    ups_x = (cx.mean() - beta * np.mean(sx ** 2)) / N
    ups_y = (cy.mean() - beta * np.mean(sy ** 2)) / N
    # jackknife error (20 blocks) for Upsilon and vortex density
    nb = 20
    blk = np.array_split(np.arange(nmeas), nb)
    ups_j, rho_j = [], []
    for b in blk:
        k = np.setdiff1d(np.arange(nmeas), b)
        ux = (cx[k].mean() - beta * np.mean(sx[k] ** 2)) / N
        uy = (cy[k].mean() - beta * np.mean(sy[k] ** 2)) / N
        ups_j.append(0.5 * (ux + uy)); rho_j.append(nv[k].mean() / N)
    ups_j = np.array(ups_j); rho_j = np.array(rho_j)
    jk = lambda a: float(np.sqrt((nb - 1) / nb * np.sum((a - a.mean()) ** 2)))
    summary = dict(L=L, T=T, seed=seed, ntherm=ntherm, nmeas=nmeas, gap=gap,
                   energy=float(e.mean()), upsilon=float(0.5 * (ups_x + ups_y)), upsilon_err=jk(ups_j),
                   vortex_density=float(nv.mean() / N), vortex_density_err=jk(rho_j),
                   magnetisation=float(np.mean(np.hypot(obs[:, 7], obs[:, 8]))),
                   tau_int_energy_meas_units=iat(e), runtime_s=time.time() - t0)
    fn = os.path.join(outdir, f"T{T:.3f}.npz")
    if save_configs:
        np.savez_compressed(fn, obs=obs, configs=cfgs, **{k: v for k, v in summary.items()})
    else:
        np.savez_compressed(fn, obs=obs, **{k: v for k, v in summary.items()})
    return summary, fn


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--L", type=int, required=True)
    ap.add_argument("--temps", type=str, required=True, help="start:stop:step inclusive, or comma list")
    ap.add_argument("--ntherm", type=int, default=2000)
    ap.add_argument("--nmeas", type=int, default=400)
    ap.add_argument("--gap", type=int, default=2)
    ap.add_argument("--save-configs", type=int, default=1)
    a = ap.parse_args()
    if ":" in a.temps:
        s, e, st = map(float, a.temps.split(":"))
        temps = np.round(np.arange(s, e + st / 2, st), 4)
    else:
        temps = [float(x) for x in a.temps.split(",")]
    outdir = os.path.join(DATA_ROOT, "xy", f"L{a.L}")
    os.makedirs(outdir, exist_ok=True)
    for T in temps:
        summ, fn = run_one(a.L, float(T), a.ntherm, a.nmeas, a.gap, a.save_configs, outdir)
        print(f"L={a.L} T={T:.3f} Ups={summ['upsilon']:.4f}+-{summ['upsilon_err']:.4f} "
              f"rho_v={summ['vortex_density']:.5f} m={summ['magnetisation']:.3f} tau_E={summ['tau_int_energy_meas_units']:.2f} "
              f"t={summ['runtime_s']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
