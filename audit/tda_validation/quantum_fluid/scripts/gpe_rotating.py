"""2-D rotating Gross-Pitaevskii ground state (imaginary time), test 2.

Units hbar = m = omega_trap = 1. Rotating-frame GPE:
  mu psi = [ -(1/2) nabla^2 + (1/2) r^2 + g |psi|^2 - Omega L_z ] psi,
  L_z = -i (x d/dy - y d/dx),  int |psi|^2 dx dy = 1.
Operator splitting (Bao, Wang & Markowich style, adapted to imaginary time):
  A_y = (1/2) p_y^2 - Omega x p_y   -> diagonal in (x, k_y):  (1/2) k_y^2 - Omega x k_y
  A_x = (1/2) p_x^2 + Omega y p_x   -> diagonal in (k_x, y):  (1/2) k_x^2 + Omega y k_x
  B   = (1/2) r^2 + g |psi|^2       -> diagonal in (x, y)
(-Omega L_z = -Omega x p_y + Omega y p_x). Strang step:
  e^{-B dt/2} e^{-A_x dt/2} e^{-A_y dt} e^{-A_x dt/2} e^{-B dt/2}, renormalise.
The imaginary-time propagator for A_y contains exp(+dt Omega x k_y), which is
bounded on the grid but can amplify high k_y at large |x|; this is harmless
because 0.5 k^2 dominates for |k| > 2 Omega |x| and the grid is chosen so
that 0.5 k_max^2 >> Omega x_max k_max. Energy/chem. potential are evaluated
with spectral derivatives.

Grid 384 x 384 on [-12, 12)^2, g = 1000, Omega in {0, 0.7, 0.8, 0.9}.
dt schedule: 0.005 (steps up to 8000) then 0.002 until |d mu| < 1e-7 per
100 steps or 12000 more steps. Initial state: Thomas-Fermi amplitude (with
omega' = sqrt(1-Omega^2)) times exp(i * smooth random phase) + small random
vortex seeds from numpy default_rng(7) (breaks the rotational symmetry so
vortices can nucleate in imaginary time).

Command:
  prlimit --as=8589934592 -- .venv-tda/bin/python gpe_rotating.py --omega 0.9
Output: DATA_ROOT/gpe/psi_Omega{Omega}.npz (psi, grid, mu history, params)
"""
import argparse
import json
import os
import time

import numpy as np

DATA_ROOT = "/mnt/disks/disk-socrateai-local-1/dualscale-data-r3/tda_validation/quantum_fluid"


def chem_potential(psi, X, Y, KX, KY, g, Om, dx):
    pk = np.fft.fft2(psi)
    lap = np.fft.ifft2(-(KX ** 2 + KY ** 2) * pk)
    dpx = np.fft.ifft2(1j * KX * pk)
    dpy = np.fft.ifft2(1j * KY * pk)
    Lz = -1j * (X * dpy - Y * dpx)
    Hpsi = -0.5 * lap + (0.5 * (X ** 2 + Y ** 2) + g * np.abs(psi) ** 2) * psi - Om * Lz
    mu = np.real(np.sum(np.conj(psi) * Hpsi)) * dx * dx
    e_int = 0.5 * g * np.sum(np.abs(psi) ** 4) * dx * dx
    return float(mu), float(mu - e_int)


def run(Om, g=1000.0, N=384, Lbox=24.0, seed=7, max1=8000, max2=12000):
    x = (np.arange(N) - N // 2) * (Lbox / N)
    dx = x[1] - x[0]
    X, Y = np.meshgrid(x, x, indexing="xy")  # X[j,i] = x_i, Y[j,i] = y_j ; axis0 = y, axis1 = x
    k = 2 * np.pi * np.fft.fftfreq(N, d=dx)
    KX, KY = np.meshgrid(k, k, indexing="xy")
    wp = np.sqrt(1 - Om ** 2)
    mu_tf = wp * np.sqrt(g / np.pi)
    n_tf = np.clip((mu_tf - 0.5 * wp ** 2 * (X ** 2 + Y ** 2)) / g, 0, None)
    rng = np.random.default_rng(seed)
    # smooth random phase + a few seeded vortices
    ph = np.real(np.fft.ifft2(np.fft.fft2(rng.normal(size=(N, N))) * np.exp(-(KX ** 2 + KY ** 2) / 2.0)))
    ph *= 2.0 / (ph.std() + 1e-12)
    psi = np.sqrt(n_tf + 1e-6) * np.exp(1j * ph)
    if Om > 0:
        R = np.sqrt(2 * mu_tf) / wp
        for _ in range(int(Om * R ** 2 * 0.5)):
            r0 = R * 0.8 * np.sqrt(rng.random()); a0 = 2 * np.pi * rng.random()
            x0, y0 = r0 * np.cos(a0), r0 * np.sin(a0)
            psi *= ((X - x0) + 1j * (Y - y0)) / np.sqrt((X - x0) ** 2 + (Y - y0) ** 2 + 0.05)
    psi /= np.sqrt(np.sum(np.abs(psi) ** 2) * dx * dx)
    V = 0.5 * (X ** 2 + Y ** 2)
    hist = []
    t0 = time.time()
    step_total = 0
    converged = False
    for dt, nmax in [(0.005, max1), (0.002, max2)]:
        # A_y in (x, k_y): fft along axis 0 (y). phase factor exp(-dt*(0.5 ky^2 - Om x ky))
        ky = k[:, None]; kx = k[None, :]
        Ay = np.exp(-dt * (0.5 * ky ** 2 - Om * x[None, :] * ky))        # [ky, x]
        Ax_half = np.exp(-0.5 * dt * (0.5 * kx ** 2 + Om * x[:, None] * kx))  # [y, kx]
        mu_prev = None
        for s in range(nmax):
            psi *= np.exp(-0.5 * dt * (V + g * np.abs(psi) ** 2))
            psi = np.fft.ifft(Ax_half * np.fft.fft(psi, axis=1), axis=1)
            psi = np.fft.ifft(Ay * np.fft.fft(psi, axis=0), axis=0)
            psi = np.fft.ifft(Ax_half * np.fft.fft(psi, axis=1), axis=1)
            psi *= np.exp(-0.5 * dt * (V + g * np.abs(psi) ** 2))
            psi /= np.sqrt(np.sum(np.abs(psi) ** 2) * dx * dx)
            step_total += 1
            if s % 100 == 99:
                mu, e_noint = chem_potential(psi, X, Y, KX, KY, g, Om, dx)
                hist.append([step_total, dt, mu])
                if dt == 0.002 and mu_prev is not None and abs(mu - mu_prev) < 1e-7:
                    converged = True
                    break
                mu_prev = mu
    mu, _ = chem_potential(psi, X, Y, KX, KY, g, Om, dx)
    out = os.path.join(DATA_ROOT, "gpe", f"psi_Omega{Om:.2f}.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez_compressed(out, psi=psi.astype(np.complex64), x=x, g=g, Omega=Om, N=N, Lbox=Lbox, seed=seed,
                        mu=mu, mu_hist=np.array(hist), converged=converged, steps=step_total)
    info = dict(Omega=Om, g=g, N=N, Lbox=Lbox, dx=float(dx), seed=seed, mu_final=mu, mu_TF_rot=float(mu_tf),
                converged=converged, steps=step_total, runtime_s=time.time() - t0,
                last_mu_changes=[h[2] for h in hist[-5:]], file=out)
    print(json.dumps(info))
    return info


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--omega", type=float, required=True)
    ap.add_argument("--max1", type=int, default=8000)
    ap.add_argument("--max2", type=int, default=12000)
    ap.add_argument("--N", type=int, default=384)
    a = ap.parse_args()
    run(a.omega, N=a.N, max1=a.max1, max2=a.max2)
