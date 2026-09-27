#!/usr/bin/env python3
"""L5 experiment: does the K3 level N reach a reported observable through the dynamics?

Spec L5 (specs/LEANFLOW_ARCHITECTURE.md): two runs at different N must give a stated
observable that differs by more than its numerical tolerance, and the difference must
survive a seed change.

This is a measurement, not an adoption. The production integrator is unchanged:
`workshopcosmo.py` keeps its hand-built double well, and nothing here alters its
defaults. The Einstein-Klein-Gordon system and the CPL extraction are the production
ones (`cosmology_rhs` equations and `run_observables_analysis`); only V is swapped for
the opt-in Gamma_0(N)+ potential of `leanflow.core.modular_potential`.

Held fixed across N, and therefore NOT derived: the overall normalisation of V (it is
a modular function in [0, 1], not a physical potential), `scale` in log mode, RHO_M0
and RHO_R0, the initial scale factor, and t_max.

"Seed change" means: the initial tau is jittered by a seeded random offset. The
integration itself is deterministic. "Numerical tolerance" is measured, not assumed,
as the change in each observable between rtol=1e-8 and rtol=1e-10.

Usage: python3 scripts/l5_level_experiment.py [--quick] [--out DIR]
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
import time
from multiprocessing import Pool

import numpy as np
from scipy.integrate import solve_ivp

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

import workshopcosmo as wc  # noqa: E402
from leanflow.core.modular_potential import modular_potential_xy  # noqa: E402

LEVELS = (7, 12)
MODES = ("ratio", "log")
FD_STEP = 1e-7
Y_FLOOR = 1e-8
# mpmath forms q = exp(2 pi i N tau) exactly (arbitrary precision), but jtheta's
# to_fixed step needs a bounded number of bits, so |Im tau| must stay finite for
# evaluation. Both modes are documented to approach their y->infinity asymptote
# (ratio: 1.0; log: w/(scale+w) -> 1.0) well before y=60, so clamping there changes
# no reachable value -- it only bounds what Radau's Jacobian probes can ask for.
Y_CAP = 60.0
CLAMP_EVENTS: list = []  # (job-local) real trajectory points that needed clamping, not FD probes


def _bounded(x: float, y: float) -> tuple[float, float]:
    xe = math.fmod(x, 1.0) if math.isfinite(x) else 0.0
    ye = y if math.isfinite(y) else Y_CAP
    ye = min(max(ye, Y_FLOOR), Y_CAP)
    return xe, ye


def potential(x: float, y: float, N: int, mode: str):
    """V and its central-difference gradient. mpmath at dps 25 makes h=1e-7 safe.

    x is reduced mod 1 before evaluation (exact: T=[[1,1],[0,1]] is in Gamma_0(N),
    so V(x+1,y)=V(x,y)) and y is clamped to [Y_FLOOR, Y_CAP] (an evaluation-domain
    bound, see Y_CAP). Both guard the same failure: with zero force in the flat
    region, Radau's numerical Jacobian enlarges its probe step without bound.
    """
    if not math.isfinite(y) or y < Y_FLOOR or y > Y_CAP:
        CLAMP_EVENTS.append((x, y))  # the *state*, not an FD probe, left the evaluation domain
    f = lambda a, b: modular_potential_xy(*_bounded(a, b), N=N, mode=mode)
    v = f(x, y)
    dvx = (f(x + FD_STEP, y) - f(x - FD_STEP, y)) / (2 * FD_STEP)
    dvy = (f(x, y + FD_STEP) - f(x, y - FD_STEP)) / (2 * FD_STEP)
    return v, dvx, dvy


def rhs(_t, s, N, mode):
    # Same equations as wc.cosmology_rhs; only the potential differs.
    a = max(float(s[0]), 1e-20)
    x, y = float(s[1]), max(float(s[2]), Y_FLOOR)
    u, v = float(s[3]), float(s[4])
    vp, dvx, dvy = potential(x, y, N, mode)
    t_kin = (u * u + v * v) / (2 * y * y)
    a3 = a ** 3
    h = math.sqrt(max((t_kin + vp + wc.RHO_M0 / a3 + wc.RHO_R0 / (a3 * a)) / 3.0, 0.0))
    return [a * h, u, v,
            (2 / y) * u * v - 3 * h * u - y * y * dvx,
            (u * u - v * v) / y - 3 * h * v - y * y * dvy]


def start_point(start: str, N: int) -> tuple[float, float]:
    if start == "common":
        return 0.2, 1.2  # a generic point, special to neither level
    return 0.001, 1.0 / math.sqrt(N) + 0.001  # the level's own self-dual point, offset as in production


def run_one(job: dict) -> dict:
    N, mode, start, seed, rtol, t_max = (job[k] for k in ("N", "mode", "start", "seed", "rtol", "t_max"))
    x0, y0 = start_point(start, N)
    if seed:
        rng = np.random.default_rng(seed)
        x0 += float(rng.uniform(-1e-3, 1e-3))
        y0 += float(rng.uniform(-1e-3, 1e-3))
    t0 = time.time()
    CLAMP_EVENTS.clear()
    try:
        sol = solve_ivp(rhs, (0.0, t_max), [1e-10, x0, y0, 0.0, 0.0], method="Radau",
                        t_eval=np.linspace(0.0, t_max, 500), rtol=rtol, atol=rtol * 1e-2,
                        first_step=1e-22, max_step=0.5, args=(N, mode))
    except Exception as exc:  # recorded as a failed run, never dropped
        return dict(job, x0=x0, y0=y0, success=False, message=f"{type(exc).__name__}: {exc}"[:300],
                    seconds=round(time.time() - t0, 1), clamp_events=len(CLAMP_EVENTS))
    out = dict(job, x0=x0, y0=y0, success=bool(sol.success), message=sol.message,
               nfev=int(sol.nfev), seconds=round(time.time() - t0, 1), clamp_events=len(CLAMP_EVENTS))
    if not sol.success or sol.y.shape[1] < 10:
        return out

    a_arr, x_arr, y_arr, u_arr, v_arr = sol.y
    w_phi, omega_phi = [], []
    for a, x, y, u, v in zip(a_arr, x_arr, y_arr, u_arr, v_arr):
        a, y = max(a, 1e-20), max(y, Y_FLOOR)
        vp = modular_potential_xy(x, y, N=N, mode=mode)
        t_kin = (u * u + v * v) / (2 * y * y)
        rho_phi = t_kin + vp
        rho_tot = rho_phi + wc.RHO_M0 / a ** 3 + wc.RHO_R0 / a ** 4
        w_phi.append((t_kin - vp) / max(rho_phi, 1e-15))
        omega_phi.append(rho_phi / max(rho_tot, 1e-15))

    obs = wc.run_observables_analysis({"a": a_arr.tolist(), "w_phi": w_phi})
    out.update(
        a_final=float(a_arr[-1]),
        tau_final=[float(x_arr[-1]), float(y_arr[-1])],
        w_final=float(w_phi[-1]),
        omega_phi_final=float(omega_phi[-1]),
        w0=obs["w0_fit"], wa=obs["wa_fit"],
        n_fit_points=int(np.sum((a_arr >= 0.2) & (a_arr <= 1.0))),
    )
    return out


OBSERVABLES = ("w0", "wa", "w_final", "omega_phi_final")


def verdicts(rows: list[dict]) -> list[dict]:
    """Per (mode, start, observable): is |obs(7) - obs(12)| > tolerance for every seed?"""
    key = lambda r: (r["mode"], r["start"], r["seed"], r["rtol"], r["N"])
    by = {key(r): r for r in rows if r.get("success") and "w0" in r}
    seeds = sorted({r["seed"] for r in rows})
    out = []
    for mode in MODES:
        for start in ("common", "self_dual"):
            for ob in OBSERVABLES:
                per_seed, ok_all, missing = [], True, False
                for s in seeds:
                    try:
                        lo = {N: by[(mode, start, s, 1e-8, N)][ob] for N in LEVELS}
                        hi = {N: by[(mode, start, s, 1e-10, N)][ob] for N in LEVELS}
                    except KeyError:
                        missing = True
                        continue
                    diff = hi[12] - hi[7]
                    tol = max(abs(hi[N] - lo[N]) for N in LEVELS)
                    passes = abs(diff) > tol
                    ok_all &= passes
                    per_seed.append(dict(seed=s, diff_12_minus_7=diff, tolerance=tol, exceeds=passes))
                signs = {math.copysign(1, p["diff_12_minus_7"]) for p in per_seed if p["diff_12_minus_7"]}
                out.append(dict(mode=mode, start=start, observable=ob, per_seed=per_seed,
                                incomplete=missing,
                                meets_L5=bool(per_seed) and ok_all and not missing and len(signs) == 1))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="1 seed, rtol 1e-8 only, t_max 20: a timing smoke test")
    ap.add_argument("--out", default=os.path.join(ROOT, "audit", "l5_level_runs"))
    ap.add_argument("--procs", type=int, default=6)
    args = ap.parse_args()

    seeds = [0] if args.quick else [0, 1, 2, 3]
    rtols = [1e-8] if args.quick else [1e-8, 1e-10]
    t_max = 20.0 if args.quick else 70.0
    jobs = [dict(N=N, mode=m, start=s, seed=sd, rtol=rt, t_max=t_max)
            for N in LEVELS for m in MODES for s in ("common", "self_dual") for sd in seeds for rt in rtols]
    print(f"{len(jobs)} runs on {args.procs} processes", flush=True)
    with Pool(args.procs) as pool:
        rows = []
        for r in pool.imap_unordered(run_one, jobs):
            rows.append(r)
            print(f"  N={r['N']:>2} {r['mode']:5} {r['start']:9} seed={r['seed']} rtol={r['rtol']:.0e} "
                  f"ok={r['success']} {r['seconds']}s w0={r.get('w0')} tau={r.get('tau_final')}", flush=True)

    os.makedirs(args.out, exist_ok=True)
    result = dict(
        spec="specs/LEANFLOW_ARCHITECTURE.md L5",
        held_fixed=["V normalisation (modular function in [0,1])", "log-mode scale=10", "RHO_M0", "RHO_R0",
                    "a_initial=1e-10", f"t_max={t_max}"],
        runs=sorted(rows, key=lambda r: (r["mode"], r["start"], r["seed"], r["rtol"], r["N"])),
        verdicts=[] if args.quick else verdicts(rows),
    )
    path = os.path.join(args.out, "quick.json" if args.quick else "results.json")
    with open(path, "w") as fh:
        json.dump(result, fh, indent=1, default=str)
    print("written", path)
    for v in result["verdicts"]:
        print(f"{v['mode']:5} {v['start']:9} {v['observable']:15} meets_L5={v['meets_L5']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
