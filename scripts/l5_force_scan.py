#!/usr/bin/env python3
"""L5, cheap route: does the level-N force (not a full ODE run) differ, beyond tolerance,
between N=7 and N=12? No integrator, no stiff-kink risk -- this is a static grid scan of
dV/dy over the opt-in Gamma_0(N)+ potential (leanflow.core.modular_potential), the same
function scripts/l5_level_experiment.py tried and failed to integrate through.

This is NOT the L5 acceptance criterion as spec'd (specs/LEANFLOW_ARCHITECTURE.md L5 asks for
two completed *runs*). It answers a narrower, cheaper question first: is there a level-dependent
force anywhere along the y-axis, robust to the choices that shouldn't matter (which x, which
numerical precision, which finite-difference step)? If not, there is no point debugging the
ODE integrator to look for it. If so, that is evidence -- not proof -- worth the further
integrator work, and it tells you where in y to aim it (away from the x-inert plateau found in
audit/l5_level_runs/README.md, and away from the stiff kink itself).

Adopts nothing: workshopcosmo.py's hand-built double well is untouched.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from leanflow.core.modular_potential import modular_potential  # noqa: E402

LEVELS = (7, 12)
MODES = ("ratio", "log")
X_PROBES = (0.03, 0.11, 0.19, 0.27, 0.34)  # stand in for "seed": arbitrary generic offsets


def dV_dy(x: float, y: float, N: int, mode: str, h: float, dps: int) -> float:
    lo = modular_potential(complex(x, y - h), N=N, mode=mode, dps=dps)
    hi = modular_potential(complex(x, y + h), N=N, mode=mode, dps=dps)
    return (hi - lo) / (2 * h)


def tolerance(x: float, y: float, N: int, mode: str) -> float:
    """Not assumed: the spread dV/dy takes under two innocuous numerical choices --
    finite-difference step (1e-6 vs 1e-5) and working precision (dps 25 vs 30) -- each
    held at a reference value while the other varies. The larger spread is the tolerance;
    anything smaller than it cannot be told from numerical noise.
    """
    ref = dV_dy(x, y, N, mode, h=1e-6, dps=25)
    spread_h = abs(dV_dy(x, y, N, mode, h=1e-5, dps=25) - ref)
    spread_dps = abs(dV_dy(x, y, N, mode, h=1e-6, dps=30) - ref)
    return max(spread_h, spread_dps)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--y-min", type=float, default=0.20)
    ap.add_argument("--y-max", type=float, default=2.0)
    ap.add_argument("--n-y", type=int, default=24)
    ap.add_argument("--out", default=os.path.join(ROOT, "audit", "l5_level_runs", "force_scan.json"))
    args = ap.parse_args()

    ys = [args.y_min * (args.y_max / args.y_min) ** (i / (args.n_y - 1)) for i in range(args.n_y)]
    rows, verdicts = [], []
    for mode in MODES:
        for y in ys:
            per_x = []
            for x in X_PROBES:
                d7 = dV_dy(x, y, 7, mode, h=1e-6, dps=25)
                d12 = dV_dy(x, y, 12, mode, h=1e-6, dps=25)
                tol = max(tolerance(x, y, 7, mode), tolerance(x, y, 12, mode))
                diff = d12 - d7
                per_x.append(dict(x=x, dVdy_N7=d7, dVdy_N12=d12, diff=diff, tolerance=tol,
                                  exceeds=abs(diff) > tol))
            rows.append(dict(mode=mode, y=y, per_x=per_x))
            signs = {math.copysign(1, p["diff"]) for p in per_x if abs(p["diff"]) > p["tolerance"]}
            verdicts.append(dict(
                mode=mode, y=y,
                all_exceed_tolerance=all(p["exceeds"] for p in per_x),
                consistent_sign=len(signs) <= 1,
                mean_diff=sum(p["diff"] for p in per_x) / len(per_x),
                mean_tolerance=sum(p["tolerance"] for p in per_x) / len(per_x),
            ))
            v = verdicts[-1]
            flag = "LEVEL-DEPENDENT" if v["all_exceed_tolerance"] and v["consistent_sign"] else ""
            print(f"{mode:5} y={y:7.4f}  mean_diff={v['mean_diff']:+.3e}  "
                  f"tol={v['mean_tolerance']:.3e}  {flag}", flush=True)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(dict(x_probes=list(X_PROBES), rows=rows, verdicts=verdicts), fh, indent=1)
    print("written", args.out)

    band = [v for v in verdicts if v["all_exceed_tolerance"] and v["consistent_sign"]]
    print(f"\n{len(band)}/{len(verdicts)} (mode, y) cells: level-dependent force, "
          f"consistent sign, above tolerance at every one of {len(X_PROBES)} x-probes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
