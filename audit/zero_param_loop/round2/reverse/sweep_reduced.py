#!/usr/bin/env python3
"""
REVERSE LOOP (round 2), step 2a: sweep the 2-parameter reduced model.
Seed 43 (task-specific override of the generic ground-rule seed 42 --
"compute it with GUDHI on a fresh sweep of the reduced model (seed 43)").

Three sweeps are written (see PREREGISTRATION.md item 3 for why the
positive controls exist):
  - sweep_reduced_2d.csv      : joint (mu_sym, c4_pta_product), N=128 Sobol
  - sweep_reduced_mu_only.csv : mu_sym swept, c4_pta_product fixed nominal
  - sweep_reduced_c4_only.csv : c4_pta_product swept, mu_sym fixed nominal

mu_sym range is restricted to [0.2, 3.0] (log10 in [-0.7, 0.48]) to stay
inside the symmetron BVP's numerically-stable region (mu_sym>=~10 is known
to diverge to NaN, per round-1/round-2 forward findings); c4_pta_product
spans the same 8 decades round 1/2 forward used, log10 in [-5, 3].
"""
import csv
import os
import sys

import numpy as np
from scipy.stats import qmc

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reduced_model as rm  # noqa: E402

SEED = 43
N_2D = 128     # 2**7
N_1D = 48

MU_LOG10_BOUNDS = (-0.7, 0.48)         # mu_sym in [0.2, 3.0]
C4_LOG10_BOUNDS = (-5.0, 3.0)          # c4_pta_product in [1e-5, 1e3]

FIELDNAMES = [
    "idx", "sweep", "mu_sym", "c4_pta_product",
    "screening_suppression_factor", "phi_center_ratio", "numerically_stable",
    "pta_max_deviation_from_hd", "H0_implied_km_s_Mpc", "Omega_m_frozen",
    "wall_time_seconds", "error",
]


def eval_row(idx, sweep, mu_sym, c4_pta_product):
    import time
    t0 = time.time()
    row = {"idx": idx, "sweep": sweep, "mu_sym": mu_sym, "c4_pta_product": c4_pta_product, "error": ""}
    try:
        res = rm.evaluate_reduced_point({"mu_sym": mu_sym, "c4_pta_product": c4_pta_product})
        row["screening_suppression_factor"] = res["screening"]["screening_suppression_factor"]
        row["phi_center_ratio"] = res["screening"]["phi_center_ratio"]
        row["numerically_stable"] = res["screening"]["numerically_stable"]
        row["pta_max_deviation_from_hd"] = res["pta"]["max_deviation_from_hd"]
        row["H0_implied_km_s_Mpc"] = res["dark_energy"]["H0_implied_km_s_Mpc"]
        row["Omega_m_frozen"] = res["dark_energy"]["Omega_m_frozen"]
    except Exception as exc:
        row["error"] = str(exc)
        row["screening_suppression_factor"] = float("nan")
        row["phi_center_ratio"] = float("nan")
        row["numerically_stable"] = False
        row["pta_max_deviation_from_hd"] = float("nan")
        row["H0_implied_km_s_Mpc"] = float("nan")
        row["Omega_m_frozen"] = float("nan")
    row["wall_time_seconds"] = time.time() - t0
    return row


def run():
    rows = []

    # nominal point in every sweep, idx 0
    nominal = rm.REDUCED_DEFAULT_PARAMS
    rows.append(eval_row(0, "2d", nominal["mu_sym"], nominal["c4_pta_product"]))

    sampler2d = qmc.Sobol(d=2, scramble=True, seed=SEED)
    u2d = sampler2d.random(N_2D)
    mu2d = 10 ** (MU_LOG10_BOUNDS[0] + u2d[:, 0] * (MU_LOG10_BOUNDS[1] - MU_LOG10_BOUNDS[0]))
    c42d = 10 ** (C4_LOG10_BOUNDS[0] + u2d[:, 1] * (C4_LOG10_BOUNDS[1] - C4_LOG10_BOUNDS[0]))
    for i in range(N_2D):
        rows.append(eval_row(i + 1, "2d", float(mu2d[i]), float(c42d[i])))

    sampler_mu = qmc.Sobol(d=1, scramble=True, seed=SEED)
    u_mu = sampler_mu.random(N_1D)[:, 0]
    mu1d = 10 ** (MU_LOG10_BOUNDS[0] + u_mu * (MU_LOG10_BOUNDS[1] - MU_LOG10_BOUNDS[0]))
    for i in range(N_1D):
        rows.append(eval_row(i, "mu_only", float(mu1d[i]), nominal["c4_pta_product"]))

    sampler_c4 = qmc.Sobol(d=1, scramble=True, seed=SEED)
    u_c4 = sampler_c4.random(N_1D)[:, 0]
    c41d = 10 ** (C4_LOG10_BOUNDS[0] + u_c4 * (C4_LOG10_BOUNDS[1] - C4_LOG10_BOUNDS[0]))
    for i in range(N_1D):
        rows.append(eval_row(i, "c4_only", nominal["mu_sym"], float(c41d[i])))

    out_path = os.path.join(HERE, "sweep_reduced.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    n_2d = sum(1 for r in rows if r["sweep"] == "2d")
    n_2d_stable = sum(1 for r in rows if r["sweep"] == "2d" and r["numerically_stable"])
    n_mu = sum(1 for r in rows if r["sweep"] == "mu_only")
    n_mu_stable = sum(1 for r in rows if r["sweep"] == "mu_only" and r["numerically_stable"])
    n_c4 = sum(1 for r in rows if r["sweep"] == "c4_only")
    n_c4_stable = sum(1 for r in rows if r["sweep"] == "c4_only" and r["numerically_stable"])
    print(f"Wrote {out_path}: {len(rows)} rows")
    print(f"  2d sweep: {n_2d_stable}/{n_2d} numerically stable")
    print(f"  mu_only sweep: {n_mu_stable}/{n_mu} numerically stable")
    print(f"  c4_only sweep: {n_c4_stable}/{n_c4} numerically stable")


if __name__ == "__main__":
    run()
